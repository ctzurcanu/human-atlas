import {appendAnatomyLabelCallouts} from './anatomy-label-renderer';
import {layoutAnatomyLabels,labelAreaForPanels,type LabelAnchor,type LabelPanel} from './label-layout';
import {displayLaterality} from './laterality';
import {useEffect,useRef} from 'react';
import sensoryTips from './data/dermatome-sensory-tips.json';
import * as T from 'three';
import {GLTFLoader} from 'three/examples/jsm/loaders/GLTFLoader.js';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import type {AnatomySex,SceneState} from './anatomy';
import {dermatomeSelection,territoryPalette,dermatomeLevels,dermatomeTerritories,dermatomeLabel} from './dermatome-colors';
import {dermatomeInkMask} from './dermatome-ink-mask';
import {nativeDermatomeDepth,type NativeDepthManifest} from './dermatome-depth';
import {enabledSections,sectionStack} from './section-stack';
import {sectionPlaneEquation,sectionFraction,type Bounds3} from './section-plane';
import {depthLayerOpacity,depthPosition} from './depth-control';
import type {SceneCapture} from './export-tools';

/** A separate native reference frame: never substitutes this surface for atlas tissue. */
export default function DermatomeSurface({sex,state,selection,onTerritory,onCapture,onProgress,onError,sectionTool,onSectionPosition,onBounds}:{sectionTool?:boolean;onSectionPosition?:(position:number)=>void;onBounds?:(bounds:Bounds3)=>void;sex:AnatomySex;state:SceneState;selection?:string;onTerritory:(id:string,toggle?:boolean)=>void;onCapture:(capture:SceneCapture|null)=>void;onProgress:(n:number)=>void;onError:(message:string)=>void}){
 const host=useRef<HTMLDivElement>(null),latest=useRef(state),selected=useRef(selection),pickCallback=useRef(onTerritory),sectionToolRef=useRef(sectionTool),sectionPositionRef=useRef(onSectionPosition);sectionToolRef.current=sectionTool;sectionPositionRef.current=onSectionPosition;latest.current=state;selected.current=selection;pickCallback.current=onTerritory;
 useEffect(()=>{
  const el=host.current!;let disposed=false,frame=0,model:T.Group|undefined,lastView='',lastReset=-1,dirty=true,lastSelection='',lastTheme='',lastHidden='uninitialized',lastDepth='',lastDermSelection='uninitialized';
  const renderer=new T.WebGLRenderer({antialias:true,alpha:true,preserveDrawingBuffer:true});
  renderer.setPixelRatio(Math.min(devicePixelRatio,2));renderer.outputColorSpace=T.SRGBColorSpace;
  renderer.toneMapping=T.ACESFilmicToneMapping;renderer.toneMappingExposure=1.05;
  renderer.domElement.dataset.dermatomeSurface=sex;
  renderer.domElement.setAttribute('aria-label',`${sex} dermatome color reference. Drag to orbit and scroll to zoom. Hover a color for its dermatome label and click to select.`);
  el.appendChild(renderer.domElement);
  const scene=new T.Scene(),camera=new T.PerspectiveCamera(34,1,.01,20),controls=new OrbitControls(camera,renderer.domElement);
  controls.enableDamping=true;controls.minDistance=.3;controls.maxDistance=8;
  // Wheel/pinch zoom updates inside OrbitControls' input handler. Its next
  // update() can return false, so invalidate on the change event as well.
  const cameraChanged=()=>{dirty=true;};controls.addEventListener('change',cameraChanged);
  scene.add(new T.HemisphereLight(0xffffff,0xffffff,1.2));
  const light=new T.DirectionalLight(0xffffff,1.4);light.position.set(-2,3,4);scene.add(light);
  const rim=new T.DirectionalLight(0xffffff,1.4);rim.position.set(2,2,-3);scene.add(rim);
  const resources=new Set<T.Material>(),textures=new Set<T.Texture>(),geometries=new Set<T.BufferGeometry>();
  const clipPlanes=[new T.Plane(new T.Vector3(0,1,0),10000),new T.Plane(new T.Vector3(0,1,0),10000)];renderer.localClippingEnabled=true;
  let sectionBounds:Bounds3={min:[-.5,-.9,-.5],max:[.5,.9,.5]},lastSections='';
  const kept=(point:T.Vector3)=>clipPlanes.every(plane=>plane.distanceToPoint(point)>=0);
  const visible=(object:T.Object3D):boolean=>object.visible&&(!object.parent||visible(object.parent));
  let depthContext:ReturnType<typeof nativeDermatomeDepth>|undefined;
  let manifest:NativeDepthManifest|undefined,skinHolder:T.Group|undefined,nativeCenter:T.Vector3|undefined;
  const skinAlpha={value:1},baseSkin={value:new T.Color('#c99a7c')},baseAlbedo:{value:T.Texture|null}={value:null};
  renderer.domElement.dataset.dermatomeOpacity='0.8';
  const initializeDepth=()=>{if(manifest&&skinHolder&&nativeCenter&&!depthContext){depthContext=nativeDermatomeDepth(skinHolder,nativeCenter,manifest,()=>{dirty=true;},onError,clipPlanes);lastDepth='';}};
  fetch(`/local-models/dermatomes/${sex}-depth.json`).then(response=>{if(!response.ok)throw new Error(`Depth manifest HTTP ${response.status}`);return response.json();}).then(raw=>{const value=raw as NativeDepthManifest;if(disposed)return;manifest=value;new T.TextureLoader().load(value.nativeBodyAlbedo,texture=>{if(disposed){texture.dispose();return;}texture.flipY=false;texture.colorSpace=T.SRGBColorSpace;baseAlbedo.value=texture;textures.add(texture);dirty=true;});initializeDepth();}).catch(error=>{if(!disposed)onError(String(error));});
  const inkMaps=new Map<T.Texture,T.DataTexture>();
  const activeLevel={value:-1},activeRegion={value:-1},activeSide={value:0};
  const selectedTerritories={value:new Float32Array(64)},hasSelectedTerritories={value:0};
  const levels=dermatomeLevels,territoryVisibility={value:new Float32Array(dermatomeTerritories.length).fill(1)};
  const pixels=new Map<T.Texture,{data:Uint8ClampedArray;width:number;height:number}>(),raycaster=new T.Raycaster();let pointerDown=[0,0];
  const tip=document.createElement('div');tip.className='dermatome-hover-tip';tip.setAttribute('role','tooltip');tip.hidden=true;
  Object.assign(tip.style,{position:'absolute',pointerEvents:'none',zIndex:'8',maxWidth:'260px',padding:'8px 12px',borderRadius:'8px',background:'#202a35',color:'#f5f6f7',border:'1px solid #52616f',fontSize:'13px',whiteSpace:'pre-line'});el.appendChild(tip);
  const anchors=new Map<string,{object:T.Object3D;local:T.Vector3}>();
  const labelLayer=document.createElementNS('http://www.w3.org/2000/svg','svg');labelLayer.classList.add('anatomy-labels');labelLayer.setAttribute('aria-hidden','true');el.appendChild(labelLayer);
  const labelMeasure=document.createElement('canvas');let lastLabelArea='',lastLayoutRead=0;
  let labelArea=labelAreaForPanels(el.clientWidth,el.clientHeight,[]);
  const readLabelArea=()=>{const viewport=el.getBoundingClientRect(),panels:LabelPanel[]=[];
   for(const panel of el.ownerDocument.querySelectorAll('.identity,.top-actions,.advanced-panel,.layers-panel,.study-panel,.section-panel,.detail-sheet,.download-panel,.search-panel,.view-controls,.share-view,.guest-menu,.anatomy-choice-menu,.anatomy-search-results')){
    if(!panel.getClientRects().length||getComputedStyle(panel).visibility==='hidden')continue;const r=panel.getBoundingClientRect();if(!r.width||!r.height)continue;
    panels.push({kind:panel.matches('.identity,.top-actions')?'top':panel.matches('.advanced-panel')?'bottom':'panel',left:r.left-viewport.left,right:r.right-viewport.left,top:r.top-viewport.top,bottom:r.bottom-viewport.top});
   }return labelAreaForPanels(viewport.width,viewport.height,panels);
  };
  const findSurfaceAnchor=(key:string)=>{
   if(!model)return;const [level,side]=key.split(':');let best=-Infinity,result:{object:T.Object3D;local:T.Vector3}|undefined;
   model.traverse(object=>{if(!(object instanceof T.Mesh))return;const material=(Array.isArray(object.material)?object.material[0]:object.material) as T.MeshStandardMaterial,map=material.map,image=map&&pixels.get(map);if(!image)return;
    const positions=object.geometry.getAttribute('position'),uvs=object.geometry.getAttribute('uv'),indices=object.geometry.index;if(!uvs)return;
    const palette=territoryPalette(sex,material.name==='head').map(([name,rgb])=>({name,color:new T.Color().setRGB(...rgb.map(v=>v/255) as [number,number,number],T.SRGBColorSpace)}));
    for(let i=0;i<(indices?.count??positions.count);i+=3){const ids=[0,1,2].map(j=>indices?indices.getX(i+j):i+j),points=ids.map(id=>new T.Vector3().fromBufferAttribute(positions,id)),local=points[0].clone().add(points[1]).add(points[2]).multiplyScalar(1/3);if((local.x>=0?'left':'right')!==side)continue;
     const u=ids.reduce((n,id)=>n+uvs.getX(id),0)/3,v=ids.reduce((n,id)=>n+uvs.getY(id),0)/3,x=Math.max(0,Math.min(image.width-1,Math.floor(u*image.width))),y=Math.max(0,Math.min(image.height-1,Math.floor(v*image.height))),offset=(y*image.width+x)*4;
     if((inkMaps.get(map!)?.image.data?.[y*image.width+x]??0)>127)continue;const rgb=Array.from(image.data.slice(offset,offset+3));if(Math.max(...rgb)<80)continue;const color=new T.Color().setRGB(...rgb.map(v=>v/255) as [number,number,number],T.SRGBColorSpace);const nearest=palette.reduce((a,b)=>color.toArray().reduce((n,c,j)=>n+(c-b.color.toArray()[j])**2,0)<color.toArray().reduce((n,c,j)=>n+(c-a.color.toArray()[j])**2,0)?b:a);if(nearest.name!==level)continue;
     const world=object.localToWorld(local.clone()),normal=points[1].clone().sub(points[0]).cross(points[2].clone().sub(points[0])).transformDirection(object.matrixWorld),facing=normal.dot(camera.position.clone().sub(world).normalize());if(facing<=0)continue;const projected=world.clone().project(camera);if(Math.abs(projected.x)>1||Math.abs(projected.y)>1||Math.abs(projected.z)>1)continue;
     const score=facing-.1*(projected.x**2+projected.y**2);if(score>best){best=score;result={object,local};}
    }
   });if(result)anchors.set(key,result);return result;
  };
  const drawLabels=(s:SceneState)=>{labelLayer.replaceChildren();if(s.labels===false||skinAlpha.value<=0)return;
   const w=el.clientWidth,h=el.clientHeight;labelLayer.setAttribute('viewBox',`0 0 ${w} ${h}`);const projected:LabelAnchor[]=[];
   for(const key of s.dermatomeSelected??[]){
    const anchor=anchors.get(key)??findSurfaceAnchor(key);if(!anchor)continue;const world=anchor.object.localToWorld(anchor.local.clone());if(!kept(world))continue;const point=world.project(camera);
    if(point&&(point.z< -1||point.z>1||Math.abs(point.x)>1||Math.abs(point.y)>1))continue;
    const name=dermatomeLabel(key),display=displayLaterality(name);
    projected.push({id:key,x:(point.x+1)*w/2,y:(1-point.y)*h/2,text:display.label,fullText:name,laterality:display.side});
   }
   const measure=(text:string)=>{const context=labelMeasure.getContext('2d');if(!context)return text.length*6.8;context.font='12px Arial';return context.measureText(text).width;};
   appendAnatomyLabelCallouts(labelLayer,layoutAnatomyLabels(projected,labelArea,measure));
  };
  const hideTip=()=>{tip.hidden=true;renderer.domElement.style.cursor='';};
  const territoryAt=(e:PointerEvent)=>{
   if(!model?.visible)return;
   const rect=renderer.domElement.getBoundingClientRect();raycaster.setFromCamera(new T.Vector2((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1),camera);
   const hit=raycaster.intersectObject(model,true).find(hit=>hit.object.visible&&kept(hit.point));
   if(!hit?.uv||!(hit.object instanceof T.Mesh))return;
   const material=(Array.isArray(hit.object.material)?hit.object.material[0]:hit.object.material) as T.MeshStandardMaterial,map=material.map;
   if(!map)return;
   const image=pixels.get(map);if(!image)return;
   const sample=(u:number,v:number)=>{const x=Math.max(0,Math.min(image.width-1,Math.floor(u*image.width))),y=Math.max(0,Math.min(image.height-1,Math.floor(v*image.height))),at=(y*image.width+x)*4;return {rgb:Array.from(image.data.slice(at,at+3)),at:y*image.width+x};};
   let pigment=sample(hit.uv.x,hit.uv.y);const mask=inkMaps.get(map)?.image.data;
   // Match the rendered annotation cleanup so hovering a removed glyph still identifies its color patch.
   if(mask&&mask[pigment.at]>127){let found=false;for(let radius=0;radius<6&&!found;radius++)for(let direction=0;direction<8&&!found;direction++){const step=2**(radius+1)/2048,angle=direction*Math.PI/4,candidate=sample(hit.uv.x+step*Math.cos(angle),hit.uv.y+step*Math.sin(angle));if(mask[candidate.at]<128&&Math.max(...candidate.rgb)>80){pigment=candidate;found=true;}}}
   if(Math.max(...pigment.rgb)<80)return;
   // The shader compares linear RGB, rather than encoded PNG values.
   const color=new T.Color().setRGB(...pigment.rgb.map(v=>v/255) as [number,number,number],T.SRGBColorSpace);
   const palette=territoryPalette(sex,material.name==='head');
   const distance=(rgb:number[])=>color.toArray().reduce((sum,v,i)=>sum+(v-new T.Color().setRGB(...rgb.map(v=>v/255) as [number,number,number],T.SRGBColorSpace).toArray()[i])**2,0);
   const nearest=palette.reduce((a,b)=>distance(b[1])<distance(a[1])?b:a);
   const local=hit.object.worldToLocal(hit.point.clone()),side=local.x>=0?'left':'right',level=nearest[0];
   if(latest.current.dermatomeHidden?.includes(`${level}:${side}`))return;
   return {level,key:`${level}:${side}`,object:hit.object,local,id:`${level.startsWith('V')?'TRIGEMINAL':'DERMATOME'}:${level}:${side}`,label:`${side==='left'?'Left':'Right'} ${level} ${level.startsWith('V')?'facial sensory territory':'dermatome'}`};
  };
  const down=(e:PointerEvent)=>{pointerDown=[e.clientX,e.clientY];hideTip();};
  const up=(e:PointerEvent)=>{if(Math.hypot(e.clientX-pointerDown[0],e.clientY-pointerDown[1])>5)return;if(sectionToolRef.current&&skinHolder){const rect=renderer.domElement.getBoundingClientRect();raycaster.setFromCamera(new T.Vector2((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1),camera);const hit=raycaster.intersectObject(skinHolder,true).find(hit=>visible(hit.object)&&kept(hit.point));if(hit){const stack=sectionStack(latest.current),section=stack[Math.min(latest.current.activeSection??0,stack.length-1)];sectionPositionRef.current?.(sectionFraction(sectionBounds,section,hit.point.toArray()));}return;}const territory=territoryAt(e);if(territory){anchors.set(territory.key,{object:territory.object,local:territory.local});pickCallback.current(territory.id,e.ctrlKey||e.metaKey);}};
  let lastHover=0;
  const move=(e:PointerEvent)=>{if(e.buttons){hideTip();return;}const now=performance.now();if(now-lastHover<40)return;lastHover=now;const territory=territoryAt(e);if(!territory){hideTip();return;}tip.textContent=`${territory.label}\n${(sensoryTips as Record<string,string>)[territory.level]??'Trigeminal facial sensation.'}\nClick to select · Ctrl/⌘-click to add or remove`;tip.hidden=false;const rect=el.getBoundingClientRect();tip.style.left=`${Math.max(4,Math.min(e.clientX-rect.left+14,el.clientWidth-tip.offsetWidth-8))}px`;tip.style.top=`${Math.max(4,Math.min(e.clientY-rect.top+14,el.clientHeight-tip.offsetHeight-8))}px`;renderer.domElement.style.cursor='pointer';};
  renderer.domElement.addEventListener('pointerdown',down);renderer.domElement.addEventListener('pointerup',up);renderer.domElement.addEventListener('pointermove',move);renderer.domElement.addEventListener('pointerleave',hideTip);
  const resize=()=>{const w=Math.max(1,el.clientWidth),h=Math.max(1,el.clientHeight);camera.aspect=w/h;camera.updateProjectionMatrix();renderer.setSize(w,h);dirty=true;};
  const observer=new ResizeObserver(resize);observer.observe(el);resize();
  onProgress(0);
  new GLTFLoader().load(`/local-models/dermatomes/${sex}-surface.gltf`,gltf=>{
   if(disposed){gltf.scene.traverse(o=>{if(o instanceof T.Mesh){o.geometry.dispose();for(const m of Array.isArray(o.material)?o.material:[o.material]){(m as T.MeshStandardMaterial).map?.dispose();m.dispose();}}});return;}
   model=gltf.scene;
   model.traverse(object=>{if(!(object instanceof T.Mesh))return;
    geometries.add(object.geometry);
    for(const original of Array.isArray(object.material)?object.material:[object.material]){
     const material=original as T.MeshStandardMaterial;resources.add(material);if(material.map){textures.add(material.map);material.map.generateMipmaps=false;material.map.minFilter=T.NearestFilter;material.map.magFilter=T.NearestFilter;material.map.needsUpdate=true;}
     if(material.map&&!inkMaps.has(material.map)){const image=material.map.image as HTMLImageElement,canvas=document.createElement('canvas');canvas.width=image.width;canvas.height=image.height;const ctx=canvas.getContext('2d',{willReadFrequently:true})!;ctx.drawImage(image,0,0);const data=ctx.getImageData(0,0,canvas.width,canvas.height).data;pixels.set(material.map,{data,width:canvas.width,height:canvas.height});const {mask,annotations}=dermatomeInkMask(data,canvas.width,canvas.height),ink=new T.DataTexture(mask,canvas.width,canvas.height,T.RedFormat);ink.flipY=false;ink.minFilter=T.NearestFilter;ink.magFilter=T.NearestFilter;ink.needsUpdate=true;inkMaps.set(material.map,ink);textures.add(ink);renderer.domElement.dataset.dermatomeInkComponents=String(annotations);}
     material.transparent=false;material.opacity=1;material.alphaTest=0;material.depthWrite=true;material.blending=T.NormalBlending;
     material.clippingPlanes=clipPlanes;material.side=T.DoubleSide;material.roughness=.85;material.metalness=0;
     material.customProgramCacheKey=()=> `dermatome-colors-skin-overlay-v5:${sex}:${material.name}`;
     if(!material.map){material.needsUpdate=true;continue;}
     material.onBeforeCompile=shader=>{
      const palette=territoryPalette(sex,material.name==='head');
      shader.uniforms.annotationInk={value:inkMaps.get(material.map!)};shader.uniforms.activeLevel=activeLevel;shader.uniforms.activeRegion=activeRegion;shader.uniforms.activeSide=activeSide;shader.uniforms.territoryVisibility=territoryVisibility;shader.uniforms.selectedTerritories=selectedTerritories;shader.uniforms.hasSelectedTerritories=hasSelectedTerritories;shader.uniforms.skinAlpha=skinAlpha;shader.uniforms.baseSkin=baseSkin;shader.uniforms.baseAlbedo=baseAlbedo;
      shader.vertexShader='varying float territorySide;\n'+shader.vertexShader;
      shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nterritorySide=position.x;');
      shader.fragmentShader='uniform float selectedTerritories[64];uniform float hasSelectedTerritories;uniform float skinAlpha;uniform vec3 baseSkin;uniform sampler2D baseAlbedo;uniform float territoryVisibility[64];uniform sampler2D annotationInk;uniform float activeLevel;uniform float activeRegion;uniform float activeSide;varying float territorySide;\n'+shader.fragmentShader;
      // Sample the supplied UV territory colors directly. Dark baked glyphs
      // and leaders use the nearest surrounding pigment; no generated map,
      // body-space boundary or new skin geometry is substituted.
      shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',`
       #ifdef USE_MAP
       vec4 territory=texture2D(map,vMapUv);
       if(texture2D(annotationInk,vMapUv).r>0.5){
        bool found=false;
        for(int radius=0;radius<6;radius++){
         float stepSize=pow(2.0,float(radius)+1.0)/2048.0;
         for(int direction=0;direction<8;direction++){
          float angle=float(direction)*0.78539816339;
          vec4 pigment=texture2D(map,vMapUv+stepSize*vec2(cos(angle),sin(angle)));
          vec2 candidateUv=vMapUv+stepSize*vec2(cos(angle),sin(angle));
          bool clean=texture2D(annotationInk,candidateUv).r<0.5 && max(pigment.r,max(pigment.g,pigment.b))>0.08;
          if(!found&&clean){territory=pigment;found=true;}
         }
        }
       }
       float nearestDistance=100.0;float territoryId=-2.0;float territoryRegion=-2.0;
       ${palette.map(([level,rgb])=>{const c=new T.Color(`rgb(${rgb.join(',')})`);return `{float distance=dot(territory.rgb-vec3(${c.r},${c.g},${c.b}),territory.rgb-vec3(${c.r},${c.g},${c.b}));if(distance<nearestDistance){nearestDistance=distance;territoryId=${levels.indexOf(level)}.0;territoryRegion=${['C','T','L','S','V'].indexOf(level[0])}.0;}}`;}).join('\n')}
       float overlayVisibility=territoryId>=0.0?territoryVisibility[int(territoryId)*2+(territorySide>=0.0?0:1)]:1.0;
       bool territorySelected=hasSelectedTerritories>0.5?(territoryId>=0.0&&selectedTerritories[int(territoryId)*2+(territorySide>=0.0?0:1)]>0.5):(activeLevel<0.0 || abs(territoryId-activeLevel)<0.1) && (activeRegion<0.0 || abs(territoryRegion-activeRegion)<0.1) && (activeSide==0.0 || territorySide*activeSide>=0.0);
       if(!territorySelected){float gray=dot(territory.rgb,vec3(.2126,.7152,.0722));territory.rgb=mix(vec3(gray*.3),territory.rgb,.12);}
       // Territory pigments and selection affect RGB only; skin remains beneath the color overlay.
       vec3 skinPigment=${material.name==='body'?'texture2D(baseAlbedo,vMapUv).rgb':'baseSkin'};
       diffuseColor.rgb*=mix(skinPigment,territory.rgb,0.8*overlayVisibility);diffuseColor.a=skinAlpha;
       if(skinAlpha<=0.0)discard;
       #endif`);
     };material.needsUpdate=true;
    }
   });
   // Normalize only this independent reference view, uniformly. The atlas
   // anatomy and its saved geometry/frame are untouched.
   const box=new T.Box3().setFromObject(model),size=box.getSize(new T.Vector3()),center=box.getCenter(new T.Vector3());
   model.position.sub(center);const holder=new T.Group();holder.add(model);holder.scale.setScalar(1.8/size.y);scene.add(holder);skinHolder=holder;nativeCenter=center;holder.updateMatrixWorld(true);const bounds=new T.Box3().setFromObject(holder);sectionBounds={min:bounds.min.toArray(),max:bounds.max.toArray()};onBounds?.(sectionBounds);lastSections='';initializeDepth();
   camera.position.set(0,.1,3.8);controls.target.set(0,0,0);controls.update();lastView='';lastReset=-1;onProgress(100);
   renderer.domElement.dataset.dermatomeNodes=String(geometries.size);lastHidden='uninitialized',lastDepth='';dirty=true;
  },undefined,error=>{if(!disposed)onError(`Could not load the local dermatome surface: ${error instanceof Error?error.message:String(error)}`);});
  const draw=()=>{if(disposed)return;frame=requestAnimationFrame(draw);const s=latest.current;
   const sectionKey=JSON.stringify([s.sections,s.section]);if(sectionKey!==lastSections){lastSections=sectionKey;clipPlanes.forEach(plane=>{plane.normal.set(0,1,0);plane.constant=10000;});for(const {section,index} of enabledSections(s)){const equation=sectionPlaneEquation(sectionBounds,section);clipPlanes[index].normal.fromArray(equation.normal);clipPlanes[index].constant=equation.constant;}dirty=true;}
   const dermKey=JSON.stringify([s.dermatomeSelected,s.labels]);if(dermKey!==lastDermSelection){lastDermSelection=dermKey;const keys=new Set(s.dermatomeSelected??[]);selectedTerritories.value.set(dermatomeTerritories.map(key=>keys.has(key)?1:0));hasSelectedTerritories.value=keys.size?1:0;
    dirty=true;
   }
   const hiddenKey=(s.dermatomeHidden??[]).join(',');if(hiddenKey!==lastHidden){lastHidden=hiddenKey;const hidden=new Set(s.dermatomeHidden);territoryVisibility.value.set(dermatomeTerritories.map(id=>hidden.has(id)?0:1));dirty=true;}
   const depthKey=JSON.stringify([s.depth,s.depthHidden,s.skinOpacity]);if(depthKey!==lastDepth){lastDepth=depthKey;skinAlpha.value=s.depthHidden?.includes('skin')?0:depthLayerOpacity('skin',s);if(model)model.visible=skinAlpha.value>0;for(const original of resources){const material=original as T.MeshStandardMaterial;const transparent=skinAlpha.value<.999;if(material.transparent!==transparent){material.transparent=transparent;material.needsUpdate=true;}material.depthWrite=!transparent;material.opacity=material.map?1:skinAlpha.value;}depthContext?.update(s);renderer.domElement.dataset.dermatomeDepth=String(depthPosition(s));dirty=true;}
   if(lastSelection!==selected.current){lastSelection=selected.current??'';dirty=true;}
   const territory=dermatomeSelection(selected.current);activeLevel.value=territory.level?levels.indexOf(territory.level):-1;
   if(territory.level&&activeLevel.value<0)activeLevel.value=99;activeRegion.value=territory.region?['C','T','L','S','V'].indexOf(territory.region):-1;if(territory.region&&activeRegion.value<0)activeRegion.value=99;activeSide.value=territory.side==='left'?1:territory.side==='right'?-1:0;
   if(s.view!==lastView||s.reset!==lastReset){if(s.reset!==lastReset&&skinHolder)skinHolder.rotation.y=0;const distance=camera.position.distanceTo(controls.target)||3.8;
    const direction=s.view==='back'?[0,0,-1]:s.view==='side'?[1,0,0]:s.view==='right'?[-1,0,0]:s.view==='superior'?[0,1,.001]:s.view==='inferior'?[0,-1,.001]:s.view==='three-quarter'?[.4,.04,1]:[0,0,1];
    camera.position.copy(controls.target).add(new T.Vector3().fromArray(direction).normalize().multiplyScalar(distance));lastView=s.view;lastReset=s.reset;dirty=true;
   }
   if(s.rotate&&skinHolder){skinHolder.rotation.y+=.003;dirty=true;}
   const theme=document.documentElement.dataset.theme??'';if(theme!==lastTheme){lastTheme=theme;renderer.setClearColor(theme==='light'?'#f2f3f3':'#141b23');dirty=true;}
   if(controls.update())dirty=true;if(performance.now()-lastLayoutRead>150){lastLayoutRead=performance.now();const area=readLabelArea(),key=JSON.stringify(area);if(key!==lastLabelArea){labelArea=area;lastLabelArea=key;dirty=true;}}if(dirty){scene.updateMatrixWorld(true);drawLabels(s);renderer.render(scene,camera);renderer.domElement.dataset.nativeDepthLayers=JSON.stringify(depthContext?.stats()??{layers:[],geometryBytes:0,drawMeshes:0});renderer.domElement.dataset.renderCalls=String(renderer.info.render.calls);dirty=false;}
  };draw();
  onCapture(async options=>{const alpha=renderer.getClearAlpha();if(options?.background===false)renderer.setClearAlpha(0);renderer.render(scene,camera);const blob=await new Promise<Blob>((resolve,reject)=>renderer.domElement.toBlob(value=>value?resolve(value):reject(new Error('Could not capture the dermatome surface.')),'image/png'));renderer.setClearAlpha(alpha);dirty=true;return blob;});
  return()=>{disposed=true;depthContext?.dispose();onCapture(null);cancelAnimationFrame(frame);observer.disconnect();controls.removeEventListener('change',cameraChanged);controls.dispose();renderer.domElement.removeEventListener('pointerdown',down);renderer.domElement.removeEventListener('pointerup',up);renderer.domElement.removeEventListener('pointermove',move);renderer.domElement.removeEventListener('pointerleave',hideTip);tip.remove();labelLayer.remove();for(const g of geometries)g.dispose();for(const m of resources)m.dispose();for(const t of textures)t.dispose();renderer.dispose();renderer.domElement.remove();};
 },[sex]);
 return <div ref={host} className="anatomy-scene dermatome-reference" style={{position:'absolute',inset:0}}/>;
}
