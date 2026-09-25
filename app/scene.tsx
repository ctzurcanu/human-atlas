import {assetUrl} from './asset-url';
import {partVisible} from './viewer-state';
import {useEffect,useRef} from 'react';
import * as T from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import {RoomEnvironment} from 'three/examples/jsm/environments/RoomEnvironment.js';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';
import {createHierarchicalExplosionLayout,type ExplodeHierarchy,type HierarchicalExplosionLayout} from './hierarchical-explosion';
import {hierarchyEntries} from './anatomy-hierarchy';
import {decodeModelResponse} from './model-download';
import {selectionCenter} from './camera-pivot';
import {PointerTap} from './pointer-tap';
import {SYSTEMS,isSurfaceSystem,structureName,type Atlas,type SceneState} from './anatomy';
interface Props {atlas:Atlas;state:SceneState;hierarchy:ExplodeHierarchy;onSelect:(id:string,add?:boolean)=>void;onCamera?:(camera:number[])=>void;onCovering?:(ids:string[])=>void;onExplosionSteps?:(steps:number)=>void;allowFadedSelection?:boolean;onCapture?:(capture:(()=>Promise<Blob>)|null)=>void;onProgress:(n:number)=>void;onError:(s:string)=>void}
export default function AnatomyScene({atlas,state,hierarchy,onSelect,onCamera,onCovering,onExplosionSteps,allowFadedSelection=false,onCapture,onProgress,onError}:Props){
 const host=useRef<HTMLDivElement>(null),latest=useRef(state),hierarchyRef=useRef(hierarchy),select=useRef(onSelect),cameraCallback=useRef(onCamera),coveringCallback=useRef(onCovering),stepsCallback=useRef(onExplosionSteps),selectFaded=useRef(allowFadedSelection);
 latest.current=state;hierarchyRef.current=hierarchy;select.current=onSelect;cameraCallback.current=onCamera;coveringCallback.current=onCovering;stepsCallback.current=onExplosionSteps;selectFaded.current=allowFadedSelection;
 useEffect(()=>{
  const el=host.current!;let disposed=false,frame=0,dirty=true,ready=false,lastView='',lastReset=-1,lastIsolate='',lastSelection='',lastRegion='',layoutKey='',amount=0;
  let lastState:SceneState|null=null;
  const abort=new AbortController();
  let renderer:T.WebGLRenderer;
  try{renderer=new T.WebGLRenderer({antialias:true,alpha:false,powerPreference:'high-performance'});}catch{onError('This browser could not start the 3D viewer. Please try a browser with WebGL enabled.');return;}
  renderer.setPixelRatio(Math.min(devicePixelRatio,innerWidth<768?1.5:2));renderer.setClearColor('#f2f3f3');renderer.outputColorSpace=T.SRGBColorSpace;renderer.toneMapping=T.ACESFilmicToneMapping;renderer.toneMappingExposure=1.12;el.appendChild(renderer.domElement);
  renderer.domElement.setAttribute('aria-label',atlas.scope==='cell'?'Interactive human cell. Drag to orbit, pinch or scroll to zoom, and tap a component to inspect it.':'Interactive human anatomy. Drag to orbit, pinch or scroll to zoom, and tap a structure to inspect it.');
  const scene=new T.Scene(),camera=new T.PerspectiveCamera(34,1,.01,100),controls=new OrbitControls(camera,renderer.domElement);
  camera.position.set(1.4,1.05,3.6);controls.target.set(0,.85,0);controls.zoomToCursor=true;controls.enableDamping=true;controls.dampingFactor=.085;controls.minDistance=.003;controls.maxDistance=6;controls.zoomSpeed=1.25;controls.maxPolarAngle=Math.PI-.001;controls.addEventListener('change',()=>{dirty=true;});
  const pmrem=new T.PMREMGenerator(renderer),room=new RoomEnvironment(),env=pmrem.fromScene(room,.04);scene.environment=env.texture;room.dispose();pmrem.dispose();
  scene.add(new T.HemisphereLight(0xffffff,0xa7acb2,1.05));
  const key=new T.DirectionalLight(0xfffaf4,2.3);key.position.set(-2,4,3);scene.add(key);
  const rim=new T.DirectionalLight(0xe9f0ff,1.8);rim.position.set(2,2,-3);scene.add(rim);
  const theme=window.matchMedia('(prefers-color-scheme: dark)');
  const applyTheme=()=>{renderer.setClearColor(theme.matches?'#141b23':'#f2f3f3');dirty=true;};
  applyTheme();theme.addEventListener('change',applyTheme);
  const width=T.MathUtils.ceilPowerOfTwo(atlas.parts.length),data=new Float32Array(width*4),partTexture=new T.DataTexture(data,width,1,T.RGBAFormat,T.FloatType);partTexture.needsUpdate=true;
  const selectedData=new Uint8Array(width*4),selectionTexture=new T.DataTexture(selectedData,width,1);selectionTexture.needsUpdate=true;
  const rotationData=new Uint8Array(width*4),rotationTexture=new T.DataTexture(rotationData,width,1);rotationTexture.needsUpdate=true;
  const selectedRotation=new T.Quaternion(),rotationUniform={value:new T.Vector4(0,0,0,1)},pivotUniform={value:new T.Vector3()},rotationPartIndices=new Set<number>(),partIndices=new Map(atlas.parts.map((part,i)=>[part.id,i]));
  const materials:T.Material[]=[],geometries:T.BufferGeometry[]=[],pickers:(T.Mesh|undefined)[]=[],centers=atlas.parts.map(p=>new T.Vector3().fromArray(p.bounds[0]).add(new T.Vector3().fromArray(p.bounds[1])).multiplyScalar(.5));
  const bounds=atlas.parts.map(p=>new T.Box3(new T.Vector3().fromArray(p.bounds[0]),new T.Vector3().fromArray(p.bounds[1])));
  const anatomyEntries=hierarchyEntries(atlas);let explosionLayout:HierarchicalExplosionLayout|null=null;
  let packingWidth=1,packingHeight=1;
  const hover=document.createElement('div');hover.className='part-hover';hover.setAttribute('role','tooltip');hover.hidden=true;el.appendChild(hover);
  type Target={index:number;x:number;y:number;left:number;right:number;top:number;bottom:number};let targets:Target[]=[];
  const projected=new T.Vector3();
  const findTarget=(x:number,y:number,radius:number,eligible:(index:number)=>boolean=()=>true)=>{
   let best=-1,score=Infinity;
   for(const t of targets){if(!eligible(t.index))continue;const dx=Math.max(t.left-x,0,x-t.right),dy=Math.max(t.top-y,0,y-t.bottom),distance=Math.hypot(dx,dy);if(distance>radius)continue;const candidate=distance+Math.hypot(t.x-x,t.y-y)*.025;if(candidate<score){score=candidate;best=t.index;}}
   return best;
  };
  const clipPlane=new T.Plane(new T.Vector3(0,1,0),10000);renderer.localClippingEnabled=true;
  const contextUniform={value:1},skinUniform={value:.1};
  const textures=new Map<string,T.Texture>(),textureLoader=new T.TextureLoader();
  const textureFor=(url:string)=>{
   let texture=textures.get(url);
   if(!texture){
    texture=textureLoader.load(assetUrl(url),()=>{if(!disposed)dirty=true;},undefined,()=>{if(!disposed)onError(`Could not load local anatomy texture: ${url}`);});
    texture.colorSpace=T.SRGBColorSpace;texture.wrapS=T.RepeatWrapping;texture.wrapT=T.RepeatWrapping;
    texture.anisotropy=Math.min(4,renderer.capabilities.getMaxAnisotropy());
    textures.set(url,texture);
   }
   return texture;
  };
  const materialFor=(system:string,ghost=false,source?:{color:number[];map?:string;opacity?:number;vertexColors?:boolean})=>{
   const color=source?.color?`#${source.color.map(value=>Math.round(Math.max(0,Math.min(255,value))).toString(16).padStart(2,'0')).join('')}`:system==='tendon'?'#d6c8b0':system==='cartilage'?'#afc0cb':SYSTEMS.find(s=>s.id===system)?.color??'#aebbb8';
   const m=new T.MeshStandardMaterial({color,map:source?.map?textureFor(source.map):null,vertexColors:!!source?.vertexColors,metalness:.03,roughness:.6,side:T.DoubleSide,transparent:ghost,opacity:1,depthWrite:!ghost,clippingPlanes:[clipPlane]});
   m.customProgramCacheKey=()=>system+':'+ghost;
   m.onBeforeCompile=shader=>{
    shader.uniforms.contextOpacity=contextUniform;shader.uniforms.skinOpacity=skinUniform;shader.uniforms.partState={value:partTexture};shader.uniforms.selectionState={value:selectionTexture};shader.uniforms.rotationState={value:rotationTexture};shader.uniforms.stateWidth={value:width};shader.uniforms.selectionRotation=rotationUniform;shader.uniforms.selectionPivot=pivotUniform;
    shader.vertexShader='attribute float partIndex; uniform sampler2D partState; uniform sampler2D selectionState; uniform sampler2D rotationState; uniform float stateWidth; uniform vec4 selectionRotation; uniform vec3 selectionPivot; varying float partVisible; varying float partSelected; vec3 rotateSelected(vec3 v){return v+2.0*cross(selectionRotation.xyz,cross(selectionRotation.xyz,v)+selectionRotation.w*v);}\n'+shader.vertexShader;
    shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nif (texture2D(rotationState,vec2((partIndex+0.5)/stateWidth,0.5)).r>0.5) objectNormal=rotateSelected(objectNormal);');
    shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvec2 stateUv = vec2((partIndex + 0.5) / stateWidth, 0.5); vec4 state = texture2D(partState, stateUv); transformed += state.xyz; partVisible = state.w; partSelected = texture2D(selectionState, stateUv).r; if(texture2D(rotationState,stateUv).r>0.5) transformed=selectionPivot+rotateSelected(transformed-selectionPivot);');
    shader.fragmentShader='uniform float contextOpacity; uniform float skinOpacity; varying float partVisible; varying float partSelected;\n'+shader.fragmentShader;
    shader.fragmentShader=shader.fragmentShader.replace('#include <clipping_planes_fragment>','#include <clipping_planes_fragment>\nif (partVisible < 0.5) discard;');
    shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>',`#include <color_fragment>
     float tissueAlpha = partSelected > 0.5 ? 1.0 : contextOpacity * ${system==='integumentary'||system==='cell-boundary'?'skinOpacity':'1.0'} * ${(source?.opacity??1).toFixed(6)};
     if (${ghost?'tissueAlpha >= 0.999 || tissueAlpha < 0.001':'tissueAlpha < 0.999'}) discard;
     diffuseColor.a = tissueAlpha;`);
    shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>','#include <map_fragment>\ndiffuseColor.rgb = mix(diffuseColor.rgb, vec3(0.42, 0.85, 0.78), partSelected * 0.35);');
   };materials.push(m);return m;
  };
  const materialSpecs=new Map<string,{system:string;source?:{color:number[];map?:string;opacity?:number;vertexColors?:boolean}}>([...SYSTEMS.map(s=>s.id),'tendon','cartilage'].map(id=>[id,{system:id}]));
  for(const part of atlas.parts)if(part.material&&atlas.materials?.[part.material])materialSpecs.set(`source:${part.system}:${part.material}`,{system:part.system,source:atlas.materials[part.material]});
  const materialKeyFor=(part:typeof atlas.parts[number])=>part.material&&atlas.materials?.[part.material]?`source:${part.system}:${part.material}`:part.tissue==='tendon'||part.tissue==='cartilage'?part.tissue:part.system;
  const mats=new Map([...materialSpecs].map(([key,{system,source}])=>[key,materialFor(system,false,source)])),ghostMats=new Map([...materialSpecs].map(([key,{system,source}])=>[key,materialFor(system,true,source)]));
  const labelLayer=document.createElementNS('http://www.w3.org/2000/svg','svg');labelLayer.classList.add('anatomy-labels');labelLayer.setAttribute('aria-hidden','true');el.appendChild(labelLayer);
  const labelAnchors=new Map<number,T.Vector3>();
  const cameraValues=()=>[...camera.position.toArray(),...controls.target.toArray(),camera.view?.enabled?camera.view.offsetX/camera.view.fullWidth:0,camera.view?.enabled?camera.view.offsetY/camera.view.fullHeight:0];
  controls.addEventListener('end',()=>cameraCallback.current?.(cameraValues()));
  let loaded=0;
  const loadChunk=async(ci:number)=>{
   const chunk=atlas.chunks[ci],compressed=!!chunk.gzip&&typeof DecompressionStream!=='undefined';const response=await fetch(assetUrl(compressed?chunk.gzip!:chunk.url),{signal:abort.signal});const buffer=await decodeModelResponse(response,chunk.bytes,compressed);if(disposed)return;
   const groups=new Map<string,T.BufferGeometry[]>();
   atlas.parts.forEach((p,i)=>{
    if(p.chunk!==ci)return;
    const g=new T.BufferGeometry();g.setAttribute('position',new T.BufferAttribute(new Float32Array(buffer,p.positions,p.vertexCount*3),3));
    // GPU normalized signed-short normals keep the complete atlas compact in memory.
    g.setAttribute('normal',new T.BufferAttribute(new Int16Array(buffer,p.normals,p.vertexCount*3),3,true));
    if(p.uvs!==undefined)g.setAttribute('uv',new T.BufferAttribute(new Float32Array(buffer,p.uvs,p.vertexCount*2),2));
    if(p.colors!==undefined)g.setAttribute('color',new T.BufferAttribute(new Uint8Array(buffer,p.colors,p.vertexCount*3),3,true));
    g.setIndex(new T.BufferAttribute(new Uint32Array(buffer,p.indices,p.indexCount),1));
    if(p.sourceOffset)g.translate(p.sourceOffset[0],p.sourceOffset[1],p.sourceOffset[2]);
    g.boundingBox=bounds[i].clone();g.computeBoundingSphere();const pick=new T.Mesh(g,new T.MeshBasicMaterial({side:T.DoubleSide}));materials.push(pick.material);pick.matrixAutoUpdate=false;pickers[i]=pick;geometries.push(g);
    g.setAttribute('partIndex',new T.BufferAttribute(new Float32Array(p.vertexCount).fill(i),1));
    const materialKey=materialKeyFor(p);const list=groups.get(materialKey)??[];list.push(g);groups.set(materialKey,list);
   });
   groups.forEach((gs,key)=>{const geometry=mergeGeometries(gs,false);if(!geometry)throw new Error('Could not assemble anatomy geometry.');geometries.push(geometry);const mesh=new T.Mesh(geometry,mats.get(key));mesh.frustumCulled=false;mesh.renderOrder=2;scene.add(mesh);const ghost=new T.Mesh(geometry,ghostMats.get(key));ghost.frustumCulled=false;ghost.renderOrder=1;scene.add(ghost);});
   lastState=null;loaded++;onProgress(Math.round(loaded/atlas.chunks.length*100));dirty=true;
  };
  (async()=>{try{let cursor=0;await Promise.all(Array.from({length:3},async()=>{while(cursor<atlas.chunks.length){const i=cursor++;await loadChunk(i);}}));if(!disposed){ready=true;dirty=true;}}catch(e){if(!disposed)onError(e instanceof Error?e.message:'Could not load the anatomy.');}})();
  const directionFor=(view:string)=>view==='front'?new T.Vector3(0,.02,1):view==='back'?new T.Vector3(0,.02,-1):view==='side'?new T.Vector3(1,.02,0):view==='right'?new T.Vector3(-1,.02,0):view==='superior'?new T.Vector3(0,1,.001):view==='inferior'?new T.Vector3(0,-1,.001):new T.Vector3(.35,.06,1).normalize();
  const viewArea=()=>{
   const w=el.clientWidth,h=el.clientHeight,mobile=w<768,landscape=w>h&&h<=600;
   let left=20,right=w-70,top=90,bottom=h-24;
   if(mobile||landscape){bottom=h-150;right=w-60;}
   for(const selector of ['.layers-panel','.study-panel','.detail-sheet']){
    const panel=document.querySelector(selector);if(!panel||!panel.getClientRects().length)continue;
    const r=panel.getBoundingClientRect();
    if(mobile){bottom=Math.min(bottom,r.top-16);}else if(r.left<w/2){left=Math.max(left,r.right+24);}else{right=Math.min(right,r.left-24);}
   }
   return {w,h,left,right,top,bottom};
  };
  const fit=(view:string,extent=0)=>{
   const area=viewArea();camera.setViewOffset(area.w,area.h,area.w/2-(area.left+area.right)/2,area.h/2-(area.top+area.bottom)/2,area.w,area.h);
   const aspect=camera.aspect,mobile=el.clientWidth<768,normalDistance=mobile?Math.max(3.5,1.85*el.clientHeight/Math.max(160,el.clientHeight-400)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))):4;
   const reservedHeight=mobile?400:270;const availableAspect=Math.max(.35,(el.clientWidth-(mobile?40:340))/Math.max(160,el.clientHeight-reservedHeight));const atlasDistance=Math.max(packingHeight,packingWidth/availableAspect)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))*(el.clientHeight/Math.max(160,el.clientHeight-reservedHeight))*1.08;
   const distance=T.MathUtils.lerp(normalDistance,Math.max(.2,atlasDistance),extent);
   const direction=directionFor(view);
   controls.target.set(0,.85,0);const pivot=latest.current.isolate||amount<.05?selectionCenter(atlas.parts,latest.current.selected,data):null;if(pivot)controls.target.copy(pivot);let fittedDistance=distance;if((atlas.scope==='embryo'||latest.current.region&&latest.current.region!=='all')&&!pivot){const box=new T.Box3();atlas.parts.forEach((p,i)=>{if(!isSurfaceSystem(p.system)&&partVisible(p,latest.current))box.union(bounds[i]);});if(!box.isEmpty()){box.getCenter(controls.target);const size=box.getSize(new T.Vector3());const {w,h,left,right,top,bottom}=viewArea();camera.setViewOffset(w,h,w/2-(left+right)/2,h/2-(top+bottom)/2,w,h);fittedDistance=Math.max(size.y*h/Math.max(100,bottom-top),size.x*w/Math.max(100,right-left)/camera.aspect,size.z)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))*1.1;}}controls.maxDistance=Math.max(fittedDistance*1.25,.2);camera.position.copy(controls.target).addScaledVector(direction,fittedDistance);controls.update();dirty=true;
  };
  const resize=()=>{lastCamera=undefined;lastFocus=-1;layoutKey='';lastIsolate='';lastSelection='';lastState=null;renderer.setPixelRatio(Math.min(devicePixelRatio,el.clientWidth<768||el.clientHeight<600?1.5:2));camera.aspect=el.clientWidth/el.clientHeight;camera.updateProjectionMatrix();renderer.setSize(el.clientWidth,el.clientHeight);fit(amount>.04?'front':latest.current.view,Math.min(1,amount*6));};const observer=new ResizeObserver(resize);observer.observe(el);
  const raycaster=new T.Raycaster(),pointer=new T.Vector2(),tap=new PointerTap(),worldBox=new T.Box3(),hitPoint=new T.Vector3();
  let lastCoveringKey='',lastCoveringScan=0,lastCoveringHidden:SceneState['hidden'],lastCoveringDepth:SceneState['depthHidden'];
  const scanCovering=(s:SceneState)=>{
   if(!s.selected.length){coveringCallback.current?.([]);return;}
   const selected=new Set(s.selected),box=new T.Box3();
   atlas.parts.forEach((part,i)=>{if(selected.has(part.id)){const mesh=pickers[i];if(mesh)box.union(new T.Box3().setFromObject(mesh));}});
   if(box.isEmpty()){coveringCallback.current?.([]);return;}
   const center=box.getCenter(new T.Vector3()),size=box.getSize(new T.Vector3());
   const points=[center,center.clone().add(new T.Vector3(size.x*.18,size.y*.18,0)),center.clone().add(new T.Vector3(-size.x*.18,-size.y*.18,0))];
   const inspectionState={...s,selected:[],hidden:[],isolate:false};
   const candidates=new Map<number,number>();
   for(const point of points){
    const direction=point.clone().sub(camera.position),distance=direction.length();if(distance<.001)continue;
    raycaster.set(camera.position,direction.normalize());raycaster.near=0;raycaster.far=distance;
    let targetDistance=distance;
    atlas.parts.forEach((part,i)=>{if(!selected.has(part.id)||!pickers[i])return;const mesh=pickers[i]!;const hit=raycaster.intersectObject(mesh,false).find(hit=>clipPlane.distanceToPoint(hit.point)>=0);if(hit)targetDistance=Math.min(targetDistance,hit.distance);});
    raycaster.far=Math.max(0,targetDistance-.001);
    atlas.parts.forEach((part,i)=>{const mesh=pickers[i];if(!mesh||selected.has(part.id)||!partVisible(part,inspectionState)||isSurfaceSystem(part.system)&&(s.skinOpacity??0)<.01)return;if(!rotationPartIndices.has(i)){worldBox.copy(bounds[i]).translate(mesh.position);if(!raycaster.ray.intersectsBox(worldBox))return;}const hit=raycaster.intersectObject(mesh,false).find(hit=>clipPlane.distanceToPoint(hit.point)>=0);if(hit)candidates.set(i,Math.min(candidates.get(i)??Infinity,hit.distance));});
   }
   coveringCallback.current?.([...candidates].sort((a,b)=>a[1]-b[1]).map(([i])=>atlas.parts[i].id));
  };
  const selectionPivot=()=>amount>.05&&!latest.current.isolate?null:selectionCenter(atlas.parts,latest.current.selected,data);
  const anchorSelection=()=>{const pivot=selectionPivot();if(pivot){controls.target.copy(pivot);controls.update();}};
  const syncSelectedRotation=()=>{
   const box=new T.Box3();for(const i of rotationPartIndices)box.union(bounds[i].clone().translate(new T.Vector3(data[i*4],data[i*4+1],data[i*4+2])));
   const pivot=box.isEmpty()?null:box.getCenter(new T.Vector3());pivotUniform.value.copy(pivot??new T.Vector3());rotationUniform.value.set(selectedRotation.x,selectedRotation.y,selectedRotation.z,selectedRotation.w);
   if(pivot)for(const i of rotationPartIndices){const mesh=pickers[i];if(!mesh)continue;
    mesh.quaternion.copy(selectedRotation);mesh.position.set(data[i*4],data[i*4+1],data[i*4+2]).sub(pivot).applyQuaternion(selectedRotation).add(pivot);mesh.updateMatrix();mesh.updateMatrixWorld(true);
   }
   dirty=true;
  };
  const pressPoints=new Map<number,{x:number;y:number;lastX:number;lastY:number;anchored:boolean;rotating:boolean}>();
  const solidContext=()=> (latest.current.skinOpacity??.1)<.95&&atlas.parts.some((p,i)=>!isSurfaceSystem(p.system)&&data[i*4+3]>.5);
  const canPick=(i:number,hasSolid:boolean)=>{
   const part=atlas.parts[i],selected=latest.current.selected.includes(part.id);
   if(data[i*4+3]<.5||hasSolid&&isSurfaceSystem(part.system)&&!selected)return false;
   const alpha=selected?1:contextUniform.value*(isSurfaceSystem(part.system)?skinUniform.value:1)*(atlas.materials?.[part.material??'']?.opacity??1);
   return alpha>.001&&(selected||!latest.current.selected.length||alpha>=.999||selectFaded.current);
  };
  const down=(e:PointerEvent)=>{const rotating=amount>.04&&!latest.current.isolate&&rotationPartIndices.size>0&&e.button===0;pressPoints.set(e.pointerId,{x:e.clientX,y:e.clientY,lastX:e.clientX,lastY:e.clientY,anchored:false,rotating});hover.hidden=true;tap.down(e.pointerId,e.clientX,e.clientY,e.pointerType==='touch'?16:10);if(rotating){e.preventDefault();e.stopImmediatePropagation();renderer.domElement.setPointerCapture(e.pointerId);}};
  const move=(e:PointerEvent)=>{tap.move(e.pointerId,e.clientX,e.clientY);const press=pressPoints.get(e.pointerId);if(press?.rotating){e.preventDefault();e.stopImmediatePropagation();const dx=e.clientX-press.lastX,dy=e.clientY-press.lastY;press.lastX=e.clientX;press.lastY=e.clientY;if(dx||dy){const up=new T.Vector3(0,1,0).applyQuaternion(camera.quaternion),right=new T.Vector3(1,0,0).applyQuaternion(camera.quaternion);selectedRotation.premultiply(new T.Quaternion().setFromAxisAngle(up,dx*.009)).premultiply(new T.Quaternion().setFromAxisAngle(right,dy*.009)).normalize();syncSelectedRotation();}hover.hidden=true;renderer.domElement.style.cursor='grabbing';return;}if(e.buttons&&press&&!press.anchored&&Math.hypot(e.clientX-press.x,e.clientY-press.y)>8){press.anchored=true;anchorSelection();}if(e.buttons||!ready||amount<.5||e.pointerType==='touch'){hover.hidden=true;renderer.domElement.style.cursor=e.buttons?'grabbing':'grab';return;}const rect=el.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top,hasSolid=solidContext(),index=findTarget(x,y,12,i=>canPick(i,hasSolid));hover.hidden=index<0;renderer.domElement.style.cursor=index<0?'grab':'pointer';if(index>=0){hover.textContent=structureName(atlas.parts[index].name);hover.style.left=`${Math.max(8,Math.min(x+14,el.clientWidth-260))}px`;hover.style.top=`${Math.max(8,Math.min(y+18,el.clientHeight-55))}px`;}};
  const cancel=(e:PointerEvent)=>{pressPoints.delete(e.pointerId);tap.cancel(e.pointerId);};
  const up=(e:PointerEvent)=>{
   const press=pressPoints.get(e.pointerId);pressPoints.delete(e.pointerId);if(press?.rotating){e.preventDefault();e.stopImmediatePropagation();}const validTap=tap.up(e.pointerId,e.clientX,e.clientY);if(!validTap||!ready)return;const rect=renderer.domElement.getBoundingClientRect();pointer.set((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1);raycaster.setFromCamera(pointer,camera);
   // Covering-tissue scans shorten raycaster.far; each click must trace the
   // entire visible anatomy, including structures beyond the last scan point.
   raycaster.near=0;raycaster.far=Infinity;
   let nearest=Infinity,found=-1;const hasSolid=solidContext();
   pickers.forEach((mesh,i)=>{if(!mesh||!canPick(i,hasSolid))return;if(!rotationPartIndices.has(i)){worldBox.copy(bounds[i]).translate(mesh.position);if(!raycaster.ray.intersectBox(worldBox,hitPoint))return;}const hit=raycaster.intersectObject(mesh,false).find(h=>clipPlane.distanceToPoint(h.point)>=0);const priority=latest.current.selected.includes(atlas.parts[i].id)&&contextUniform.value<1&&!selectFaded.current?1000:0;if(hit&&hit.distance-priority<nearest){nearest=hit.distance-priority;found=i;}});
   if(found<0&&amount>.45)found=findTarget(e.clientX-rect.left,e.clientY-rect.top,e.pointerType==='touch'?24:16,i=>canPick(i,hasSolid));if(found>=0){hover.hidden=true;renderer.domElement.style.cursor='grab';select.current(atlas.parts[found].id,e.shiftKey);}
  };
  renderer.domElement.addEventListener('pointerdown',down,true);renderer.domElement.addEventListener('pointermove',move);renderer.domElement.addEventListener('pointerup',up);renderer.domElement.addEventListener('pointercancel',cancel);
  const drawLabels=(s:SceneState)=>{
   labelLayer.replaceChildren();if(s.labels===false)return;const w=el.clientWidth,h=el.clientHeight;labelLayer.setAttribute('viewBox',`0 0 ${w} ${h}`);
   if(amount>.04&&explosionLayout&&!s.isolate){
    const progress=amount*explosionLayout.steps,from=Math.min(explosionLayout.steps,Math.floor(progress)),to=Math.min(explosionLayout.steps,from+1),mix=progress-from;
    const depth=Math.min(4,Math.max(1,Math.floor(progress+.05))),first=new Map(explosionLayout.stages[from].groups.map(group=>[group.id,group]));
    const candidates=explosionLayout.stages[to].groups.filter(group=>group.depth<=depth),names=new Set(candidates.map(group=>group.id)),parents=new Set<string>();
    for(const group of candidates){let parent=group.id.lastIndexOf('/');while(parent>0){const id=group.id.slice(0,parent);if(names.has(id))parents.add(id);parent=id.lastIndexOf('/');}}
    const area=viewArea(),occupied:{left:number;right:number;top:number;bottom:number}[]=[];
    for(const group of candidates.filter(group=>!parents.has(group.id)).sort((a,b)=>b.count-a.count).slice(0,w<768?50:180)){
     const previous=first.get(group.id)??group;
     const x=T.MathUtils.lerp(previous.x,group.x,mix),y=T.MathUtils.lerp(previous.y,group.y,mix),z=T.MathUtils.lerp(previous.z,group.z,mix),height=T.MathUtils.lerp(previous.height,group.height,mix);
     projected.set(x,y+height/2,z).project(camera);if(projected.z< -1||projected.z>1)continue;
     const screenX=(projected.x+1)*w/2,screenY=(1-projected.y)*h/2-10;
     const name=group.name.length>32?`${group.name.slice(0,31)}…`:group.name,labelWidth=Math.min(220,name.length*6.2+14);
     const box={left:screenX-labelWidth/2,right:screenX+labelWidth/2,top:screenY-17,bottom:screenY+2};
     if(box.left<area.left||box.right>area.right||box.top<area.top||box.bottom>area.bottom||occupied.some(other=>box.left<other.right+4&&box.right>other.left-4&&box.top<other.bottom+3&&box.bottom>other.top-3))continue;
     occupied.push(box);
     const element=document.createElementNS('http://www.w3.org/2000/svg','g');element.setAttribute('class','explosion-group-label');
     const background=document.createElementNS('http://www.w3.org/2000/svg','rect');background.setAttribute('x',String(box.left));background.setAttribute('y',String(box.top));background.setAttribute('width',String(labelWidth));background.setAttribute('height','19');background.setAttribute('rx','5');element.appendChild(background);
     const label=document.createElementNS('http://www.w3.org/2000/svg','text');label.setAttribute('x',String(screenX));label.setAttribute('y',String(screenY-3));label.setAttribute('text-anchor','middle');label.textContent=name;element.appendChild(label);labelLayer.appendChild(element);
    }
   }
   let row=0;for(const id of s.selected.slice(0,20)){const i=atlas.parts.findIndex(p=>p.id===id),mesh=pickers[i];if(i<0||!mesh)continue;
    let anchor=labelAnchors.get(i);if(!anchor){const g=mesh.geometry,positions=g.getAttribute('position'),indices=g.index!;let best=Infinity;const center=centers[i];for(let j=0;j<indices.count;j+=Math.max(1,Math.floor(indices.count/1800))*3){const point=new T.Vector3();for(let k=0;k<3;k++)point.add(new T.Vector3().fromBufferAttribute(positions,indices.getX(Math.min(j+k,indices.count-1))));point.multiplyScalar(1/3);if(clipPlane.distanceToPoint(point.clone().add(mesh.position))<0)continue;const d=point.distanceToSquared(center);if(d<best){best=d;anchor=point;}}if(!anchor)continue;labelAnchors.set(i,anchor);}
    const world=anchor.clone().applyMatrix4(mesh.matrixWorld);if(clipPlane.distanceToPoint(world)<0)continue;const point=world.clone().project(camera);
    if(point.z< -1||point.z>1||Math.abs(point.x)>1||Math.abs(point.y)>1)continue;
    const x=(point.x+1)*w/2,y=(1-point.y)*h/2,labelX=Math.min(w-165,Math.max(w<768?18:310,x+35)),labelY=Math.max(w<768?235:145,Math.min(h-200,y-25+row*23));row++;
    const make=(name:string,attrs:Record<string,string>)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',name);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,v));labelLayer.appendChild(e);return e;};
    make('line',{x1:String(x),y1:String(y),x2:String(labelX),y2:String(labelY)});make('circle',{cx:String(x),cy:String(y),r:'3'});const text=make('text',{x:String(labelX+5),y:String(labelY-5)});text.textContent=structureName(atlas.parts[i].name);
   }
  };
  const clock=new T.Clock();let lastExtent=-1,lastCamera:number[]|undefined,lastFocus=0,lastAreaKey='',lastRotationKey='';
  const animate=()=>{
   if(disposed)return;frame=requestAnimationFrame(animate);const dt=Math.min(clock.getDelta(),.05),s=latest.current;
   const nextContext=s.selected.length?(s.contextOpacity??1):1;if(contextUniform.value!==nextContext||skinUniform.value!==(s.skinOpacity??.1)){for(const material of mats.values()){material.transparent=nextContext<.999;material.needsUpdate=true;}contextUniform.value=nextContext;skinUniform.value=s.skinOpacity??.1;dirty=true;}
   const changed=lastState!==s,selectionCleared=!!lastState?.selected.length&&!s.selected.length;
   if(changed){const a=viewArea(),key=[a.w,a.h,a.left,a.right,a.top,a.bottom].join(':');if(key!==lastAreaKey){if(!selectionCleared)camera.setViewOffset(a.w,a.h,a.w/2-(a.left+a.right)/2,a.h/2-(a.top+a.bottom)/2,a.w,a.h);if(lastAreaKey&&s.focus&&!selectionCleared)lastFocus=-1;lastAreaKey=key;dirty=true;}}

   if(changed){labelAnchors.clear();const box=new T.Box3();atlas.parts.forEach((p,i)=>{if(p.system!=='schematic'&&!isSurfaceSystem(p.system)&&partVisible(p,s))box.union(bounds[i]);});const section=s.section;if(section?.enabled&&!box.isEmpty()){const axis=section.axis==='axial'?1:section.axis==='sagittal'?0:2;const sign=section.flip?1:-1;clipPlane.normal.set(0,0,0).setComponent(axis,sign);clipPlane.constant=-sign*T.MathUtils.lerp(box.min.getComponent(axis),box.max.getComponent(axis),section.position);}else clipPlane.constant=10000;dirty=true;}
   const moving=Math.abs(amount-s.explode)>.0001;
   if(moving){amount=T.MathUtils.damp(amount,s.explode,8,dt);dirty=true;}
   if(changed||moving||lastExtent<0){
    const selection=new Set(s.selected);
    const visibleParts=atlas.parts.filter(p=>partVisible(p,s));
    const nextLayoutKey=visibleParts.map(p=>p.id).join(',')+':'+camera.aspect.toFixed(3)+':'+hierarchyRef.current;
    const layoutChanged=nextLayoutKey!==layoutKey;
    if(layoutChanged){explosionLayout=createHierarchicalExplosionLayout(atlas,anatomyEntries,new Set(visibleParts.map(part=>part.id)),hierarchyRef.current,camera.aspect);layoutKey=nextLayoutKey;stepsCallback.current?.(explosionLayout.steps);}
    const stages=explosionLayout!.stages,progress=amount*explosionLayout!.steps,from=Math.min(stages.length-1,Math.floor(progress)),to=Math.min(stages.length-1,from+1),blend=progress-from;
    packingWidth=T.MathUtils.lerp(stages[from].width,stages[to].width,blend);packingHeight=T.MathUtils.lerp(stages[from].height,stages[to].height,blend);
    if(layoutChanged&&amount>.001&&!s.isolate&&!selectionCleared)fit('front',Math.min(1,amount*6));
    const clusterStageIndex=Math.min(stages.length-1,Math.floor(progress+.05)),clusterStage=stages[clusterStageIndex],selectedClusters=new Set<number>();
    if(clusterStageIndex>0)for(const id of s.selected){const i=partIndices.get(id),cluster=i===undefined?-1:clusterStage.clusterIds[i];if(cluster>=0)selectedClusters.add(cluster);}
    const rotationKey=s.selected.join(',')+':'+s.reset+':'+clusterStageIndex+':'+[...selectedClusters].join(',');
    if(layoutChanged||rotationKey!==lastRotationKey||s.explode<=.04)selectedRotation.identity();lastRotationKey=rotationKey;
    rotationPartIndices.clear();rotationData.fill(0);
    for(const cluster of selectedClusters)for(const i of clusterStage.clusters[cluster]){rotationPartIndices.add(i);rotationData[i*4]=255;}

    atlas.parts.forEach((p,i)=>{
     const c=centers[i],j=i*3,dx=T.MathUtils.lerp(stages[from].positions[j],stages[to].positions[j],blend)-c.x,dy=T.MathUtils.lerp(stages[from].positions[j+1],stages[to].positions[j+1],blend)-c.y,dz=T.MathUtils.lerp(stages[from].positions[j+2],stages[to].positions[j+2],blend)-c.z;
     const selected=selection.has(p.id);data.set([dx,dy,dz,partVisible(p,s)?1:0],i*4);selectedData[i*4]=selected?255:0;
     const mesh=pickers[i];if(mesh){mesh.quaternion.identity();mesh.position.set(dx,dy,dz);mesh.updateMatrix();mesh.updateMatrixWorld(true);}
    });syncSelectedRotation();partTexture.needsUpdate=true;selectionTexture.needsUpdate=true;rotationTexture.needsUpdate=true;lastState=s;lastExtent=amount;dirty=true;
   }
   if(s.view!==lastView||s.reset!==lastReset||s.region!==lastRegion){if(s.focus)lastFocus=-1;fit(amount>.04?'front':s.view,Math.min(1,amount*6));lastView=s.view;lastReset=s.reset;lastRegion=s.region??'all';}
   if(moving&&!s.isolate)fit(amount>.04?'front':s.view,Math.min(1,amount*6));
   const isolateKey=s.isolate?s.selected.join(',')+':'+s.reset+':'+s.inspectorOpen+':'+camera.aspect:'';
   const focusChanged=(s.focus??0)!==lastFocus;lastFocus=s.focus??0;
   if(isolateKey!==lastIsolate||(s.isolate&&moving)||focusChanged&&!!s.focus){
    if(s.isolate||focusChanged&&!!s.focus){const box=new T.Box3();atlas.parts.forEach((p,i)=>{if(s.selected.includes(p.id)&&pickers[i])box.union(new T.Box3().setFromObject(pickers[i]!));});
     if(!box.isEmpty()){const center=box.getCenter(new T.Vector3()),size=box.getSize(new T.Vector3());const {w,h,left,right,top,bottom}=viewArea();const availableWidth=Math.max(150,right-left),availableHeight=Math.max(40,bottom-top);camera.setViewOffset(w,h,w/2-(left+right)/2,h/2-(top+bottom)/2,w,h);const distance=Math.max(.006,Math.max(size.y*h/availableHeight,size.x*w/availableWidth/camera.aspect,size.z)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))*1.35);controls.maxDistance=Math.max(.12,distance*1.6);controls.target.copy(center);camera.position.copy(center).add(directionFor(s.view).normalize().multiplyScalar(distance));controls.update();dirty=true;}
    }else if(lastIsolate&&!selectionCleared){camera.clearViewOffset();fit(s.view,Math.min(1,amount*6));}
    lastIsolate=isolateKey;
   }
   const selectionKey=s.selected.join(',')+':'+s.reset;if(selectionKey!==lastSelection||moving&&s.selected.length>0){const pivot=selectionPivot();if(pivot){controls.target.copy(pivot);controls.update();}lastSelection=selectionKey;dirty=true;}
   const canOrbit=amount<=.04||s.isolate;controls.enableRotate=canOrbit;controls.mouseButtons.LEFT=canOrbit?T.MOUSE.ROTATE:T.MOUSE.PAN;controls.touches.ONE=canOrbit?T.TOUCH.ROTATE:T.TOUCH.PAN;controls.autoRotate=s.rotate&&!s.isolate&&canOrbit;controls.autoRotateSpeed=.65;controls.update();if(controls.autoRotate)dirty=true;
   if(!s.camera)lastCamera=undefined;
   if(s.camera&&s.camera!==lastCamera){camera.position.fromArray(s.camera);controls.target.fromArray(s.camera,3);camera.clearViewOffset();if(s.camera.length===8)camera.setViewOffset(el.clientWidth,el.clientHeight,s.camera[6]*el.clientWidth,s.camera[7]*el.clientHeight,el.clientWidth,el.clientHeight);controls.maxDistance=Math.max(controls.maxDistance,camera.position.distanceTo(controls.target)*1.05);controls.update();lastCamera=s.camera;dirty=true;}
   // Keep depth precision near the anatomy as the camera moves. A fixed tiny near
   // plane causes close skin and muscle surfaces to fight at whole-body distances.
   const near=Math.max(.00001,Math.min(.05,camera.position.distanceTo(controls.target)/500));
   if(Math.abs(camera.near-near)>near*.01){camera.near=near;camera.updateProjectionMatrix();dirty=true;}
   if(dirty){
    renderer.render(scene,camera);cameraCallback.current?.(cameraValues());drawLabels(s);targets=[];
    if(amount>.45){
     const hasSolid=(latest.current.skinOpacity??.1)<.95&&atlas.parts.some((p,i)=>!isSurfaceSystem(p.system)&&data[i*4+3]>.5);
     atlas.parts.forEach((p,i)=>{
      if(data[i*4+3]<.5||(!s.selected.includes(p.id)&&contextUniform.value<.001)||(hasSolid&&isSurfaceSystem(p.system)&&!s.selected.includes(p.id)))return;
      const mesh=pickers[i];if(!mesh)return;let left=Infinity,right=-Infinity,top=Infinity,bottom=-Infinity;
      for(let corner=0;corner<8;corner++){
       projected.set(p.bounds[(corner&1)?1:0][0],p.bounds[(corner&2)?1:0][1],p.bounds[(corner&4)?1:0][2]).applyMatrix4(mesh.matrixWorld).project(camera);
       const x=(projected.x+1)*el.clientWidth/2,y=(1-projected.y)*el.clientHeight/2;left=Math.min(left,x);right=Math.max(right,x);top=Math.min(top,y);bottom=Math.max(bottom,y);
      }
      projected.copy(centers[i]).applyMatrix4(mesh.matrixWorld).project(camera);if(projected.z< -1||projected.z>1)return;
      targets.push({index:i,x:(projected.x+1)*el.clientWidth/2,y:(1-projected.y)*el.clientHeight/2,left,right,top,bottom});
     });
    }
    dirty=false;
   }
   if(ready){const key=[s.selected.join(','),s.visible.join(','),s.region,s.skinOpacity,s.explode,...selectedRotation.toArray().map(v=>v.toFixed(2)),...camera.position.toArray().map(v=>v.toFixed(2)),...controls.target.toArray().map(v=>v.toFixed(2))].join(':');if((key!==lastCoveringKey||s.hidden!==lastCoveringHidden||s.depthHidden!==lastCoveringDepth)&&performance.now()-lastCoveringScan>250){scanCovering(s);lastCoveringKey=key;lastCoveringHidden=s.hidden;lastCoveringDepth=s.depthHidden;lastCoveringScan=performance.now();}}

  };animate();
  onCapture?.(async()=>{
   if(!ready||disposed)throw new Error('The 3D view is still loading.');
   renderer.render(scene,camera);
   const png=renderer.domElement.toDataURL('image/png');
   return (await fetch(png)).blob();
  });
  const contextLost=(e:Event)=>{e.preventDefault();onError('The 3D session was paused by your device. Reload to continue.');};renderer.domElement.addEventListener('webglcontextlost',contextLost);
  return()=>{onCapture?.(null);disposed=true;abort.abort();cancelAnimationFrame(frame);observer.disconnect();theme.removeEventListener('change',applyTheme);controls.dispose();geometries.forEach(g=>g.dispose());materials.forEach(m=>m.dispose());textures.forEach(t=>t.dispose());scene.traverse(o=>{if(o instanceof T.Mesh&&!geometries.includes(o.geometry)){o.geometry.dispose();const ms=Array.isArray(o.material)?o.material:[o.material];ms.forEach(m=>m.dispose());}});env.dispose();partTexture.dispose();selectionTexture.dispose();rotationTexture.dispose();hover.remove();labelLayer.remove();renderer.dispose();renderer.domElement.remove();};
 },[atlas]);
 return <div className="scene" ref={host}/>;
}
