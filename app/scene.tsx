import {assetUrl} from './asset-url';
import {inRegion,partVisible,sectionPartVisible} from './viewer-state';
import {useEffect,useMemo,useRef} from 'react';
import * as T from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import {RoomEnvironment} from 'three/examples/jsm/environments/RoomEnvironment.js';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';
import {createHierarchicalExplosionLayout,type ExplodeHierarchy,type HierarchicalExplosionLayout} from './hierarchical-explosion';
import {hierarchyEntries} from './anatomy-hierarchy';
import {resolveGuestHierarchy,type GuestHierarchy} from './guest-hierarchy';
import {decodeModelResponse} from './model-download';
import {selectionCenter} from './camera-pivot';
import {PointerTap} from './pointer-tap';
import {sectionDot,sectionFraction,sectionNormal,sectionPlaneEquation,sectionPoint,sectionRange,type Bounds3,type Point3} from './section-plane';
import {sectionCapGeometry} from './section-cap';
import {sectionCapProfile} from './section-cap-profile';
import {sectionCapOpacity} from './section-opacity';
import {sectionSetBounds} from './section-set';
import {enabledSections} from './section-stack';
import {SYSTEMS,structureName,type Atlas,type SceneState,type SystemId} from './anatomy';
import {isSkinPart} from './depth-layers';
interface Props {atlas:Atlas;state:SceneState;hierarchy:ExplodeHierarchy;guestHierarchy?:GuestHierarchy;sectionTool?:boolean;onSectionPosition?:(position:number)=>void;onSelect:(id:string,add?:boolean)=>void;onCamera?:(camera:number[])=>void;onCovering?:(ids:string[])=>void;onExplosionSteps?:(steps:number)=>void;onCapture?:(capture:(()=>Promise<Blob>)|null)=>void;onProgress:(n:number)=>void;onError:(s:string)=>void}
type SceneMaterial={color:number[];map?:string;normalMap?:string;normalScale?:number;roughness?:number;metalness?:number;opacity?:number;vertexColors?:boolean;proceduralMuscle?:boolean};
export default function AnatomyScene({atlas,state,hierarchy,guestHierarchy,sectionTool=false,onSectionPosition,onSelect,onCamera,onCovering,onExplosionSteps,onCapture,onProgress,onError}:Props){
 const guestNodes=useMemo(()=>guestHierarchy?resolveGuestHierarchy(atlas,guestHierarchy):undefined,[atlas,guestHierarchy]);
 const host=useRef<HTMLDivElement>(null),latest=useRef(state),hierarchyRef=useRef(hierarchy),guestNodesRef=useRef(guestNodes),guestKeyRef=useRef(''),sectionToolRef=useRef(sectionTool),sectionPositionRef=useRef(onSectionPosition),select=useRef(onSelect),cameraCallback=useRef(onCamera),coveringCallback=useRef(onCovering),stepsCallback=useRef(onExplosionSteps);
 latest.current=state;hierarchyRef.current=hierarchy;guestNodesRef.current=guestNodes;guestKeyRef.current=guestHierarchy?JSON.stringify(guestHierarchy):'';sectionToolRef.current=sectionTool;sectionPositionRef.current=onSectionPosition;select.current=onSelect;cameraCallback.current=onCamera;coveringCallback.current=onCovering;stepsCallback.current=onExplosionSteps;
 useEffect(()=>{
  const el=host.current!;let disposed=false,frame=0,dirty=true,ready=false,lastView='',lastReset=-1,lastIsolate='',lastSelection='',lastRegion='',layoutKey='',amount=0;
  let lastState:SceneState|null=null;
  const abort=new AbortController();
  let renderer:T.WebGLRenderer;
  try{renderer=new T.WebGLRenderer({antialias:true,alpha:false,stencil:true,powerPreference:'high-performance'});}catch{onError('This browser could not start the 3D viewer. Please try a browser with WebGL enabled.');return;}
  const referenceModel=atlas.version==='Anatomy Atlas adapted GLB local import';
  const stageInitialCut=referenceModel&&enabledSections(latest.current).length>0;
  let stagedCapsReady=!stageInitialCut;
  renderer.setPixelRatio(Math.min(devicePixelRatio,innerWidth<768?1.5:2));renderer.setClearColor('#f2f3f3');renderer.outputColorSpace=T.SRGBColorSpace;renderer.toneMapping=T.ACESFilmicToneMapping;renderer.toneMappingExposure=referenceModel ? .78 : .96;el.appendChild(renderer.domElement);
  renderer.domElement.setAttribute('aria-label',atlas.scope==='cell'?'Interactive human cell. Drag to orbit, pinch or scroll to zoom, and tap a component to inspect it.':'Interactive human anatomy. Drag to orbit, pinch or scroll to zoom, and tap a structure to inspect it.');
  const scene=new T.Scene(),camera=new T.PerspectiveCamera(34,1,.01,100),controls=new OrbitControls(camera,renderer.domElement);
  camera.position.set(1.4,1.05,3.6);controls.target.set(0,.85,0);controls.zoomToCursor=true;controls.enableDamping=true;controls.dampingFactor=.085;controls.minDistance=.003;controls.maxDistance=6;controls.zoomSpeed=1.25;controls.maxPolarAngle=Math.PI-.001;controls.addEventListener('change',()=>{dirty=true;});
  const pmrem=new T.PMREMGenerator(renderer),room=new RoomEnvironment(),env=pmrem.fromScene(room,.04);scene.environment=env.texture;scene.environmentIntensity=referenceModel ? .7 : 1;room.dispose();pmrem.dispose();
  scene.add(new T.HemisphereLight(0xffffff,0x8c97a0,referenceModel ? .55 : .8));
  const key=new T.DirectionalLight(0xfffaf4,referenceModel?1.65:2.05);key.position.set(-2,4,3);scene.add(key);
  const rim=new T.DirectionalLight(0xe9f0ff,referenceModel ? .65 : 1.1);rim.position.set(2,2,-3);scene.add(rim);
  const theme=window.matchMedia('(prefers-color-scheme: dark)');
  const applyTheme=()=>{renderer.setClearColor(theme.matches?'#141b23':'#f2f3f3');dirty=true;};
  applyTheme();theme.addEventListener('change',applyTheme);
  const width=T.MathUtils.ceilPowerOfTwo(atlas.parts.length),data=new Float32Array(width*4),partTexture=new T.DataTexture(data,width,1,T.RGBAFormat,T.FloatType);partTexture.needsUpdate=true;
  const selectedData=new Uint8Array(width*4),selectionTexture=new T.DataTexture(selectedData,width,1);selectionTexture.needsUpdate=true;
  const previewData=new Uint8Array(width*4),previewTexture=new T.DataTexture(previewData,width,1);previewTexture.needsUpdate=true;
  const rotationData=new Uint8Array(width*4),rotationTexture=new T.DataTexture(rotationData,width,1);rotationTexture.needsUpdate=true;
  const selectedRotation=new T.Quaternion(),rotationUniform={value:new T.Vector4(0,0,0,1)},pivotUniform={value:new T.Vector3()},rotationPartIndices=new Set<number>(),partIndices=new Map(atlas.parts.map((part,i)=>[part.id,i]));
  const materials:T.Material[]=[],geometries:T.BufferGeometry[]=[],pickers:(T.Mesh|undefined)[]=[],centers=atlas.parts.map(p=>new T.Vector3().fromArray(p.bounds[0]).add(new T.Vector3().fromArray(p.bounds[1])).multiplyScalar(.5));
  const bounds=atlas.parts.map(p=>new T.Box3(new T.Vector3().fromArray(p.bounds[0]),new T.Vector3().fromArray(p.bounds[1])));
  const anatomyEntries=hierarchyEntries(atlas);let explosionLayout:HierarchicalExplosionLayout|null=null;
  let packingWidth=1,packingHeight=1;
  const hover=document.createElement('div');hover.className='part-hover';hover.setAttribute('role','tooltip');hover.hidden=true;
  const hoverLabel=document.createElement('span');hover.appendChild(hoverLabel);el.appendChild(hover);
  type Target={index:number;x:number;y:number;left:number;right:number;top:number;bottom:number};let targets:Target[]=[];
  const projected=new T.Vector3();
  const findTarget=(x:number,y:number,radius:number,eligible:(index:number)=>boolean=()=>true)=>{
   let best=-1,score=Infinity;
   for(const t of targets){if(!eligible(t.index))continue;const dx=Math.max(t.left-x,0,x-t.right),dy=Math.max(t.top-y,0,y-t.bottom),distance=Math.hypot(dx,dy);if(distance>radius)continue;const candidate=distance+Math.hypot(t.x-x,t.y-y)*.025;if(candidate<score){score=candidate;best=t.index;}}
   return best;
  };
  const clipPlanes=[new T.Plane(new T.Vector3(0,1,0),10000),new T.Plane(new T.Vector3(0,1,0),10000)];renderer.localClippingEnabled=true;
  const keptBySections=(point:T.Vector3)=>clipPlanes.every(plane=>plane.distanceToPoint(point)>=-.00001);
  const guideGeometry=new T.PlaneGeometry(1,1),guideFill=new T.Mesh(guideGeometry,new T.MeshBasicMaterial({color:0x5b9d9c,side:T.DoubleSide,transparent:true,opacity:.09,depthWrite:false}));
  const guideBorder=new T.LineSegments(new T.EdgesGeometry(guideGeometry),new T.LineBasicMaterial({color:0x478f8d,transparent:true,opacity:.72,depthTest:false}));
  guideFill.renderOrder=5;guideBorder.renderOrder=6;scene.add(guideFill,guideBorder);
  const previewScene=new T.Scene(),previewCamera=new T.PerspectiveCamera(28,1,.01,50),previewTarget=new T.WebGLRenderTarget(1,1);
  previewScene.background=new T.Color('#1b252b');previewScene.add(new T.HemisphereLight(0xffffff,0x77808a,2));
  const previewLight=new T.DirectionalLight(0xfff8ed,2.2);previewLight.position.set(-2,3,4);previewScene.add(previewLight);
  const previewMaterials=new Map<string,T.MeshStandardMaterial>();
  const hasSurfacePreview=atlas.parts.some(isSkinPart);
  let previewAzimuth=.48,previewElevation=.08,previewCanvas:HTMLCanvasElement|null=null,previewDirty=true,previewPointer=-1,previewX=0,previewY=0;
  const previewPlanes=Array.from({length:2},(_,index)=>{
   const geometry=new T.PlaneGeometry(1,1),color=index?0xf2aa66:0x49b6b1;
   const fill=new T.Mesh(geometry,new T.MeshBasicMaterial({color,side:T.DoubleSide,transparent:true,opacity:.32,depthWrite:false}));
   const edge=new T.LineSegments(new T.EdgesGeometry(geometry),new T.LineBasicMaterial({color,depthTest:false}));
   fill.renderOrder=5+index*2;edge.renderOrder=6+index*2;previewScene.add(fill,edge);return {geometry,fill,edge};
  });
  const previewDown=(event:PointerEvent)=>{previewPointer=event.pointerId;previewX=event.clientX;previewY=event.clientY;previewCanvas?.setPointerCapture(event.pointerId);};
  const previewMove=(event:PointerEvent)=>{if(event.pointerId!==previewPointer)return;previewAzimuth+=(event.clientX-previewX)*.012;previewElevation=Math.max(-1.2,Math.min(1.2,previewElevation+(event.clientY-previewY)*.009));previewX=event.clientX;previewY=event.clientY;previewDirty=true;};
  const previewUp=(event:PointerEvent)=>{if(event.pointerId===previewPointer)previewPointer=-1;};
  const attachPreview=()=>{const next=document.querySelector<HTMLCanvasElement>('.section-preview canvas');if(next===previewCanvas)return;
   previewCanvas?.removeEventListener('pointerdown',previewDown);previewCanvas?.removeEventListener('pointermove',previewMove);previewCanvas?.removeEventListener('pointerup',previewUp);previewCanvas?.removeEventListener('pointercancel',previewUp);
   previewCanvas=next;if(next){next.addEventListener('pointerdown',previewDown);next.addEventListener('pointermove',previewMove);next.addEventListener('pointerup',previewUp);next.addEventListener('pointercancel',previewUp);}previewDirty=true;
  };
  const drawPreview=(s:SceneState)=>{if(!previewCanvas)return;const canvas=previewCanvas,w=Math.max(1,Math.round(canvas.clientWidth*Math.min(devicePixelRatio,2))),h=Math.max(1,Math.round(canvas.clientHeight*Math.min(devicePixelRatio,2)));
   if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}
   if(previewTarget.width!==w||previewTarget.height!==h)previewTarget.setSize(w,h);
   const previewBox=new T.Box3(new T.Vector3().fromArray(sectionBounds.min),new T.Vector3().fromArray(sectionBounds.max));
   const previewCenter=previewBox.getCenter(new T.Vector3()),previewSize=previewBox.getSize(new T.Vector3());
   previewCamera.aspect=w/h;
   const tangent=Math.tan(T.MathUtils.degToRad(previewCamera.fov/2));
   const previewDistance=Math.max(.04,previewSize.y/(2*tangent),previewSize.x/(2*tangent*previewCamera.aspect),previewSize.z/(2*tangent))*1.35;
   previewCamera.near=Math.max(.00001,previewDistance/10000);previewCamera.far=Math.max(50,previewDistance*10);
   previewCamera.position.copy(previewCenter).add(new T.Vector3(Math.sin(previewAzimuth)*Math.cos(previewElevation),Math.sin(previewElevation),Math.cos(previewAzimuth)*Math.cos(previewElevation)).multiplyScalar(previewDistance));previewCamera.lookAt(previewCenter);previewCamera.updateProjectionMatrix();
   const stack=(s.sections?.length?s.sections:[s.section]).filter((value):value is NonNullable<typeof value>=>!!value).slice(0,2);
   previewPlanes.forEach(({fill,edge},index)=>{const section=stack[index],visible=!!section?.enabled;fill.visible=visible;edge.visible=visible;if(!visible)return;
    const point=sectionPoint(sectionBounds,section),normal=sectionNormal(section),quaternion=new T.Quaternion().setFromUnitVectors(new T.Vector3(0,0,1),new T.Vector3().fromArray(normal)),span=Math.max(previewSize.x,previewSize.z,previewSize.y*.72,.04)*1.16;
    fill.position.fromArray(point);edge.position.fromArray(point);fill.quaternion.copy(quaternion);edge.quaternion.copy(quaternion);fill.scale.set(span,span,1);edge.scale.set(span,span,1);
   });
   renderer.setRenderTarget(previewTarget);renderer.render(previewScene,previewCamera);
   const pixels=new Uint8Array(w*h*4),flipped=new Uint8ClampedArray(w*h*4);renderer.readRenderTargetPixels(previewTarget,0,0,w,h,pixels);
   for(let y=0;y<h;y++)flipped.set(pixels.subarray((h-y-1)*w*4,(h-y)*w*4),y*w*4);
   canvas.getContext('2d')?.putImageData(new ImageData(flipped,w,h),0,0);renderer.setRenderTarget(null);previewDirty=false;
  };
  let sectionBounds:Bounds3={min:[-.5,0,-.5],max:[.5,1.8,.5]};
  const sectionBox=(s:SceneState)=>{
   sectionBounds=sectionSetBounds(atlas,s);
   return new T.Box3(new T.Vector3().fromArray(sectionBounds.min),new T.Vector3().fromArray(sectionBounds.max));
  };
  const contextUniform={value:1},skinUniform={value:.1},sectionActiveUniform={value:0};
  const textures=new Map<string,T.Texture>(),textureLoader=new T.TextureLoader();
  const textureFor=(url:string,normal=false)=>{
   const key=`${normal?'normal':'color'}:${url}`;
   let texture=textures.get(key);
   if(!texture){
    texture=textureLoader.load(assetUrl(url),()=>{if(!disposed)dirty=true;},undefined,()=>{if(!disposed)onError(`Could not load local anatomy texture: ${url}`);});
    texture.colorSpace=normal?T.NoColorSpace:T.SRGBColorSpace;texture.wrapS=T.RepeatWrapping;texture.wrapT=T.RepeatWrapping;
    if(url.startsWith('/local-models/reference/'))texture.flipY=false;
    texture.anisotropy=Math.min(4,renderer.capabilities.getMaxAnisotropy());
    textures.set(key,texture);
   }
   return texture;
  };
  type InnerFace='skin'|'fascia'|'serosa';
  const innerFaceFor=(part:typeof atlas.parts[number]):InnerFace|undefined=>{
   if(!referenceModel)return;
   // The derived mesh contains an outer and an inner shell. Topographic
   // region overlays are exterior-only and must retain their source surface.
   if(part.name.startsWith('Body surface (derived)'))return 'skin';
   if(/\b(pleura|peritone\w*|pericardi\w*|omentum|mesenter\w*|mesocolon|serosa)\b/i.test(part.name))return 'serosa';
   if(part.system==='fascia'||/\b(fascia|aponeurosis|capsule|membrane)\b/i.test(part.name))return 'fascia';
  };
  const materialFor=(system:string,ghost=false,source?:SceneMaterial,skinSensitive=system==='integumentary'||system==='regions'||system==='cell-boundary',innerFace?:InnerFace)=>{
   const color=source?.color?`#${source.color.map(value=>Math.round(Math.max(0,Math.min(255,value))).toString(16).padStart(2,'0')).join('')}`:system==='tendon'?'#d6c8b0':system==='cartilage'?'#afc0cb':SYSTEMS.find(s=>s.id===system)?.color??'#aebbb8';
   const m=new T.MeshStandardMaterial({color,map:source?.map?textureFor(source.map):null,normalMap:source?.normalMap?textureFor(source.normalMap,true):null,normalScale:new T.Vector2(source?.normalScale??1,source?.normalScale??1),vertexColors:!!source?.vertexColors,metalness:referenceModel?Math.min(.08,source?.metalness??.02):source?.metalness??.02,roughness:referenceModel?Math.max(.63,source?.roughness??.76):source?.roughness??.76,side:T.DoubleSide,transparent:ghost,opacity:1,depthWrite:!ghost,clippingPlanes:clipPlanes});
   if(innerFace&&!ghost)m.alphaHash=true;
   m.customProgramCacheKey=()=>system+':'+ghost+':'+skinSensitive+':'+!!source?.proceduralMuscle+':'+innerFace;
   m.onBeforeCompile=shader=>{
    shader.uniforms.contextOpacity=contextUniform;shader.uniforms.skinOpacity=skinUniform;shader.uniforms.partState={value:partTexture};shader.uniforms.selectionState={value:selectionTexture};shader.uniforms.rotationState={value:rotationTexture};shader.uniforms.stateWidth={value:width};shader.uniforms.selectionRotation=rotationUniform;shader.uniforms.selectionPivot=pivotUniform;
    shader.vertexShader='attribute float partIndex; uniform sampler2D partState; uniform sampler2D selectionState; uniform sampler2D rotationState; uniform float stateWidth; uniform vec4 selectionRotation; uniform vec3 selectionPivot; varying float partVisible; varying float partSelected; vec3 rotateSelected(vec3 v){return v+2.0*cross(selectionRotation.xyz,cross(selectionRotation.xyz,v)+selectionRotation.w*v);}\n'+shader.vertexShader;
    shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nif (texture2D(rotationState,vec2((partIndex+0.5)/stateWidth,0.5)).r>0.5) objectNormal=rotateSelected(objectNormal);');
    shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvec2 stateUv = vec2((partIndex + 0.5) / stateWidth, 0.5); vec4 state = texture2D(partState, stateUv); transformed += state.xyz; partVisible = state.w; partSelected = texture2D(selectionState, stateUv).r; if(texture2D(rotationState,stateUv).r>0.5) transformed=selectionPivot+rotateSelected(transformed-selectionPivot);');
    shader.fragmentShader='uniform float contextOpacity; uniform float skinOpacity; varying float partVisible; varying float partSelected;\n'+shader.fragmentShader;
    shader.fragmentShader=shader.fragmentShader.replace('#include <clipping_planes_fragment>','#include <clipping_planes_fragment>\nif (partVisible < 0.5) discard;');
    shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>',`#include <color_fragment>
     float tissueAlpha = (partSelected > 0.5 ? 1.0 : contextOpacity) * ${skinSensitive?'skinOpacity':'1.0'} * ${(source?.opacity??1).toFixed(6)};
     if (${ghost?'tissueAlpha >= 0.999 || tissueAlpha < 0.001':'tissueAlpha < 0.999'+(innerFace?' && sectionActive < 0.5':'')}) discard;
     diffuseColor.a = tissueAlpha;`);
    if(innerFace){
     shader.uniforms.sectionActive=sectionActiveUniform;
     shader.vertexShader='attribute float innerSide; varying float vInnerSide; varying vec3 innerPosition; varying vec3 innerNormal;\n'+shader.vertexShader;
     shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvInnerSide=innerSide;innerPosition=position;innerNormal=normal;');
     shader.fragmentShader=`uniform float sectionActive; varying float vInnerSide; varying vec3 innerPosition; varying vec3 innerNormal;
      float innerHash(vec2 p){p=mod(p,vec2(11.0));return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}
      float innerCell(vec2 p){vec2 cell=floor(p),f=fract(p);float nearest=2.0;
       for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++){vec2 n=vec2(float(x),float(y)),q=cell+n;
        vec2 center=n+vec2(innerHash(q),innerHash(q+vec2(19.0,7.0)))*.55+.225;
        nearest=min(nearest,length(center-f));}return nearest;}
      float innerFiber(vec2 p){return .5+.5*sin(6.2831853*(p.x*7.0+.12*sin(6.2831853*p.y)));}
      vec3 innerFaceColor(vec3 p,vec3 n){vec3 w=pow(abs(normalize(n)),vec3(4.0));w/=max(w.x+w.y+w.z,1e-5);
       float cell=innerCell(p.zy*120.0)*w.x+innerCell(p.xz*120.0)*w.y+innerCell(p.xy*120.0)*w.z;
       float fiber=innerFiber(p.zy*27.0)*w.x+innerFiber(p.xz*27.0)*w.y+innerFiber(p.xy*27.0)*w.z;
       ${innerFace==='skin'?'float lobe=smoothstep(.29,.51,cell);return mix(vec3(.71,.51,.47),vec3(.91,.73,.59),lobe)*(.96+.04*fiber);':innerFace==='fascia'?'return mix(vec3(.55,.48,.43),vec3(.83,.77,.67),fiber)*(.91+.09*cell);':'return mix(vec3(.62,.40,.38),vec3(.82,.64,.58),.65*cell+.35*fiber);'} }
     `+shader.fragmentShader;
     if(ghost)shader.fragmentShader=shader.fragmentShader.replace('#include <clipping_planes_fragment>','#include <clipping_planes_fragment>\nif(sectionActive>0.5) discard;');
     shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\nif(sectionActive>0.5&&vInnerSide>0.5)diffuseColor.rgb=innerFaceColor(innerPosition,innerNormal);');
    }
    // Selection changes opacity only; tinting diffuseColor washes out source textures.
    if(source?.proceduralMuscle){
     shader.uniforms.muscleSurface={value:textureFor('/models/o3m-image-20.jpg')};
     shader.vertexShader='varying vec3 musclePosition; varying vec3 muscleNormal;\n'+shader.vertexShader;
     shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nmuscleNormal=objectNormal;');
     shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nmusclePosition=position;');
     shader.fragmentShader=`uniform sampler2D muscleSurface;
      varying vec3 musclePosition; varying vec3 muscleNormal;
      vec2 muscleUv(vec2 p){return vec2(.22)+vec2(.43,.48)*abs(fract(p*5.5)*2.0-1.0);}
      vec3 muscleColor(){
       vec3 weight=pow(abs(normalize(muscleNormal)),vec3(4.0));weight/=max(weight.x+weight.y+weight.z,1e-5);
       return texture2D(muscleSurface,muscleUv(musclePosition.zy)).rgb*weight.x+
        texture2D(muscleSurface,muscleUv(musclePosition.xz)).rgb*weight.y+
        texture2D(muscleSurface,muscleUv(musclePosition.xy)).rgb*weight.z;
      }
     `+shader.fragmentShader;
     shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\ndiffuseColor.rgb*=mix(vec3(1.0),clamp(muscleColor()*1.65,vec3(.55),vec3(1.25)),.68);');
    }
   };materials.push(m);return m;
  };
  const materialSpecs=new Map<string,{system:string;source?:SceneMaterial;skinSensitive?:boolean;innerFace?:InnerFace}>([...SYSTEMS.map(s=>s.id),'tendon','cartilage'].map(id=>[id,{system:id}]));
  materialSpecs.set('integumentary:internal',{system:'integumentary',skinSensitive:false});
  for(const part of atlas.parts)if(part.material&&atlas.materials?.[part.material]){
   const source=atlas.materials[part.material];
   const internal=part.system==='integumentary'&&!isSkinPart(part);
   const innerFace=innerFaceFor(part);
   materialSpecs.set(`source:${part.system}:${part.material}${internal?':internal':''}${innerFace?`:${innerFace}`:''}`,{system:part.system,skinSensitive:isSkinPart(part),innerFace,source:part.id.startsWith('ZA:')&&part.system==='muscular'&&!source.map?{...source,proceduralMuscle:true}:source});
  }
  const materialKeyFor=(part:typeof atlas.parts[number])=>{
   const innerFace=innerFaceFor(part);
   return part.material&&atlas.materials?.[part.material]?`source:${part.system}:${part.material}${part.system==='integumentary'&&!isSkinPart(part)?':internal':''}${innerFace?`:${innerFace}`:''}`:part.tissue==='tendon'||part.tissue==='cartilage'?part.tissue:part.system==='integumentary'&&!isSkinPart(part)?'integumentary:internal':part.system;
  };
  const mats=new Map([...materialSpecs].map(([key,{system,source,skinSensitive,innerFace}])=>[key,materialFor(system,false,source,skinSensitive,innerFace)])),ghostMats=new Map([...materialSpecs].map(([key,{system,source,skinSensitive,innerFace}])=>[key,materialFor(system,true,source,skinSensitive,innerFace)]));
  // Tissue cut colors follow anatomy-atlas.brianp.chatgpt.site's CAPTISSUE palette. Its source
  // textures live in GLBs; our Z-Anatomy meshes have no UV maps, so use a fine
  // procedural cut-face grain instead of projecting unrelated surface textures.
  const capPalette:Record<string,number>={muscular:0xB0605A,cardiac:0xAE3239,skeletal:0xCBA286,arterial:0xC8393F,venous:0x466FBA,nervous:0xE7C55C,lymphatic:0x7FB894,connective:0xB7B4AB,fascia:0xBEA99C,cartilage:0x86AEBD,tendon:0xC3B79B,digestive:0xBE8763,respiratory:0xBE8763,urinary:0xBE8763,endocrine:0xBE8763,integumentary:0xBBAA9F,sensory:0xB58BC4,attachments:0x8E5049};
  const capMaterials=new Map<string,T.MeshStandardMaterial>(),capMeshes:T.Mesh[]=[];
  let capGeneration=0,capTimer=0,capKey='';
  const clearCaps=()=>{for(const mesh of capMeshes){scene.remove(mesh);mesh.geometry.dispose();}capMeshes.length=0;for(const material of capMaterials.values())material.dispose();capMaterials.clear();dirty=true;};
  const capMaterialFor=(part:typeof atlas.parts[number],cutIndex:number,opacity:number,cortex=false)=>{
   const role=cortex?'skeletal-cortex':part.system==='connective'?/cartilage/i.test(part.material??part.name)?'cartilage':/tendon/i.test(part.material??part.name)?'tendon':'connective':part.system;
   const profile=sectionCapProfile(part),priority=profile.thinShell||profile.hollowWall?0:role==='skeletal-cortex'?4:role==='skeletal'||role==='cartilage'||role==='arterial'||role==='venous'||role==='nervous'?3:role==='muscular'||role==='tendon'?1:2;
   const key=`${role}:${cutIndex}:${priority}:${opacity.toFixed(4)}`;let material=capMaterials.get(key);
   if(!material){material=new T.MeshStandardMaterial({color:new T.Color(role==='skeletal'?0xC6A989:role==='skeletal-cortex'?0xE1D3BA:capPalette[role]??0x8B9099).multiplyScalar(role==='skeletal'?1:.82),roughness:.92,metalness:0,side:T.DoubleSide,transparent:opacity<.999,opacity,depthWrite:opacity>=.999,polygonOffset:true,polygonOffsetFactor:-1,polygonOffsetUnits:-2-priority*2,clippingPlanes:[clipPlanes[1-cutIndex]]});
    const strength=role==='muscular'||role==='cardiac'?.5:role==='skeletal-cortex'||role==='nervous'?.18:role==='skeletal'?.04:.24;
    material.customProgramCacheKey=()=>`section-cap-${key}`;
    material.onBeforeCompile=shader=>{
     if(role==='skeletal')shader.uniforms.boneTile={value:textureFor('/models/cancellous-bone-seamless.png')};
     shader.vertexShader='attribute vec2 cutUv; varying vec2 vCutUv;\n'+shader.vertexShader;
     shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvCutUv=cutUv;');
     shader.fragmentShader=(role==='skeletal'?'uniform sampler2D boneTile;\n':'')+'varying vec2 vCutUv;\nfloat cutHash(vec2 p){return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}\nfloat cutNoise(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(cutHash(i),cutHash(i+vec2(1.0,0.0)),f.x),mix(cutHash(i+vec2(0.0,1.0)),cutHash(i+vec2(1.0,1.0)),f.x),f.y);}\n'+shader.fragmentShader;
     shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\nfloat tissueGrain=.65*cutNoise(vCutUv*260.0)+.35*cutNoise(vCutUv*920.0);diffuseColor.rgb*=1.0+'+strength.toFixed(3)+'*(tissueGrain-.5);'+(role==='skeletal'?'diffuseColor.rgb*=mix(vec3(1.0),texture2D(boneTile,vCutUv*26.0).rgb*1.4,.32);':''));
    };capMaterials.set(key,material);}
   return material;
  };
  const queueCaps=(s:SceneState)=>{
   if(referenceModel&&!ready)return;
   const cuts=enabledSections(s);
   const key=cuts.length?JSON.stringify([cuts,s.visible,s.hidden,s.depthHidden,s.region,s.isolate,s.selected,s.skinOpacity,s.contextOpacity,s.explode,loaded]):'';
   if(key===capKey)return;capKey=key;
   const generation=++capGeneration;clearTimeout(capTimer);clearCaps();
   if(!cuts.length){stagedCapsReady=true;return;}
   const candidates=cuts.flatMap(({index:cutIndex})=>{
    const plane=clipPlanes[cutIndex].clone(),offset=-plane.constant,normal=plane.normal.toArray() as Point3;
    return atlas.parts.map((part,index)=>({part,index,cutIndex,plane})).filter(({part,index})=>{
     if(!pickers[index]||!sectionPartVisible(part,s)||sectionCapOpacity(part,s,atlas.materials?.[part.material??'']?.opacity??1)<.001)return false;
     if(s.explode>.001)return true;
     const range=sectionRange({min:part.bounds[0] as Point3,max:part.bounds[1] as Point3},normal);
     return range.min<=offset+.0002&&range.max>=offset-.0002;
    });
   });
   const jobs=new Map<string,typeof candidates>();
   for(const candidate of candidates){
    const group=candidate.part.id.startsWith('REF:')&&candidate.part.system==='skeletal'&&s.explode<.001?candidate.part.conceptId:candidate.part.id;
    const key=`${candidate.cutIndex}:${group}`;const job=jobs.get(key);if(job)job.push(candidate);else jobs.set(key,[candidate]);
   }
   const cutJobs=[...jobs.values()];
   capTimer=window.setTimeout(()=>{let cursor=0;const step=()=>{
    if(disposed||generation!==capGeneration)return;
    const deadline=performance.now()+12;
    while(cursor<cutJobs.length&&performance.now()<deadline){const job=cutJobs[cursor++],{part,index,cutIndex,plane}=job[0],picker=pickers[index];if(!picker)continue;
     const fragments=job.length>1?job.map(({index})=>{const geometry=pickers[index]!.geometry,partGeometry=new T.BufferGeometry();partGeometry.setAttribute('position',geometry.getAttribute('position'));partGeometry.setIndex(geometry.getIndex());return partGeometry;}):[];
     const joined=fragments.length?mergeGeometries(fragments):null;for(const fragment of fragments)fragment.dispose();
     const profile=sectionCapProfile(part);
     const geometry=sectionCapGeometry(joined??picker.geometry,plane,picker.matrixWorld,{...profile,outlineWidth:part.system==='skeletal'?.003:undefined,closureWidth:part.system==='skeletal'?.012:part.system==='muscular'?.016:undefined,closureFraction:part.system==='skeletal'||part.system==='muscular'?.3:undefined});joined?.dispose();if(!geometry)continue;
     const opacity=sectionCapOpacity(part,s,atlas.materials?.[part.material??'']?.opacity??1);
     const material=capMaterialFor(part,cutIndex,opacity),mesh=new T.Mesh(geometry,part.system==='skeletal'?[material,capMaterialFor(part,cutIndex,opacity,true)]:material);mesh.renderOrder=3;mesh.frustumCulled=false;mesh.visible=stagedCapsReady;mesh.userData.partIndex=index;mesh.userData.cutIndex=cutIndex;scene.add(mesh);capMeshes.push(mesh);
    }
    dirty=true;if(cursor<cutJobs.length)capTimer=window.setTimeout(step,0);else if(!stagedCapsReady){stagedCapsReady=true;for(const mesh of capMeshes)mesh.visible=true;lastState=null;dirty=true;}
   };step();},0);
  };
  const labelLayer=document.createElementNS('http://www.w3.org/2000/svg','svg');labelLayer.classList.add('anatomy-labels');labelLayer.setAttribute('aria-hidden','true');el.appendChild(labelLayer);
  const labelAnchors=new Map<number,T.Vector3>();
  const cameraValues=()=>[...camera.position.toArray(),...controls.target.toArray(),camera.view?.enabled?camera.view.offsetX/camera.view.fullWidth:0,camera.view?.enabled?camera.view.offsetY/camera.view.fullHeight:0];
  controls.addEventListener('end',()=>cameraCallback.current?.(cameraValues()));
  let loaded=0;
  const chunkParts=Array.from({length:atlas.chunks.length},()=>[] as number[]);
  atlas.parts.forEach((part,index)=>chunkParts[part.chunk].push(index));
  const initiallySelected=new Set(latest.current.selected);
  const chunkOrder=chunkParts.map((_,index)=>index).sort((a,b)=>{
   const selected=(index:number)=>chunkParts[index].some(i=>initiallySelected.has(atlas.parts[i].id));
   const surface=(index:number)=>chunkParts[index].filter(i=>isSkinPart(atlas.parts[i])).length;
   const distal=(index:number)=>chunkParts[index].reduce((sum,i)=>{const c=centers[i];return sum+Math.hypot(c.x,(c.y-.9)*.55);},0)/Math.max(1,chunkParts[index].length);
   return Number(selected(b))-Number(selected(a))||surface(b)-surface(a)||distal(b)-distal(a)||a-b;
  });
  const loadChunk=async(ci:number)=>{
   const chunk=atlas.chunks[ci],compressed=!!chunk.gzip&&typeof DecompressionStream!=='undefined';const response=await fetch(assetUrl(compressed?chunk.gzip!:chunk.url),{signal:abort.signal});const buffer=await decodeModelResponse(response,chunk.bytes,compressed);if(disposed)return;
   const groups=new Map<string,T.BufferGeometry[]>();
   chunkParts[ci].forEach(i=>{const p=atlas.parts[i];
    const g=new T.BufferGeometry();g.setAttribute('position',new T.BufferAttribute(new Float32Array(buffer,p.positions,p.vertexCount*3),3));
    // GPU normalized signed-short normals keep the complete atlas compact in memory.
    g.setAttribute('normal',new T.BufferAttribute(new Int16Array(buffer,p.normals,p.vertexCount*3),3,true));
    if(p.uvs!==undefined)g.setAttribute('uv',new T.BufferAttribute(new Float32Array(buffer,p.uvs,p.vertexCount*2),2));
    if(p.colors!==undefined)g.setAttribute('color',new T.BufferAttribute(new Uint8Array(buffer,p.colors,p.vertexCount*3),3,true));
    g.setIndex(new T.BufferAttribute(new Uint32Array(buffer,p.indices,p.indexCount),1));
    if(p.sourceOffset)g.translate(p.sourceOffset[0],p.sourceOffset[1],p.sourceOffset[2]);
    if(innerFaceFor(p)){
     const position=g.getAttribute('position'),normal=g.getAttribute('normal'),innerSide=new Uint8Array(p.vertexCount);
     const centerX=(p.bounds[0][0]+p.bounds[1][0])*.5,centerZ=(p.bounds[0][2]+p.bounds[1][2])*.5;
     for(let v=0;v<p.vertexCount;v++){
      const dx=position.getX(v)-centerX,dz=position.getZ(v)-centerZ;
      innerSide[v]=normal.getX(v)*dx+normal.getZ(v)*dz<-.002?255:0;
     }
     g.setAttribute('innerSide',new T.BufferAttribute(innerSide,1,true));
    }
    g.boundingBox=bounds[i].clone();g.computeBoundingSphere();const pick=new T.Mesh(g,new T.MeshBasicMaterial({side:T.DoubleSide}));materials.push(pick.material);pick.matrixAutoUpdate=false;pickers[i]=pick;geometries.push(g);
    g.setAttribute('partIndex',new T.BufferAttribute(new Float32Array(p.vertexCount).fill(i),1));
    const materialKey=materialKeyFor(p);const list=groups.get(materialKey)??[];list.push(g);groups.set(materialKey,list);
   });
   groups.forEach((gs,key)=>{const geometry=mergeGeometries(gs,false);if(!geometry)throw new Error('Could not assemble anatomy geometry.');geometries.push(geometry);const mesh=new T.Mesh(geometry,mats.get(key));mesh.frustumCulled=false;mesh.renderOrder=2;scene.add(mesh);const ghost=new T.Mesh(geometry,ghostMats.get(key));ghost.frustumCulled=false;ghost.renderOrder=1;scene.add(ghost);
    let material=previewMaterials.get(key);
    if(!material){const spec=materialSpecs.get(key),source=spec?.source,system=spec?.system??key;
     const color=source?.color?new T.Color().setRGB(source.color[0]/255,source.color[1]/255,source.color[2]/255):new T.Color(SYSTEMS.find(item=>item.id===system)?.color??(atlas.scope==='cell'?'#e6a8cf':'#d8c3b2'));
     material=new T.MeshStandardMaterial({color,map:source?.map?textureFor(source.map):null,normalMap:source?.normalMap?textureFor(source.normalMap,true):null,normalScale:new T.Vector2(source?.normalScale??1,source?.normalScale??1),vertexColors:!!source?.vertexColors,roughness:source?.roughness??.8,metalness:source?.metalness??0,side:T.DoubleSide,transparent:false,opacity:1,depthWrite:true});
     material.customProgramCacheKey=()=>`section-preview-${key}`;
     material.onBeforeCompile=shader=>{
      shader.uniforms.previewState={value:previewTexture};shader.uniforms.previewWidth={value:width};
      shader.vertexShader='attribute float partIndex; uniform sampler2D previewState; uniform float previewWidth; varying float previewShown;\n'+shader.vertexShader;
      shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\npreviewShown=texture2D(previewState,vec2((partIndex+0.5)/previewWidth,0.5)).r;');
      shader.fragmentShader='varying float previewShown;\n'+shader.fragmentShader;
      shader.fragmentShader=shader.fragmentShader.replace('#include <clipping_planes_fragment>','#include <clipping_planes_fragment>\nif(previewShown<0.5) discard;');
     };
     previewMaterials.set(key,material);materials.push(material);}
    const previewMesh=new T.Mesh(geometry,material);previewMesh.frustumCulled=false;previewMesh.renderOrder=1;previewScene.add(previewMesh);
   });
   lastState=null;previewDirty=true;loaded++;onProgress(Math.round(loaded/atlas.chunks.length*100));dirty=true;
  };
  (async()=>{try{let cursor=0;await Promise.all(Array.from({length:3},async()=>{while(cursor<chunkOrder.length){const i=chunkOrder[cursor++];await loadChunk(i);}}));if(!disposed){ready=true;lastState=null;dirty=true;}}catch(e){if(!disposed)onError(e instanceof Error?e.message:'Could not load the anatomy.');}})();
  const directionFor=(view:string)=>view==='front'?new T.Vector3(0,.02,1):view==='back'?new T.Vector3(0,.02,-1):view==='side'?new T.Vector3(1,.02,0):view==='right'?new T.Vector3(-1,.02,0):view==='superior'?new T.Vector3(0,1,.001):view==='inferior'?new T.Vector3(0,-1,.001):new T.Vector3(.35,.06,1).normalize();
  const viewArea=()=>{
   const w=el.clientWidth,h=el.clientHeight,mobile=w<768,landscape=w>h&&h<=600;
   let left=20,right=w-70,top=90,bottom=h-24;
   if(mobile||landscape){bottom=h-150;right=w-60;}
   for(const selector of ['.layers-panel','.study-panel','.section-panel','.detail-sheet']){
    const panel=document.querySelector(selector);if(!panel||!panel.getClientRects().length)continue;
    const r=panel.getBoundingClientRect();
    if(mobile&&landscape){if(r.left<w/2)left=Math.max(left,r.right+12);else right=Math.min(right,r.left-12);}
    else if(mobile){bottom=Math.min(bottom,r.top-16);}
    else if(r.left<w/2){left=Math.max(left,r.right+24);}else{right=Math.min(right,r.left-24);}
   }
   return {w,h,left,right,top,bottom};
  };
  const fit=(view:string,extent=0)=>{
   const area=viewArea();camera.setViewOffset(area.w,area.h,area.w/2-(area.left+area.right)/2,area.h/2-(area.top+area.bottom)/2,area.w,area.h);
   const aspect=camera.aspect,mobile=el.clientWidth<768,normalDistance=mobile?Math.max(3.5,1.85*el.clientHeight/Math.max(160,el.clientHeight-400)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))):4;
   const reservedHeight=mobile?400:270;const availableAspect=Math.max(.35,(el.clientWidth-(mobile?40:340))/Math.max(160,el.clientHeight-reservedHeight));const atlasDistance=Math.max(packingHeight,packingWidth/availableAspect)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))*(el.clientHeight/Math.max(160,el.clientHeight-reservedHeight))*1.08;
   const distance=T.MathUtils.lerp(normalDistance,Math.max(.2,atlasDistance),extent);
   const direction=directionFor(view);
   controls.target.set(0,.85,0);const pivot=latest.current.isolate||amount<.05?selectionCenter(atlas.parts,latest.current.selected,data):null;if(pivot)controls.target.copy(pivot);let fittedDistance=distance;if((atlas.scope==='embryo'||latest.current.region&&latest.current.region!=='all')&&!pivot){const box=new T.Box3();atlas.parts.forEach((p,i)=>{if(!isSkinPart(p)&&(enabledSections(latest.current).length?sectionPartVisible(p,latest.current):partVisible(p,latest.current)))box.union(bounds[i]);});if(!box.isEmpty()){box.getCenter(controls.target);const size=box.getSize(new T.Vector3());const {w,h,left,right,top,bottom}=viewArea();camera.setViewOffset(w,h,w/2-(left+right)/2,h/2-(top+bottom)/2,w,h);fittedDistance=Math.max(size.y*h/Math.max(100,bottom-top),size.x*w/Math.max(100,right-left)/camera.aspect,size.z)/(2*Math.tan(T.MathUtils.degToRad(camera.fov/2)))*1.1;}}controls.maxDistance=Math.max(fittedDistance*1.25,.2);camera.position.copy(controls.target).addScaledVector(direction,fittedDistance);controls.update();if(enabledSections(latest.current).length){lastSectionOrientation='';lastSectionPoint=undefined;}dirty=true;
  };
  const resize=()=>{lastCamera=undefined;lastFocus=-1;layoutKey='';lastIsolate='';lastSelection='';lastState=null;lastSectionOrientation='';lastSectionPoint=undefined;renderer.setPixelRatio(Math.min(devicePixelRatio,el.clientWidth<768||el.clientHeight<600?1.5:2));camera.aspect=el.clientWidth/el.clientHeight;camera.updateProjectionMatrix();renderer.setSize(el.clientWidth,el.clientHeight);fit(amount>.04?'front':latest.current.view,Math.min(1,amount*6));};const observer=new ResizeObserver(resize);observer.observe(el);
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
    atlas.parts.forEach((part,i)=>{if(!selected.has(part.id)||!pickers[i])return;const mesh=pickers[i]!;const hit=raycaster.intersectObject(mesh,false).find(hit=>keptBySections(hit.point));if(hit)targetDistance=Math.min(targetDistance,hit.distance);});
    raycaster.far=Math.max(0,targetDistance-.001);
    atlas.parts.forEach((part,i)=>{const mesh=pickers[i];if(!mesh||selected.has(part.id)||!partVisible(part,inspectionState)||isSkinPart(part)&&(s.skinOpacity??0)<.01)return;if(!rotationPartIndices.has(i)){worldBox.copy(bounds[i]).translate(mesh.position);if(!raycaster.ray.intersectsBox(worldBox))return;}const hit=raycaster.intersectObject(mesh,false).find(hit=>keptBySections(hit.point));if(hit)candidates.set(i,Math.min(candidates.get(i)??Infinity,hit.distance));});
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
  const solidContext=()=> (latest.current.skinOpacity??.1)<.95&&atlas.parts.some((p,i)=>!isSkinPart(p)&&data[i*4+3]>.5);
  const canPick=(i:number,hasSolid:boolean)=>{
   const part=atlas.parts[i],selected=latest.current.selected.includes(part.id);
   if(data[i*4+3]<.5||hasSolid&&isSkinPart(part)&&!selected)return false;
   const alpha=selected?1:contextUniform.value*(isSkinPart(part)?skinUniform.value:1)*(atlas.materials?.[part.material??'']?.opacity??1);
   return alpha>.001;
  };
  const pickAt=(clientX:number,clientY:number,radius=16)=>{
   const rect=renderer.domElement.getBoundingClientRect(),hasSolid=solidContext();
   pointer.set((clientX-rect.left)/rect.width*2-1,-(clientY-rect.top)/rect.height*2+1);
   raycaster.setFromCamera(pointer,camera);
   // Covering scans may leave a short ray range; a pointer pick always checks
   // the full visible depth so hover and click resolve the same structure.
   raycaster.near=0;raycaster.far=Infinity;
   let nearest=Infinity,found=-1;
   // The visible cut face is its own mesh; raycasting only the source surface
   // misses the interior of a section and cannot identify the tissue there.
   for(const mesh of capMeshes){
    const index=mesh.userData.partIndex as number,cutIndex=mesh.userData.cutIndex as number;
    if(!mesh.visible||!canPick(index,hasSolid))continue;
    const hit=raycaster.intersectObject(mesh,false).find(candidate=>clipPlanes.every((plane,i)=>i===cutIndex||plane.distanceToPoint(candidate.point)>=-.00001));
    if(hit&&hit.distance<nearest){nearest=hit.distance;found=index;}
   }
   pickers.forEach((mesh,i)=>{
    if(!mesh||!canPick(i,hasSolid))return;
    if(referenceModel&&enabledSections(latest.current).length&&atlas.parts[i].name.startsWith('Body surface (derived)'))return;
    if(!rotationPartIndices.has(i)){
     worldBox.copy(bounds[i]).translate(mesh.position);
     if(!raycaster.ray.intersectBox(worldBox,hitPoint))return;
    }
    const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));
    if(hit&&hit.distance<nearest){nearest=hit.distance;found=i;}
   });
   if(found<0&&amount>.45&&!enabledSections(latest.current).length)found=findTarget(clientX-rect.left,clientY-rect.top,radius,i=>canPick(i,hasSolid));
   return found;
  };
  const down=(e:PointerEvent)=>{const rotating=!sectionToolRef.current&&amount>.04&&!latest.current.isolate&&rotationPartIndices.size>0&&e.button===0;pressPoints.set(e.pointerId,{x:e.clientX,y:e.clientY,lastX:e.clientX,lastY:e.clientY,anchored:false,rotating});hover.hidden=true;tap.down(e.pointerId,e.clientX,e.clientY,e.pointerType==='touch'?16:10);if(rotating){e.preventDefault();e.stopImmediatePropagation();renderer.domElement.setPointerCapture(e.pointerId);}};
  const move=(e:PointerEvent)=>{
   tap.move(e.pointerId,e.clientX,e.clientY);
   const press=pressPoints.get(e.pointerId);
   if(press?.rotating){
    e.preventDefault();e.stopImmediatePropagation();
    const dx=e.clientX-press.lastX,dy=e.clientY-press.lastY;press.lastX=e.clientX;press.lastY=e.clientY;
    if(dx||dy){const up=new T.Vector3(0,1,0).applyQuaternion(camera.quaternion),right=new T.Vector3(1,0,0).applyQuaternion(camera.quaternion);selectedRotation.premultiply(new T.Quaternion().setFromAxisAngle(up,dx*.009)).premultiply(new T.Quaternion().setFromAxisAngle(right,dy*.009)).normalize();syncSelectedRotation();}
    hover.hidden=true;renderer.domElement.style.cursor='grabbing';return;
   }
   if(e.buttons&&press&&!press.anchored&&Math.hypot(e.clientX-press.x,e.clientY-press.y)>8){press.anchored=true;anchorSelection();}
   if(e.buttons||!ready||sectionToolRef.current&&!enabledSections(latest.current).length||e.pointerType==='touch'){hover.hidden=true;renderer.domElement.style.cursor=e.buttons?'grabbing':'default';return;}
   const rect=el.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top,index=pickAt(e.clientX,e.clientY);
   hover.hidden=index<0;renderer.domElement.style.cursor='default';
   if(index>=0){hoverLabel.textContent=structureName(atlas.parts[index].name);hover.style.left=`${Math.max(8,Math.min(x+14,el.clientWidth-260))}px`;hover.style.top=`${Math.max(8,Math.min(y+18,el.clientHeight-55))}px`;}
  };
  const cancel=(e:PointerEvent)=>{pressPoints.delete(e.pointerId);tap.cancel(e.pointerId);hover.hidden=true;renderer.domElement.style.cursor='default';};
  const leave=()=>{hover.hidden=true;renderer.domElement.style.cursor='default';};
  const up=(e:PointerEvent)=>{
   const press=pressPoints.get(e.pointerId);pressPoints.delete(e.pointerId);if(press?.rotating){e.preventDefault();e.stopImmediatePropagation();}const validTap=tap.up(e.pointerId,e.clientX,e.clientY);if(!validTap||!ready)return;const rect=renderer.domElement.getBoundingClientRect();pointer.set((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1);raycaster.setFromCamera(pointer,camera);
   // Covering-tissue scans shorten raycaster.far; each click must trace the
   // entire visible anatomy, including structures beyond the last scan point.
   raycaster.near=0;raycaster.far=Infinity;
   if(sectionToolRef.current&&!enabledSections(latest.current).length){let nearest=Infinity,point:T.Vector3|undefined;pickers.forEach((mesh,i)=>{if(!mesh||!sectionPartVisible(atlas.parts[i],latest.current))return;const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));if(hit&&hit.distance<nearest){nearest=hit.distance;point=hit.point;}});if(point){const section=latest.current.section??{enabled:false,axis:'axial',position:.38,flip:true};sectionPositionRef.current?.(sectionFraction(sectionBounds,section,point.toArray() as Point3));}return;}
   const found=pickAt(e.clientX,e.clientY,e.pointerType==='touch'?24:16);
   if(found>=0){hover.hidden=true;renderer.domElement.style.cursor='default';select.current(atlas.parts[found].id,e.shiftKey);}
   else if(sectionToolRef.current){let nearest=Infinity,point:T.Vector3|undefined;pickers.forEach((mesh,i)=>{if(!mesh||!sectionPartVisible(atlas.parts[i],latest.current))return;const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));if(hit&&hit.distance<nearest){nearest=hit.distance;point=hit.point;}});if(point){const section=latest.current.section??{enabled:false,axis:'axial',position:.38,flip:true};sectionPositionRef.current?.(sectionFraction(sectionBounds,section,point.toArray() as Point3));}}
  };
  renderer.domElement.addEventListener('pointerdown',down,true);renderer.domElement.addEventListener('pointermove',move);renderer.domElement.addEventListener('pointerup',up);renderer.domElement.addEventListener('pointercancel',cancel);renderer.domElement.addEventListener('pointerleave',leave);
  const sectionWheel=(event:WheelEvent)=>{const section=latest.current.section;if(!sectionToolRef.current||!section?.enabled||event.ctrlKey||event.metaKey)return;event.preventDefault();event.stopImmediatePropagation();sectionPositionRef.current?.(Math.max(0,Math.min(1,section.position+Math.sign(event.deltaY)*Math.min(.04,Math.abs(event.deltaY)/2500))));};
  renderer.domElement.addEventListener('wheel',sectionWheel,{capture:true,passive:false});
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
   let row=0;for(const id of s.selected.slice(0,20)){const i=atlas.parts.findIndex(p=>p.id===id),mesh=pickers[i];if(i<0||!mesh||data[i*4+3]<.5)continue;
    let anchor=labelAnchors.get(i);if(!anchor){const g=mesh.geometry,positions=g.getAttribute('position'),indices=g.index!;let best=Infinity;const center=centers[i];for(let j=0;j<indices.count;j+=Math.max(1,Math.floor(indices.count/1800))*3){const point=new T.Vector3();for(let k=0;k<3;k++)point.add(new T.Vector3().fromBufferAttribute(positions,indices.getX(Math.min(j+k,indices.count-1))));point.multiplyScalar(1/3);if(!keptBySections(point.clone().add(mesh.position)))continue;const d=point.distanceToSquared(center);if(d<best){best=d;anchor=point;}}if(!anchor)continue;labelAnchors.set(i,anchor);}
    const world=anchor.clone().applyMatrix4(mesh.matrixWorld);if(!keptBySections(world))continue;const point=world.clone().project(camera);
    if(point.z< -1||point.z>1||Math.abs(point.x)>1||Math.abs(point.y)>1)continue;
    const x=(point.x+1)*w/2,y=(1-point.y)*h/2,labelX=Math.min(w-165,Math.max(w<768?18:310,x+35)),labelY=Math.max(w<768?235:145,Math.min(h-200,y-25+row*23));row++;
    const make=(name:string,attrs:Record<string,string>)=>{const e=document.createElementNS('http://www.w3.org/2000/svg',name);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,v));labelLayer.appendChild(e);return e;};
    make('line',{x1:String(x),y1:String(y),x2:String(labelX),y2:String(labelY)});make('circle',{cx:String(x),cy:String(y),r:'3'});const text=make('text',{x:String(labelX+5),y:String(labelY-5)});text.textContent=structureName(atlas.parts[i].name);
   }
  };
  const clock=new T.Clock();let lastExtent=-1,lastCamera:number[]|undefined,lastFocus=0,lastAreaKey='',lastRotationKey='',lastSectionOrientation='',lastSectionPoint:T.Vector3|undefined,selectionStencil=false;
  const animate=()=>{
   if(disposed)return;frame=requestAnimationFrame(animate);const dt=Math.min(clock.getDelta(),.05),s=latest.current;attachPreview();
   const previewShell=(!ready||!stagedCapsReady)&&!s.selected.length&&(referenceModel||!enabledSections(s).length)&&atlas.parts.some(isSkinPart);
   const nextSectionActive=enabledSections(s).length?1:0;if(sectionActiveUniform.value!==nextSectionActive){sectionActiveUniform.value=nextSectionActive;dirty=true;}
   const nextContext=previewShell?1:s.selected.length?(s.contextOpacity??1):1,nextSkin=previewShell?1:(s.skinOpacity??.1);
   if(contextUniform.value!==nextContext||skinUniform.value!==nextSkin){contextUniform.value=nextContext;skinUniform.value=nextSkin;dirty=true;}
   const nextStencil=s.selected.length>0&&nextContext<.999;
   if(nextStencil!==selectionStencil){selectionStencil=nextStencil;
    for(const material of mats.values()){material.stencilWrite=nextStencil;material.stencilRef=1;material.stencilFunc=T.AlwaysStencilFunc;material.stencilZPass=T.ReplaceStencilOp;}
    for(const material of ghostMats.values()){material.stencilWrite=nextStencil;material.stencilRef=1;material.stencilFunc=T.NotEqualStencilFunc;material.stencilZPass=T.KeepStencilOp;}
    dirty=true;
   }
   const changed=lastState!==s,selectionCleared=!!lastState?.selected.length&&!s.selected.length;
   if(changed){const a=viewArea(),key=[a.w,a.h,a.left,a.right,a.top,a.bottom].join(':');if(key!==lastAreaKey){if(!selectionCleared)camera.setViewOffset(a.w,a.h,a.w/2-(a.left+a.right)/2,a.h/2-(a.top+a.bottom)/2,a.w,a.h);if(lastAreaKey&&s.focus&&!selectionCleared)lastFocus=-1;if(enabledSections(s).length){lastSectionOrientation='';lastSectionPoint=undefined;}lastAreaKey=key;dirty=true;}}

   if(changed){previewDirty=true;labelAnchors.clear();const box=sectionBox(s),section=s.section??{enabled:false,axis:'axial',position:.38,flip:true};
    clipPlanes.forEach(plane=>{plane.constant=10000;});
    for(const cut of enabledSections(s)){const equation=sectionPlaneEquation(sectionBounds,cut.section),plane=clipPlanes[cut.index];plane.normal.fromArray(equation.normal);plane.constant=equation.constant;}
    const point=sectionPoint(sectionBounds,section),normal=sectionNormal(section),quaternion=new T.Quaternion().setFromUnitVectors(new T.Vector3(0,0,1),new T.Vector3().fromArray(normal)),size=box.getSize(new T.Vector3());
    const span=Math.max(section.axis==='axial'?size.x:size.y,section.axis==='sagittal'?size.z:size.x,size.z,.2)*1.12;
    guideFill.visible=sectionToolRef.current&&!enabledSections(s).length;guideBorder.visible=guideFill.visible;guideFill.position.fromArray(point);guideBorder.position.fromArray(point);guideFill.quaternion.copy(quaternion);guideBorder.quaternion.copy(quaternion);guideFill.scale.set(span,span,1);guideBorder.scale.set(span,span,1);
    dirty=true;}
   const moving=Math.abs(amount-s.explode)>.0001;
   if(moving){amount=T.MathUtils.damp(amount,s.explode,8,dt);dirty=true;}
   if(changed||moving||lastExtent<0){
    const selection=new Set(s.selected),activeCut=enabledSections(s).length>0;
    const visibleParts=atlas.parts.filter(p=>activeCut?sectionPartVisible(p,s):partVisible(p,s));
    const nextLayoutKey=visibleParts.map(p=>p.id).join(',')+':'+camera.aspect.toFixed(3)+':'+hierarchyRef.current+':'+(hierarchyRef.current==='guest'?guestKeyRef.current:'');
    const layoutChanged=nextLayoutKey!==layoutKey;
    if(layoutChanged){explosionLayout=createHierarchicalExplosionLayout(atlas,anatomyEntries,new Set(visibleParts.map(part=>part.id)),hierarchyRef.current,camera.aspect,guestNodesRef.current);layoutKey=nextLayoutKey;stepsCallback.current?.(explosionLayout.steps);}
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
     const selected=selection.has(p.id),checked=sectionPartVisible(p,s);
     data.set([dx,dy,dz,(previewShell?isSkinPart(p)&&inRegion(p,s.region):activeCut?checked:partVisible(p,s))?1:0],i*4);selectedData[i*4]=selected?255:0;previewData[i*4]=(hasSurfacePreview?isSkinPart(p):checked)?255:0;
     const mesh=pickers[i];if(mesh){mesh.quaternion.identity();mesh.position.set(dx,dy,dz);mesh.updateMatrix();mesh.updateMatrixWorld(true);}
    });syncSelectedRotation();partTexture.needsUpdate=true;selectionTexture.needsUpdate=true;previewTexture.needsUpdate=true;rotationTexture.needsUpdate=true;lastState=s;lastExtent=amount;dirty=true;if(changed)queueCaps(s);
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
   const section=s.section?.enabled?s.section:enabledSections(s)[0]?.section,orientation=section?.enabled?[section.axis,section.azimuth??'',section.elevation??'',section.flip].join(':'):'';
   if(orientation){const point=new T.Vector3().fromArray(sectionPoint(sectionBounds,section!));
    if(orientation!==lastSectionOrientation&&!s.camera){
     const direction=new T.Vector3().fromArray(sectionNormal(section!)).multiplyScalar(section!.flip?-1:1);
     camera.up.copy(Math.abs(direction.y)>.9?new T.Vector3(0,0,-1):new T.Vector3(0,1,0));
     const right=new T.Vector3().crossVectors(camera.up,direction).normalize(),up=new T.Vector3().crossVectors(direction,right).normalize();
     const cutBox=new T.Box3(),cutNormal=sectionNormal(section!),cutOffset=sectionDot(point.toArray() as Point3,cutNormal);
     atlas.parts.forEach((part,i)=>{if(!sectionPartVisible(part,s))return;const partBounds={min:part.bounds[0] as Point3,max:part.bounds[1] as Point3},range=sectionRange(partBounds,cutNormal);if(range.min>cutOffset+.001||range.max<cutOffset-.001)return;
      if(clipPlanes.some(plane=>sectionRange(partBounds,plane.normal.toArray() as Point3).max+plane.constant<-.001))return;
      cutBox.union(bounds[i]);});
     const size=cutBox.isEmpty()?new T.Vector3().fromArray(sectionBounds.max).sub(new T.Vector3().fromArray(sectionBounds.min)):cutBox.getSize(new T.Vector3());
     const spanX=Math.abs(right.x)*size.x+Math.abs(right.y)*size.y+Math.abs(right.z)*size.z;
     const spanY=Math.abs(up.x)*size.x+Math.abs(up.y)*size.y+Math.abs(up.z)*size.z;
     const {w,h,left,right:areaRight,top,bottom}=viewArea(),tan=Math.tan(T.MathUtils.degToRad(camera.fov/2));
     const distance=Math.max(.12,spanX*w/Math.max(160,areaRight-left)/camera.aspect,spanY*h/Math.max(160,bottom-top))/(2*tan)*(referenceModel ? .9 : 1.12);
     controls.maxDistance=Math.max(controls.maxDistance,distance*1.5);
     controls.target.copy(point);camera.position.copy(point).addScaledVector(direction,distance);controls.update();dirty=true;
    }
    else if(lastSectionPoint&&orientation===lastSectionOrientation&&!point.equals(lastSectionPoint)){const delta=point.clone().sub(lastSectionPoint);controls.target.add(delta);camera.position.add(delta);controls.update();dirty=true;}
    lastSectionPoint=point;
   }else if(lastSectionOrientation){camera.up.set(0,1,0);if(!s.camera)fit(s.view,Math.min(1,amount*6));lastSectionPoint=undefined;dirty=true;}
   lastSectionOrientation=orientation;
   // Keep depth precision near the anatomy as the camera moves. A fixed tiny near
   // plane causes close skin and muscle surfaces to fight at whole-body distances.
   const near=Math.max(.00001,Math.min(.05,camera.position.distanceTo(controls.target)/500));
   if(Math.abs(camera.near-near)>near*.01){camera.near=near;camera.updateProjectionMatrix();dirty=true;}
   if(dirty){
    renderer.render(scene,camera);cameraCallback.current?.(cameraValues());drawLabels(s);targets=[];
    if(amount>.45){
     const hasSolid=(latest.current.skinOpacity??.1)<.95&&atlas.parts.some((p,i)=>!isSkinPart(p)&&data[i*4+3]>.5);
     atlas.parts.forEach((p,i)=>{
      if(data[i*4+3]<.5||(!s.selected.includes(p.id)&&contextUniform.value<.001)||(hasSolid&&isSkinPart(p)&&!s.selected.includes(p.id)))return;
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
   if(previewCanvas&&(previewDirty||changed))drawPreview(s);
   if(ready){const key=[s.selected.join(','),s.visible.join(','),s.region,s.skinOpacity,s.explode,...selectedRotation.toArray().map(v=>v.toFixed(2)),...camera.position.toArray().map(v=>v.toFixed(2)),...controls.target.toArray().map(v=>v.toFixed(2))].join(':');if((key!==lastCoveringKey||s.hidden!==lastCoveringHidden||s.depthHidden!==lastCoveringDepth)&&performance.now()-lastCoveringScan>250){scanCovering(s);lastCoveringKey=key;lastCoveringHidden=s.hidden;lastCoveringDepth=s.depthHidden;lastCoveringScan=performance.now();}}

  };animate();
  onCapture?.(async()=>{
   if(!ready||disposed)throw new Error('The 3D view is still loading.');
   renderer.render(scene,camera);
   const png=renderer.domElement.toDataURL('image/png');
   return (await fetch(png)).blob();
  });
  const contextLost=(e:Event)=>{e.preventDefault();onError('The 3D session was paused by your device. Reload to continue.');};renderer.domElement.addEventListener('webglcontextlost',contextLost);
  return()=>{onCapture?.(null);disposed=true;abort.abort();cancelAnimationFrame(frame);clearTimeout(capTimer);clearCaps();observer.disconnect();theme.removeEventListener('change',applyTheme);renderer.domElement.removeEventListener('wheel',sectionWheel,true);previewCanvas?.removeEventListener('pointerdown',previewDown);previewCanvas?.removeEventListener('pointermove',previewMove);previewCanvas?.removeEventListener('pointerup',previewUp);previewCanvas?.removeEventListener('pointercancel',previewUp);controls.dispose();geometries.forEach(g=>g.dispose());materials.forEach(m=>m.dispose());textures.forEach(t=>t.dispose());guideGeometry.dispose();guideBorder.geometry.dispose();guideBorder.material.dispose();guideFill.material.dispose();previewPlanes.forEach(({geometry,fill,edge})=>{geometry.dispose();edge.geometry.dispose();fill.material.dispose();edge.material.dispose();});previewTarget.dispose();scene.traverse(o=>{if(o instanceof T.Mesh&&!geometries.includes(o.geometry)&&o!==guideFill){o.geometry.dispose();const ms=Array.isArray(o.material)?o.material:[o.material];ms.forEach(m=>m.dispose());}});env.dispose();partTexture.dispose();selectionTexture.dispose();previewTexture.dispose();rotationTexture.dispose();hover.remove();labelLayer.remove();renderer.dispose();renderer.domElement.remove();};
 },[atlas]);
 return <div className="scene" ref={host}/>;
}
