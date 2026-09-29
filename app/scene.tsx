import {assetUrl} from './asset-url';
import {fetchAsset,loadModelBundle,rememberModelBundle} from './asset-cache';
import {inRegion,partVisible,sectionPartVisible} from './viewer-state';
import {useEffect,useMemo,useRef} from 'react';
import * as T from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import {RoomEnvironment} from 'three/examples/jsm/environments/RoomEnvironment.js';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';
import {createHierarchicalExplosionLayout,explosionBoundsAt,nextExplosionBounds,type ExplodeHierarchy,type HierarchicalExplosionLayout} from './hierarchical-explosion';
import {dampMotion,smoothStep} from './transition-motion';
import {hierarchyEntries} from './anatomy-hierarchy';
import {resolveGuestHierarchy,type GuestHierarchy} from './guest-hierarchy';
import {loadModelBuffer} from './model-download';
import {selectionCenter} from './camera-pivot';
import {orbitAroundPointer} from './pointer-orbit';
import {advanceControllerOrbit} from './xr-orbit';
import {labelAreaForPanels,layoutAnatomyLabels,type LabelAnchor,type LabelArea,type LabelPanel} from './label-layout';
import {displayLaterality,lateralityClass} from './laterality';
import {PointerTap} from './pointer-tap';
import {sectionDot,sectionFraction,sectionNormal,sectionPlaneEquation,sectionPoint,sectionRange,type Bounds3,type Point3} from './section-plane';
import {sectionCapGeometry} from './section-cap';
import {analyzeMeshTopology} from './mesh-topology';
import {sectionCapProfile} from './section-cap-profile';
import {sectionTissue} from './section-tissue';
import {createStencilCaps} from './stencil-caps';
import {sectionCapOpacity} from './section-opacity';
import {sectionSetBounds} from './section-set';
import {enabledSections} from './section-stack';
import {SYSTEMS,structureName,type Atlas,type SceneState,type SystemId} from './anatomy';
import {DEPTH_LAYERS,isSkinPart} from './depth-layers';
import {depthLayerOpacity,depthPosition,partLayerOpacity,setDepthPosition} from './depth-control';
import type {SceneCapture} from './export-tools';
import type {CaptureSize,SceneFrameCapture} from './page-capture';
import type {VrAnatomy} from './vr-anatomy';
import type {ScenePose} from './connect-state';
import {applyFrame,captureCameraPose,captureFrame,cameraPoseForFrame,decodeFrame,frameDistance,immersiveArea,interpolateCameraPose,preferredAnchor,screenAnchor,setFrameOffset,usableViewArea,zoomInToFill,zoomOutToFit,type CameraPose,type ViewArea,type ViewRect} from './view-framing';
import {validCamera} from '../shared/camera-frame.mjs';
interface Props {atlas:Atlas;vrModelUrl:string;state:SceneState;hierarchy:ExplodeHierarchy;guestHierarchy?:GuestHierarchy;remotePoseRef?:{current:ScenePose|null};remoteModel?:string;onPose?:(pose:Omit<ScenePose,'model'>)=>void;cameraTransition?:{id:number;from:number[]};sectionTool?:boolean;onSectionPosition?:(position:number)=>void;onSelect:(id:string,toggle?:boolean)=>void;onHidePart:(id:string,peel?:boolean)=>void;onCamera?:(camera:number[])=>void;onCovering?:(ids:string[])=>void;onExplosionSteps?:(steps:number)=>void;onCapture?:(capture:SceneCapture|null)=>void;onFrameCapture?:(capture:SceneFrameCapture|null)=>void;onProgress:(n:number)=>void;onError:(s:string)=>void}
type SceneMaterial={color:number[];map?:string;normalMap?:string;normalScale?:number;roughness?:number;metalness?:number;opacity?:number;vertexColors?:boolean};
export default function AnatomyScene({atlas,vrModelUrl,state,hierarchy,guestHierarchy,remotePoseRef,remoteModel,onPose,cameraTransition,sectionTool=false,onSectionPosition,onSelect,onHidePart,onCamera,onCovering,onExplosionSteps,onCapture,onFrameCapture,onProgress,onError}:Props){
 const guestNodes=useMemo(()=>guestHierarchy?resolveGuestHierarchy(atlas,guestHierarchy):undefined,[atlas,guestHierarchy]);
 const host=useRef<HTMLDivElement>(null),latest=useRef(state),transitionRef=useRef(cameraTransition),hierarchyRef=useRef(hierarchy),guestNodesRef=useRef(guestNodes),guestKeyRef=useRef(''),sectionToolRef=useRef(sectionTool),sectionPositionRef=useRef(onSectionPosition),select=useRef(onSelect),hidePart=useRef(onHidePart),cameraCallback=useRef(onCamera),coveringCallback=useRef(onCovering),stepsCallback=useRef(onExplosionSteps);
 latest.current=state;transitionRef.current=cameraTransition;hierarchyRef.current=hierarchy;guestNodesRef.current=guestNodes;guestKeyRef.current=guestHierarchy?`${guestHierarchy.id}:${guestHierarchy.revision??0}:${guestHierarchy.roots?.join(',')??guestHierarchy.source??''}`:'';sectionToolRef.current=sectionTool;sectionPositionRef.current=onSectionPosition;select.current=onSelect;hidePart.current=onHidePart;cameraCallback.current=onCamera;coveringCallback.current=onCovering;stepsCallback.current=onExplosionSteps;
 const remoteRef=useRef(remotePoseRef),remoteModelRef=useRef(remoteModel),poseCallback=useRef(onPose);remoteRef.current=remotePoseRef;remoteModelRef.current=remoteModel;poseCallback.current=onPose;
 useEffect(()=>{
  const el=host.current!;let disposed=false,frame=0,dirty=true,ready=false,lastView='',lastReset=-1,lastIsolate='',lastSelection='',lastRegion='',layoutKey='',amount=0;
  let lastState:SceneState|null=null,depthAmount:number|undefined,depthVelocity=0;
  let lastTransitionId=0;
  let cameraTween:{from:CameraPose;to:CameraPose;start:number;duration:number;source:'saved'|'action'}|null=null;
  const cancelCameraTween=()=>{if(cameraTween){cameraTween=null;controls.enableDamping=true;}};
  const abort=new AbortController();
  let renderer:T.WebGLRenderer;
  try{renderer=new T.WebGLRenderer({antialias:true,alpha:true,stencil:true,powerPreference:'high-performance'});}catch{onError('This browser could not start the 3D viewer. Please try a browser with WebGL enabled.');return;}
  const referenceModel=atlas.version==='Anatomy Atlas adapted GLB local import';
  const pairedSkinPart=(part:typeof atlas.parts[number])=>referenceModel&&part.name.startsWith('Body surface (derived)')||atlas.version==='BodyParts3D 4.0'&&part.id==='FJ2810';
  renderer.setPixelRatio(Math.min(devicePixelRatio,innerWidth<768?1.5:2));renderer.setClearColor('#f2f3f3');renderer.outputColorSpace=T.SRGBColorSpace;renderer.toneMapping=T.ACESFilmicToneMapping;renderer.toneMappingExposure=referenceModel ? .78 : .96;el.appendChild(renderer.domElement);
  renderer.domElement.setAttribute('aria-label',atlas.scope==='cell'?'Interactive human cell. Drag to orbit, pinch or scroll to zoom, tap a component to inspect it, Control-click to add or remove it, or Shift-click to hide it.':'Interactive human anatomy. Drag to orbit, pinch or scroll to zoom, tap a structure to inspect it, Control-click to add or remove it, or Shift-click to hide it.');
  renderer.domElement.setAttribute('data-scene-canvas','');
  let pageCaptureCanvas:HTMLCanvasElement|null=null,captureSize:CaptureSize|undefined;
  const pageFrameListeners=new Set<Parameters<SceneFrameCapture['subscribe']>[0]>();
  const syncCaptureResolution=()=>{const normal=Math.min(devicePixelRatio,el.clientWidth<768||el.clientHeight<600?1.5:2),ratio=captureSize?Math.min(captureSize.width/Math.max(1,el.clientWidth),captureSize.height/Math.max(1,el.clientHeight)):normal;if(renderer.getPixelRatio()!==ratio){renderer.setPixelRatio(ratio);dirty=true;return true;}return false;};
  const rememberPageFrame=()=>{if(!pageCaptureCanvas)return;const source=renderer.domElement;if(pageCaptureCanvas.width!==source.width||pageCaptureCanvas.height!==source.height){pageCaptureCanvas.width=source.width;pageCaptureCanvas.height=source.height;}pageCaptureCanvas.getContext('2d')?.drawImage(source,0,0);};
  const scene=new T.Scene(),camera=new T.PerspectiveCamera(34,1,.01,100),controls=new OrbitControls(camera,renderer.domElement);
  const desktopAnatomy=new T.Group();scene.add(desktopAnatomy);
  const questVr=/OculusBrowser|Quest 2/i.test(navigator.userAgent),vrPreview=import.meta.env.DEV&&new URLSearchParams(location.search).has('vr-preview');
  let vrAnatomy:VrAnatomy|null=null,vrPreparation:Promise<void>|null=null,xrSupported=false,xrArSupported=false,lightweightActive=false;
  const xrButton=document.createElement('button');xrButton.type='button';xrButton.className='enter-vr glass';xrButton.textContent='Enter VR';xrButton.setAttribute('aria-label','Enter immersive VR');
  let xrSession:XRSession|null=null,xrActiveMode:'immersive-vr'|'immersive-ar'|null=null,xrBaseSpace:XRReferenceSpace|null=null,xrLoopActive=false,xrZoom=1,xrAnchorMatrix:T.Matrix4|null=null;
  let xrOrbitHeld=false,xrLastLeftAim:{yaw:number;pitch:number}|null=null,xrOrbit={yaw:0,pitch:0};
  const xrRig=new T.Group(),xrRenderCamera=new T.PerspectiveCamera(),xrPresentationCamera=new T.PerspectiveCamera(),xrHeadCamera=new T.PerspectiveCamera();xrRig.add(xrRenderCamera);scene.add(xrRig);
  const xrArButton=document.createElement('button');xrArButton.type='button';xrArButton.className='enter-ar glass';xrArButton.textContent='Enter AR';xrArButton.setAttribute('aria-label','Enter immersive AR');
  const xrRays=[0,1].map(()=>{const ray=new T.Line(new T.BufferGeometry().setFromPoints([new T.Vector3(),new T.Vector3(0,0,-3)]),new T.LineBasicMaterial({color:0x52d5ca,transparent:true,opacity:.85,depthTest:false}));ray.visible=false;scene.add(ray);return ray;});
  const xrInfoCanvas=document.createElement('canvas');xrInfoCanvas.width=1024;xrInfoCanvas.height=256;
  const xrInfoTexture=new T.CanvasTexture(xrInfoCanvas),xrInfo=new T.Sprite(new T.SpriteMaterial({map:xrInfoTexture,transparent:true,depthTest:false}));xrInfo.scale.set(.8,.2,1);xrInfo.position.set(0,1.95,.08);xrInfo.visible=false;scene.add(xrInfo);
  let xrInfoText='';
  const updateXrInfo=(message:string)=>{if(message===xrInfoText)return;xrInfoText=message;const ctx=xrInfoCanvas.getContext('2d');if(!ctx)return;ctx.clearRect(0,0,1024,256);ctx.fillStyle='rgba(24,39,51,.9)';ctx.beginPath();ctx.roundRect(0,0,1024,256,28);ctx.fill();ctx.textAlign='center';ctx.fillStyle='#d9f5f1';ctx.font='bold 46px Arial';ctx.fillText(message.length>38?`${message.slice(0,35)}…`:message,512,116);ctx.fillStyle='#a8bbc0';ctx.font='28px Arial';ctx.fillText(questVr&&xrActiveMode==='immersive-vr'?'L trigger: peel   R trigger + left turn: orbit   L stick: zoom':'Trigger: select     Left stick: move closer or farther',512,186);xrInfoTexture.needsUpdate=true;};
  updateXrInfo('Human Atlas');
  camera.position.set(1.4,1.05,3.6);controls.target.set(0,.85,0);controls.zoomToCursor=true;controls.enableDamping=true;controls.dampingFactor=.085;controls.minDistance=.003;controls.maxDistance=6;controls.zoomSpeed=1.25;controls.maxPolarAngle=Math.PI-.001;controls.addEventListener('change',()=>{dirty=true;});
  const pmrem=new T.PMREMGenerator(renderer),room=new RoomEnvironment(),env=pmrem.fromScene(room,.04);scene.environment=env.texture;scene.environmentIntensity=referenceModel ? .7 : 1;room.dispose();pmrem.dispose();
  const hemisphere=new T.HemisphereLight(0xffffff,0x8c97a0,referenceModel ? .55 : .8);scene.add(hemisphere);
  const key=new T.DirectionalLight(0xfffaf4,referenceModel?1.65:2.05);key.position.set(-2,4,3);scene.add(key);
  const rim=new T.DirectionalLight(0xe9f0ff,referenceModel ? .65 : 1.1);rim.position.set(2,2,-3);scene.add(rim);
  const theme=window.matchMedia('(prefers-color-scheme: dark)');
  const applyTheme=()=>{renderer.setClearColor(theme.matches?'#141b23':'#f2f3f3',xrSession?.environmentBlendMode==='opaque'||!xrSession?1:0);dirty=true;};
  applyTheme();theme.addEventListener('change',applyTheme);
  const width=T.MathUtils.ceilPowerOfTwo(atlas.parts.length),data=new Float32Array(width*4),partTexture=new T.DataTexture(data,width,1,T.RGBAFormat,T.FloatType);partTexture.needsUpdate=true;
  const selectedData=new Uint8Array(width*4),selectionTexture=new T.DataTexture(selectedData,width,1);selectionTexture.needsUpdate=true;
  const previewData=new Uint8Array(width*4),previewTexture=new T.DataTexture(previewData,width,1);previewTexture.needsUpdate=true;
  const rotationData=new Uint8Array(width*4),rotationTexture=new T.DataTexture(rotationData,width,1);rotationTexture.needsUpdate=true;
  const selectedRotation=new T.Quaternion(),rotationUniform={value:new T.Vector4(0,0,0,1)},pivotUniform={value:new T.Vector3()},rotationPartIndices=new Set<number>(),partIndices=new Map(atlas.parts.map((part,i)=>[part.id,i]));
  const materials:T.Material[]=[],geometries:T.BufferGeometry[]=[],pickers:(T.Mesh|undefined)[]=[],centers=atlas.parts.map(p=>new T.Vector3().fromArray(p.bounds[0]).add(new T.Vector3().fromArray(p.bounds[1])).multiplyScalar(.5));
  const bounds=atlas.parts.map(p=>new T.Box3(new T.Vector3().fromArray(p.bounds[0]),new T.Vector3().fromArray(p.bounds[1])));
  const modelBounds=new T.Box3();atlas.parts.forEach((part,index)=>{if(!part.suppressed)modelBounds.union(bounds[index]);});
  let currentArea:ViewArea=usableViewArea([]),currentLabelArea:LabelArea=labelAreaForPanels(0,0,[]),layoutDirty=true,areaInitialized=false;
  let initialSelectionFitPending=!!(latest.current.isolate||latest.current.focus)&&!latest.current.camera;
  controls.addEventListener('start',()=>{initialSelectionFitPending=false;});
  const viewBounds=modelBounds.clone(),framingPoint=new T.Vector3();
  const anatomyEntries=hierarchyEntries(atlas);let explosionLayout:HierarchicalExplosionLayout|null=null;
  let packingWidth=1,packingHeight=1,explosionAspect=0,lastExplode=0,explosionVelocity=0,zoomVelocity=0;
  const explosionFitBounds=new T.Box3(),explosionVisibleBounds=new T.Box3();let explosionFitPending=false,implosionFillActive=false;
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
  const contextUniform={value:1},sectionActiveUniform={value:0};
  const textures=new Map<string,T.Texture>(),textureLoader=new T.ImageLoader(),textureUrls=new Set<string>();
  const textureFor=(url:string,normal=false)=>{
   const key=`${normal?'normal':'color'}:${url}`;
   let texture=textures.get(key);
   if(!texture){
    const pendingTexture=new T.Texture();texture=pendingTexture;
    void (async()=>{
     const response=await fetchAsset(assetUrl(url),{signal:abort.signal});if(!response.ok)throw new Error('Texture unavailable.');
     const blob=await response.blob();if(disposed)return;
     const objectUrl=URL.createObjectURL(blob);textureUrls.add(objectUrl);
     try{const image=await textureLoader.loadAsync(objectUrl);if(!disposed){pendingTexture.image=image;pendingTexture.needsUpdate=true;dirty=true;}}
     finally{URL.revokeObjectURL(objectUrl);textureUrls.delete(objectUrl);}
    })().catch(()=>{if(!disposed)onError(`Could not load local anatomy texture: ${url}`);});
    texture.colorSpace=normal?T.NoColorSpace:T.SRGBColorSpace;texture.wrapS=T.RepeatWrapping;texture.wrapT=T.RepeatWrapping;
    if(url.startsWith('/local-models/reference/'))texture.flipY=false;
    texture.anisotropy=Math.min(4,renderer.capabilities.getMaxAnisotropy());
    textures.set(key,texture);
   }
   return texture;
  };
  type InnerFace='skin'|'fascia'|'serosa';
  const innerFaceFor=(part:typeof atlas.parts[number]):InnerFace|undefined=>{
   if(pairedSkinPart(part))return 'skin';
   if(!referenceModel)return;
   // Topographic overlays are exterior-only and retain their source surface.
   if(/\b(pleura|peritone\w*|pericardi\w*|omentum|mesenter\w*|mesocolon|serosa)\b/i.test(part.name))return 'serosa';
   if(part.system==='fascia'||/\b(fascia|aponeurosis|capsule|membrane)\b/i.test(part.name))return 'fascia';
  };
  const materialFor=(system:string,ghost=false,source?:SceneMaterial,skinSensitive=system==='integumentary'||system==='regions'||system==='cell-boundary',innerFace?:InnerFace)=>{
   const color=source?.color?`#${source.color.map(value=>Math.round(Math.max(0,Math.min(255,value))).toString(16).padStart(2,'0')).join('')}`:system==='tendon'?'#d6c8b0':system==='cartilage'?'#afc0cb':SYSTEMS.find(s=>s.id===system)?.color??'#aebbb8';
   const m=new T.MeshStandardMaterial({color,map:source?.map?textureFor(source.map):null,normalMap:source?.normalMap?textureFor(source.normalMap,true):null,normalScale:new T.Vector2(source?.normalScale??1,source?.normalScale??1),vertexColors:!!source?.vertexColors,metalness:referenceModel?Math.min(.08,source?.metalness??.02):source?.metalness??.02,roughness:referenceModel?Math.max(.63,source?.roughness??.76):source?.roughness??.76,side:T.DoubleSide,transparent:ghost,opacity:1,depthWrite:!ghost,clippingPlanes:clipPlanes});
   const muscleSurface=system==='muscular'&&!source?.map;
   m.customProgramCacheKey=()=>system+':'+ghost+':'+skinSensitive+':'+muscleSurface+':'+innerFace;
   m.onBeforeCompile=shader=>{
    shader.uniforms.contextOpacity=contextUniform;shader.uniforms.partState={value:partTexture};shader.uniforms.selectionState={value:selectionTexture};shader.uniforms.rotationState={value:rotationTexture};shader.uniforms.stateWidth={value:width};shader.uniforms.selectionRotation=rotationUniform;shader.uniforms.selectionPivot=pivotUniform;
    shader.vertexShader='attribute float partIndex; uniform sampler2D partState; uniform sampler2D selectionState; uniform sampler2D rotationState; uniform float stateWidth; uniform vec4 selectionRotation; uniform vec3 selectionPivot; varying float partVisible; varying float partSelected; varying float partOpacity; vec3 rotateSelected(vec3 v){return v+2.0*cross(selectionRotation.xyz,cross(selectionRotation.xyz,v)+selectionRotation.w*v);}\n'+shader.vertexShader;
    shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nif (texture2D(rotationState,vec2((partIndex+0.5)/stateWidth,0.5)).r>0.5) objectNormal=rotateSelected(objectNormal);');
    shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvec2 stateUv = vec2((partIndex + 0.5) / stateWidth, 0.5); vec4 state = texture2D(partState, stateUv); transformed += state.xyz; partVisible = state.w; vec4 selection=texture2D(selectionState,stateUv); partSelected=selection.r; partOpacity=selection.g; if(texture2D(rotationState,stateUv).r>0.5) transformed=selectionPivot+rotateSelected(transformed-selectionPivot);');
    shader.fragmentShader='uniform float contextOpacity; varying float partVisible; varying float partSelected; varying float partOpacity;\n'+shader.fragmentShader;
    shader.fragmentShader=shader.fragmentShader.replace('#include <clipping_planes_fragment>','#include <clipping_planes_fragment>\nif (partVisible < 0.5) discard;');
    shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>',`#include <color_fragment>
     float tissueAlpha = (partSelected > 0.5 ? 1.0 : contextOpacity) * partOpacity * ${(source?.opacity??1).toFixed(6)};
     if (${ghost?'tissueAlpha >= 0.999 || tissueAlpha < 0.001':'tissueAlpha < 0.999'}) discard;
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
     shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\nif(sectionActive>0.5&&vInnerSide>0.5)diffuseColor.rgb=innerFaceColor(innerPosition,innerNormal);');
    }
    // Selection changes opacity only; tinting diffuseColor washes out source textures.
    if(muscleSurface){
     shader.uniforms.muscleSurface={value:textureFor('/models/muscle-surface-seamless.png')};
     shader.vertexShader='varying vec3 musclePosition; varying vec3 muscleNormal;\n'+shader.vertexShader;
     shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nmuscleNormal=objectNormal;');
     shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nmusclePosition=position;');
     shader.fragmentShader=`uniform sampler2D muscleSurface;
      varying vec3 musclePosition; varying vec3 muscleNormal;
      vec2 muscleUv(vec2 p){return p*36.0;}
      vec3 muscleColor(){
       vec3 weight=pow(abs(normalize(muscleNormal)),vec3(4.0));weight/=max(weight.x+weight.y+weight.z,1e-5);
       return texture2D(muscleSurface,muscleUv(musclePosition.zy)).rgb*weight.x+
        texture2D(muscleSurface,muscleUv(musclePosition.xz)).rgb*weight.y+
        texture2D(muscleSurface,muscleUv(musclePosition.xy)).rgb*weight.z;
      }
     `+shader.fragmentShader;
     shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\ndiffuseColor.rgb=mix(diffuseColor.rgb,muscleColor(),.83);');
    }
   };materials.push(m);return m;
  };
  const materialSpecs=new Map<string,{system:string;source?:SceneMaterial;skinSensitive?:boolean;innerFace?:InnerFace}>([...SYSTEMS.map(s=>s.id),'tendon','cartilage'].map(id=>[id,{system:id}]));
  materialSpecs.set('integumentary:internal',{system:'integumentary',skinSensitive:false});
  if(atlas.version==='BodyParts3D 4.0')materialSpecs.set('integumentary:paired-skin',{system:'integumentary',skinSensitive:true,innerFace:'skin'});
  const spleenColor:[number,number,number]=[203,168,94];
  for(const part of atlas.parts)if(part.material&&atlas.materials?.[part.material]){
   const source=atlas.materials[part.material];
   const internal=part.system==='integumentary'&&!isSkinPart(part);
   const innerFace=innerFaceFor(part);
   const spleen=/^spleen(?:\s|$)/i.test(part.name);
   materialSpecs.set(`source:${part.system}:${part.material}${internal?':internal':''}${innerFace?`:${innerFace}`:''}${spleen?':spleen':''}`,{system:part.system,skinSensitive:isSkinPart(part),innerFace,source:spleen?{...source,color:spleenColor,map:undefined,normalMap:undefined}:source});
  }
  const materialKeyFor=(part:typeof atlas.parts[number])=>{
   const innerFace=innerFaceFor(part);
   const spleen=/^spleen(?:\s|$)/i.test(part.name);
   if(pairedSkinPart(part)&&!part.material)return 'integumentary:paired-skin';
   return part.material&&atlas.materials?.[part.material]?`source:${part.system}:${part.material}${part.system==='integumentary'&&!isSkinPart(part)?':internal':''}${innerFace?`:${innerFace}`:''}${spleen?':spleen':''}`:part.tissue==='tendon'||part.tissue==='cartilage'?part.tissue:part.system==='integumentary'&&!isSkinPart(part)?'integumentary:internal':part.system;
  };
  const mats=new Map([...materialSpecs].map(([key,{system,source,skinSensitive,innerFace}])=>[key,materialFor(system,false,source,skinSensitive,innerFace)])),ghostMats=new Map([...materialSpecs].map(([key,{system,source,skinSensitive,innerFace}])=>[key,materialFor(system,true,source,skinSensitive,innerFace)]));
  // Cut patterns use plane coordinates, so source meshes need no UVs. Keep
  // them on cap geometry; a cut texture must never stain the exterior mesh.
  const capPalette:Record<string,number>={muscular:0xB0605A,cardiac:0xAE3239,skeletal:0xCBA286,arterial:0xC8393F,venous:0x466FBA,nervous:0xE7C55C,lymphatic:0x7FB894,connective:0xB7B4AB,fascia:0xBEA99C,serosa:0xBC8B82,adipose:0xD8BE91,cartilage:0x86AEBD,tendon:0xC3B79B,lung:0xD99DA3,liver:0xA85A62,spleen:0xCBA85E,cns:0xC4BBB2,organ:0xB98580,integumentary:0xBBAA9F,attachments:0x8E5049};
  const stencilCaps=createStencilCaps(desktopAnatomy,clipPlanes as [T.Plane,T.Plane],url=>textureFor(url));
  const capMaterials=new Map<string,T.MeshStandardMaterial>(),capMeshes:T.Mesh[]=[];
  const topologyCache=new WeakMap<T.BufferGeometry,ReturnType<typeof analyzeMeshTopology>>();
  const topologyFor=(geometry:T.BufferGeometry)=>{let value=topologyCache.get(geometry);if(!value){value=analyzeMeshTopology(geometry);topologyCache.set(geometry,value);}return value;};
  type CapCandidate={part:typeof atlas.parts[number];index:number;cutIndex:number;plane:T.Plane};
  type CapJob={key:string;job:CapCandidate[];signature:string};
  let capGeneration=0,capTimer=0,capKey='';
  const capQueue:CapJob[]=[],pendingCapSignatures=new Map<string,string>();
  const capSignatures=new Map<string,string>();
  const clearCaps=()=>{for(const mesh of capMeshes){desktopAnatomy.remove(mesh);mesh.geometry.dispose();}capMeshes.length=0;capSignatures.clear();for(const material of capMaterials.values())material.dispose();capMaterials.clear();renderer.domElement.dataset.sectionCapCount='0';dirty=true;};
  const capMaterialFor=(part:typeof atlas.parts[number],cutIndex:number,opacity:number,cortex=false)=>{
   const role=cortex?'skeletal-cortex':sectionTissue(part);
   const profile=sectionCapProfile(part),priority=profile.thinShell||profile.hollowWall?0:role==='skeletal-cortex'?4:role==='skeletal'||role==='cartilage'||role==='arterial'||role==='venous'||role==='nervous'?3:role==='muscular'||role==='tendon'?1:2;
   const key=`${role}:${cutIndex}:${priority}:${opacity.toFixed(4)}`;let material=capMaterials.get(key);
   const orientedMuscle=role==='muscular'||role==='cardiac';
   const tilePath:Record<string,string>={adipose:'/models/adipose-cut-v2.png',lung:'/models/lung-cut-seamless.png',liver:'/models/liver-cut-seamless.png',spleen:'/models/spleen-cut-seamless.png',cns:'/models/cns-cut-seamless.png',organ:'/models/organ-cut-seamless.png'};
   const tiled=orientedMuscle||!!tilePath[role];
   if(!material){material=new T.MeshStandardMaterial({color:new T.Color(tiled?0xFFFFFF:role==='skeletal'?0xC6A989:role==='skeletal-cortex'?0xE1D3BA:capPalette[role]??0x8B9099).multiplyScalar(tiled||role==='skeletal'?1:.82),roughness:.92,metalness:0,side:T.DoubleSide,transparent:opacity<.999,opacity,depthWrite:opacity>=.999,polygonOffset:true,polygonOffsetFactor:-1,polygonOffsetUnits:-1,clippingPlanes:[clipPlanes[1-cutIndex]]});
    const strength=role==='muscular'||role==='cardiac'?.5:role==='skeletal-cortex'||role==='nervous'?.18:role==='skeletal'?.04:.24;
    material.customProgramCacheKey=()=>`section-cap-${key}`;
    material.onBeforeCompile=shader=>{
     if(role==='skeletal')shader.uniforms.boneTile={value:textureFor('/models/cancellous-bone-seamless.png')};
     if(orientedMuscle){
      shader.uniforms.muscleAxial={value:textureFor('/models/muscle-axial-seamless.png')};
      shader.uniforms.muscleCoronal={value:textureFor('/models/muscle-coronal-seamless.png')};
      shader.uniforms.muscleSagittal={value:textureFor('/models/muscle-sagittal-seamless.png')};
      shader.uniforms.cutPlaneNormal={value:clipPlanes[cutIndex].normal};
     }else if(tilePath[role])shader.uniforms.cutTile={value:textureFor(tilePath[role])};
     shader.vertexShader='attribute vec2 cutUv; varying vec2 vCutUv;\n'+shader.vertexShader;
     shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvCutUv=cutUv;');
     shader.fragmentShader=(role==='skeletal'?'uniform sampler2D boneTile;\n':'')+(orientedMuscle?'uniform sampler2D muscleAxial; uniform sampler2D muscleCoronal; uniform sampler2D muscleSagittal; uniform vec3 cutPlaneNormal;\n':tilePath[role]?'uniform sampler2D cutTile;\n':'')+'varying vec2 vCutUv;\nfloat cutHash(vec2 p){p=mod(p,vec2(97.0));return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}\nfloat cutNoise(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(cutHash(i),cutHash(i+vec2(1.0,0.0)),f.x),mix(cutHash(i+vec2(0.0,1.0)),cutHash(i+vec2(1.0,1.0)),f.x),f.y);}\n'+shader.fragmentShader;
     const fibers='float cutFibers=.5+.5*sin(6.2831853*(vCutUv.x*160.0+.14*sin(6.2831853*vCutUv.y*20.0)));';
     const tissuePattern=role==='fascia'||role==='tendon'||role==='integumentary'?`${fibers}diffuseColor.rgb*=.95+.10*cutFibers;`:role==='serosa'?'diffuseColor.rgb*=.97+.06*cutNoise(vCutUv*180.0);':'';
     const tilePattern=orientedMuscle?'vec3 planeWeight=pow(abs(normalize(cutPlaneNormal)),vec3(4.0));planeWeight/=max(planeWeight.x+planeWeight.y+planeWeight.z,1e-5);vec3 tissueTile=texture2D(muscleSagittal,vCutUv*25.0).rgb*planeWeight.x+texture2D(muscleAxial,vCutUv*25.0).rgb*planeWeight.y+texture2D(muscleCoronal,vCutUv*25.0).rgb*planeWeight.z;diffuseColor.rgb=mix(diffuseColor.rgb,tissueTile,.88);':tilePath[role]?`diffuseColor.rgb=mix(diffuseColor.rgb,texture2D(cutTile,vCutUv*${role==='lung'?'20.0':'22.0'}).rgb,.88);`:'';
     shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\nfloat tissueGrain=.65*cutNoise(vCutUv*260.0)+.35*cutNoise(vCutUv*920.0);diffuseColor.rgb*=1.0+'+strength.toFixed(3)+'*(tissueGrain-.5);'+tissuePattern+tilePattern+(role==='skeletal'?'diffuseColor.rgb*=mix(vec3(1.0),texture2D(boneTile,vCutUv*26.0).rgb*1.4,.32);':''));
    };capMaterials.set(key,material);}
   return material;
  };
  const queueCaps=(s:SceneState)=>{
   if(lightweightActive)return;
   const cuts=enabledSections(s);
   const key=cuts.length?JSON.stringify([cuts,s.visible,s.hidden,s.depthHidden,s.depth,s.region,s.isolate,s.selected,s.skinOpacity,s.contextOpacity,s.explode]):'';
   const changed=key!==capKey;capKey=key;
   if(changed){++capGeneration;clearTimeout(capTimer);capTimer=0;capQueue.length=0;pendingCapSignatures.clear();clearCaps();}
   if(!cuts.length){stencilCaps.update([]);return;}
   const candidates:CapCandidate[]=cuts.flatMap(({index:cutIndex})=>{
    const plane=clipPlanes[cutIndex].clone(),offset=-plane.constant,normal=plane.normal.toArray() as Point3;
    return atlas.parts.map((part,index)=>({part,index,cutIndex,plane})).filter(({part,index})=>{
     if(!pickers[index]||!sectionPartVisible(part,s)||sectionCapOpacity(part,s,atlas.materials?.[part.material??'']?.opacity??1)<.001)return false;
     if(s.explode>.001)return true;
     const range=sectionRange({min:part.bounds[0] as Point3,max:part.bounds[1] as Point3},normal);
     return range.min<=offset+.0002&&range.max>=offset-.0002;
    });
   });
   const gpuRoles=new Set(['integumentary','muscular','cardiac','skeletal','cartilage','adipose','lung','liver','spleen','cns','organ','nervous','lymphatic','arterial','venous','tendon','connective']);
   const gpuCandidates=candidates.filter(({part})=>{
    const profile=sectionCapProfile(part);
    if(s.explode>.001||profile.thinShell&&!pairedSkinPart(part)||!gpuRoles.has(sectionTissue(part)))return false;
    return sectionCapOpacity(part,s,atlas.materials?.[part.material??'']?.opacity??1)>=.999;
   });
   const gpuStarted=performance.now();
   stencilCaps.update(gpuCandidates.map(({part,index,cutIndex})=>{const profile=sectionCapProfile(part);return {part,index,cutIndex,geometry:pickers[index]!.geometry,matrix:pickers[index]!.matrixWorld,mode:profile.hollowWall||pairedSkinPart(part)?'wall' as const:'solid' as const,wallWidth:profile.width};}));
   renderer.domElement.dataset.sectionGpuMs=(performance.now()-gpuStarted).toFixed(1);
   renderer.domElement.dataset.sectionCapCount=String(capMeshes.length+stencilCaps.activeCount());
   const gpuKeys=new Set(gpuCandidates.map(({index,cutIndex})=>`${cutIndex}:${index}`));
   const jobs=new Map<string,CapCandidate[]>();
   for(const candidate of candidates){
    if(gpuKeys.has(`${candidate.cutIndex}:${candidate.index}`))continue;
    // Reference GLBs often split one organ or muscle into material primitives.
    // Intersecting each primitive alone leaves artificial open cut edges.
    const combinable=['skeletal','muscular','cardiac','digestive','respiratory','urinary','reproductive'].includes(candidate.part.system);
    const group=candidate.part.id.startsWith('REF:')&&combinable&&s.explode<.001?candidate.part.conceptId:candidate.part.id;
    const key=`${candidate.cutIndex}:${candidate.part.system}:${group}`;const job=jobs.get(key);if(job)job.push(candidate);else jobs.set(key,[candidate]);
   }
   const selected=new Set(s.selected);
   const rank=(part:typeof atlas.parts[number])=>selected.has(part.id)?-1:isSkinPart(part)?0:part.system==='fascia'?1:part.system==='muscular'?2:part.system==='skeletal'?3:4;
   for(const [jobKey,job] of jobs){const signature=job.map(({index})=>index).join(',');
    if(capSignatures.get(jobKey)===signature||pendingCapSignatures.get(jobKey)===signature)continue;
    pendingCapSignatures.set(jobKey,signature);capQueue.push({key:jobKey,job,signature});
   }
   capQueue.sort((a,b)=>rank(a.job[0].part)-rank(b.job[0].part));
   renderer.domElement.dataset.sectionPendingCount=String(capQueue.length);
   renderer.domElement.dataset.sectionPendingParts=capQueue.map(({job})=>job[0].part.name).join(' | ');
   if(capTimer||!capQueue.length)return;
   const generation=capGeneration;
   const step=()=>{
    capTimer=0;if(disposed||generation!==capGeneration)return;
    const deadline=performance.now()+12;
    let processed=0;
    while(capQueue.length&&(!processed||performance.now()<deadline)){const {key:jobKey,job,signature}=capQueue.shift()!;processed++;
     if(pendingCapSignatures.get(jobKey)!==signature)continue;
     pendingCapSignatures.delete(jobKey);
     const {part,index,cutIndex,plane}=job[0],picker=pickers[index];if(!picker)continue;
     const previous=capMeshes.findIndex(mesh=>mesh.userData.capJobKey===jobKey);
     if(previous>=0){const mesh=capMeshes.splice(previous,1)[0];desktopAnatomy.remove(mesh);mesh.geometry.dispose();}
     const fragments=job.length>1?job.map(({index})=>{const geometry=pickers[index]!.geometry,partGeometry=new T.BufferGeometry();partGeometry.setAttribute('position',geometry.getAttribute('position'));partGeometry.setIndex(geometry.getIndex());return partGeometry;}):[];
     const joined=fragments.length?mergeGeometries(fragments):null;for(const fragment of fragments)fragment.dispose();
     const source=joined??picker.geometry,profile=sectionCapProfile(part,topologyFor(source));
     const geometry=sectionCapGeometry(source,plane,picker.matrixWorld,{...profile,outlineWidth:part.system==='skeletal'?.003:undefined,closureWidth:profile.closureWidth??(part.system==='skeletal'?.012:part.system==='muscular'?.016:undefined),closureFraction:profile.closureFraction??(part.system==='skeletal'||part.system==='muscular'?.3:undefined)});joined?.dispose();
     if(geometry){const opacity=sectionCapOpacity(part,s,atlas.materials?.[part.material??'']?.opacity??1);
      const material=capMaterialFor(part,cutIndex,opacity),mesh=new T.Mesh(geometry,part.system==='skeletal'?[material,capMaterialFor(part,cutIndex,opacity,true)]:material);mesh.renderOrder=0;mesh.frustumCulled=false;mesh.userData.partIndex=index;mesh.userData.cutIndex=cutIndex;mesh.userData.capJobKey=jobKey;desktopAnatomy.add(mesh);capMeshes.push(mesh);}
     capSignatures.set(jobKey,signature);
    }
    renderer.domElement.dataset.sectionCapCount=String(capMeshes.length+stencilCaps.activeCount());
    renderer.domElement.dataset.sectionPendingCount=String(capQueue.length);
    renderer.domElement.dataset.sectionPendingParts=capQueue.map(({job})=>job[0].part.name).join(' | ');
    dirty=true;if(capQueue.length)capTimer=window.setTimeout(step,0);
   };
   capTimer=window.setTimeout(step,0);
  };
  const labelLayer=document.createElementNS('http://www.w3.org/2000/svg','svg');labelLayer.classList.add('anatomy-labels');labelLayer.setAttribute('aria-hidden','true');el.appendChild(labelLayer);
  const labelMeasure=document.createElement('canvas');
  const labelAnchors=new Map<number,T.Vector3>();
  const cameraValues=()=>captureFrame(camera,controls.target,viewBounds,currentArea);
  controls.addEventListener('end',()=>cameraCallback.current?.(cameraValues()));
  let loaded=0,lastReportedProgress=0;
  const chunkParts=Array.from({length:atlas.chunks.length},()=>[] as number[]);
  atlas.parts.forEach((part,index)=>chunkParts[part.chunk].push(index));
  const chunkUrls=atlas.chunks.map(chunk=>assetUrl(chunk.gzip&&typeof DecompressionStream!=='undefined'?chunk.gzip:chunk.url));
  const chunkSizes=atlas.chunks.map(chunk=>chunk.bytes),chunkOffsets:number[]=[];let chunkTotal=0;
  for(const size of chunkSizes){chunkOffsets.push(chunkTotal);chunkTotal+=size;}
  const useBundle=!questVr&&chunkTotal<=256*1024*1024&&chunkOffsets.every(offset=>offset%4===0);
  const cachedBuffers:ArrayBuffer[]=[];let bundleBuffer:ArrayBuffer|undefined;
  const skinChunks=new Set(chunkParts.flatMap((parts,index)=>parts.some(i=>atlas.parts[i].system!=='cell-boundary'&&isSkinPart(atlas.parts[i]))?[index]:[]));
  let loadedSkinChunks=0,surfaceReady=skinChunks.size===0;
  const initiallySelected=new Set(latest.current.selected);
  const initialCuts=enabledSections(latest.current),initialBounds=initialCuts.length?sectionSetBounds(atlas,latest.current):null;
  const initialCutEquations=initialBounds?initialCuts.map(({section})=>sectionPlaneEquation(initialBounds,section)):[];
  const cutCounts=chunkParts.map(parts=>parts.filter(i=>{
   const part=atlas.parts[i];if(!sectionPartVisible(part,latest.current))return false;
   return initialCutEquations.some(equation=>{const normal=equation.normal as Point3,range=sectionRange({min:part.bounds[0] as Point3,max:part.bounds[1] as Point3},normal);return range.min<=-equation.constant+.0002&&range.max>=-equation.constant-.0002;});
  }).length);
  const hasPairedSurface=atlas.parts.some(pairedSkinPart);
  const chunkOrder=chunkParts.map((_,index)=>index).sort((a,b)=>{
   const selected=(index:number)=>chunkParts[index].some(i=>initiallySelected.has(atlas.parts[i].id));
   const surface=(index:number)=>Number(chunkParts[index].some(i=>hasPairedSurface?pairedSkinPart(atlas.parts[i]):isSkinPart(atlas.parts[i])));
   const distal=(index:number)=>chunkParts[index].reduce((sum,i)=>{const c=centers[i];return sum+Math.hypot(c.x,(c.y-.9)*.55);},0)/Math.max(1,chunkParts[index].length);
   // The visible outer boundary must arrive before deeper cut geometry, even
   // when it contains fewer pieces than a neurovascular chunk.
   return Number(selected(b))-Number(selected(a))||surface(b)-surface(a)||cutCounts[b]-cutCounts[a]||distal(b)-distal(a)||a-b;
  });
  const loadChunk=async(ci:number)=>{
   const chunk=atlas.chunks[ci],compressed=!!chunk.gzip&&typeof DecompressionStream!=='undefined';const buffer=bundleBuffer??await loadModelBuffer(chunkUrls[ci],chunk.bytes,compressed,abort.signal);if(disposed)return;
   if(!bundleBuffer)cachedBuffers[ci]=buffer;
   const offset=bundleBuffer?chunkOffsets[ci]:0;
   const groups=new Map<string,T.BufferGeometry[]>();
   chunkParts[ci].forEach(i=>{const p=atlas.parts[i];
    const positions=new Float32Array(buffer,offset+p.positions,p.vertexCount*3);
    const g=new T.BufferGeometry();g.setAttribute('position',new T.BufferAttribute(p.sourceOffset?positions.slice():positions,3));
    // GPU normalized signed-short normals keep the complete atlas compact in memory.
    g.setAttribute('normal',new T.BufferAttribute(new Int16Array(buffer,offset+p.normals,p.vertexCount*3),3,true));
    if(p.uvs!==undefined)g.setAttribute('uv',new T.BufferAttribute(new Float32Array(buffer,offset+p.uvs,p.vertexCount*2),2));
    if(p.colors!==undefined)g.setAttribute('color',new T.BufferAttribute(new Uint8Array(buffer,offset+p.colors,p.vertexCount*3),3,true));
    g.setIndex(new T.BufferAttribute(new Uint32Array(buffer,offset+p.indices,p.indexCount),1));
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
   groups.forEach((gs,key)=>{const geometry=mergeGeometries(gs,false);if(!geometry)throw new Error('Could not assemble anatomy geometry.');geometries.push(geometry);const mesh=new T.Mesh(geometry,mats.get(key));mesh.frustumCulled=false;mesh.renderOrder=2;desktopAnatomy.add(mesh);const ghost=new T.Mesh(geometry,ghostMats.get(key));ghost.frustumCulled=false;ghost.renderOrder=1;desktopAnatomy.add(ghost);
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
   if(skinChunks.has(ci)){loadedSkinChunks++;if(loadedSkinChunks===skinChunks.size){surfaceReady=true;lastState=null;}}
   previewDirty=true;loaded++;const nextProgress=Math.round(loaded/atlas.chunks.length*100);if(nextProgress===100||nextProgress-lastReportedProgress>=10){lastReportedProgress=nextProgress;onProgress(nextProgress);}dirty=true;
  };
  (async()=>{try{bundleBuffer=useBundle?await loadModelBundle(chunkUrls,chunkSizes):undefined;if(disposed)return;let cursor=0;await Promise.all(Array.from({length:questVr?3:8},async()=>{while(cursor<chunkOrder.length){const i=chunkOrder[cursor++];await loadChunk(i);}}));if(!disposed){ready=true;lastState=null;dirty=true;if(!bundleBuffer&&useBundle)void rememberModelBundle(chunkUrls,cachedBuffers);if(vrPreview||questVr&&xrSupported||xrArSupported)void prepareVr();}}catch(e){if(!disposed)onError(e instanceof Error?e.message:'Could not load the anatomy.');}})();
  const directionFor=(view:string)=>view==='front'?new T.Vector3(0,.02,1):view==='back'?new T.Vector3(0,.02,-1):view==='side'?new T.Vector3(1,.02,0):view==='right'?new T.Vector3(-1,.02,0):view==='superior'?new T.Vector3(0,1,.001):view==='inferior'?new T.Vector3(0,-1,.001):new T.Vector3(.35,.06,1).normalize();
  const panelSelectors='.identity,.top-actions,.advanced-panel,.layers-panel,.study-panel,.section-panel,.detail-sheet,.download-panel,.search-panel,.view-controls,.share-view,.guest-menu,.anatomy-choice-menu,.anatomy-search-results';
  const studio=el.closest('.studio')??el.parentElement!;
  const readLayoutAreas=()=>{
   const viewport=el.getBoundingClientRect(),panels:LabelPanel[]=[];
   for(const panel of studio.ownerDocument.querySelectorAll(panelSelectors)){
    if(!panel.getClientRects().length||getComputedStyle(panel).visibility==='hidden')continue;const r=panel.getBoundingClientRect();if(!r.width||!r.height)continue;
    const kind=panel.matches('.identity,.top-actions')?'top':panel.matches('.advanced-panel')?'bottom':'panel';
    panels.push({kind,left:r.left-viewport.left,right:r.right-viewport.left,top:r.top-viewport.top,bottom:r.bottom-viewport.top});
   }
   const rects:ViewRect[]=panels.map(panel=>({...panel,left:panel.left/viewport.width,right:panel.right/viewport.width,top:panel.top/viewport.height,bottom:panel.bottom/viewport.height}));
   return {view:usableViewArea(rects),labels:labelAreaForPanels(viewport.width,viewport.height,panels)};
  };
  const viewArea=()=>currentArea;
  const fit=(view:string,extent=0)=>{
   const area=viewArea(),direction=directionFor(view);setFrameOffset(camera,area);
   const packed=new T.Box3().setFromCenterAndSize(modelBounds.getCenter(new T.Vector3()),new T.Vector3(packingWidth,packingHeight,modelBounds.getSize(new T.Vector3()).z));
   const normalDistance=frameDistance(modelBounds,direction,camera.up,camera.aspect,area,camera.fov)*1.08,atlasDistance=frameDistance(packed,direction,camera.up,camera.aspect,area,camera.fov)*1.08;
   const distance=T.MathUtils.lerp(normalDistance,atlasDistance,extent);
   modelBounds.getCenter(controls.target);const pivot=latest.current.isolate||amount<.05?selectionCenter(atlas.parts,latest.current.selected,data):null;if(pivot)controls.target.copy(pivot);let fittedDistance=distance;
   if((atlas.scope==='embryo'||latest.current.region&&latest.current.region!=='all')&&!pivot){const box=new T.Box3();atlas.parts.forEach((p,i)=>{if(!isSkinPart(p)&&(enabledSections(latest.current).length?sectionPartVisible(p,latest.current):partVisible(p,latest.current)))box.union(bounds[i]);});if(!box.isEmpty()){box.getCenter(controls.target);fittedDistance=frameDistance(box,direction,camera.up,camera.aspect,area,camera.fov)*1.1;}}
   controls.maxDistance=Math.max(fittedDistance*1.25,.2);camera.position.copy(controls.target).addScaledVector(direction,fittedDistance);controls.update();if(enabledSections(latest.current).length){lastSectionOrientation='';lastSectionPoint=undefined;}dirty=true;
  };
  const markLayoutDirty=()=>{layoutDirty=true;};
  const observer=new ResizeObserver(markLayoutDirty);observer.observe(el);
  const observedPanels=new WeakSet<Element>();
  const observePanels=()=>{for(const panel of studio.ownerDocument.querySelectorAll(panelSelectors))if(!observedPanels.has(panel)){observedPanels.add(panel);observer.observe(panel);layoutObserver.observe(panel,{attributes:true,attributeFilter:['class','style','data-open','data-closed','data-starting-style','data-ending-style']});}};
  const layoutObserver=new MutationObserver(mutations=>{if(mutations.some(mutation=>mutation.type==='attributes'||[...mutation.addedNodes,...mutation.removedNodes].some(node=>node instanceof Element&&(node.matches(panelSelectors)||node.querySelector(panelSelectors))))){markLayoutDirty();observePanels();}});layoutObserver.observe(studio.ownerDocument.body,{childList:true,subtree:true});observePanels();
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
    atlas.parts.forEach((part,i)=>{const mesh=pickers[i];if(!mesh||selected.has(part.id)||!partVisible(part,inspectionState)||partLayerOpacity(part,s)<.01)return;if(!rotationPartIndices.has(i)){worldBox.copy(bounds[i]).translate(mesh.position);if(!raycaster.ray.intersectsBox(worldBox))return;}const hit=raycaster.intersectObject(mesh,false).find(hit=>keptBySections(hit.point));if(hit)candidates.set(i,Math.min(candidates.get(i)??Infinity,hit.distance));});
   }
   coveringCallback.current?.([...candidates].sort((a,b)=>a[1]-b[1]).map(([i])=>atlas.parts[i].id));
  };
  const selectionPivot=()=>amount>.05&&!latest.current.isolate?null:selectionCenter(atlas.parts,latest.current.selected,data);
  const syncSelectedRotation=()=>{
   const box=new T.Box3();for(const i of rotationPartIndices)box.union(bounds[i].clone().translate(new T.Vector3(data[i*4],data[i*4+1],data[i*4+2])));
   const pivot=box.isEmpty()?null:box.getCenter(new T.Vector3());pivotUniform.value.copy(pivot??new T.Vector3());rotationUniform.value.set(selectedRotation.x,selectedRotation.y,selectedRotation.z,selectedRotation.w);
   if(pivot)for(const i of rotationPartIndices){const mesh=pickers[i];if(!mesh)continue;
    mesh.quaternion.copy(selectedRotation);mesh.position.set(data[i*4],data[i*4+1],data[i*4+2]).sub(pivot).applyQuaternion(selectedRotation).add(pivot);mesh.updateMatrix();mesh.updateMatrixWorld(true);
   }
   dirty=true;
  };
  const pressPoints=new Map<number,{x:number;y:number;lastX:number;lastY:number;rotating:boolean;orbitCandidate:boolean;orbiting:boolean;pivot?:T.Vector3}>();
  const solidContext=()=> ((lastState??latest.current).depthHidden?.includes('skin')||depthLayerOpacity('skin',lastState??latest.current)<.95)&&atlas.parts.some((p,i)=>!isSkinPart(p)&&data[i*4+3]>.5);
  const canPick=(i:number,hasSolid:boolean)=>{
   const part=atlas.parts[i],selected=latest.current.selected.includes(part.id);
   const solidQuest=questVr&&xrActiveMode==='immersive-vr';
   if(data[i*4+3]<.5||!solidQuest&&hasSolid&&isSkinPart(part)&&!selected)return false;
   const alpha=(selected||solidQuest?1:contextUniform.value)*partLayerOpacity(part,lastState??latest.current)*(atlas.materials?.[part.material??'']?.opacity??1);
   return alpha>=(solidQuest?.5:.001);
  };
  const pickAt=(clientX:number,clientY:number,radius=16,pivotOut?:T.Vector3)=>{
   const rect=renderer.domElement.getBoundingClientRect(),hasSolid=solidContext();
   pivotOut?.set(NaN,NaN,NaN);
   pointer.set((clientX-rect.left)/rect.width*2-1,-(clientY-rect.top)/rect.height*2+1);
   raycaster.setFromCamera(pointer,camera);
   // Covering scans may leave a short ray range; a pointer pick always checks
   // the full visible depth so hover and click resolve the same structure.
   raycaster.near=0;raycaster.far=Infinity;
   let nearest=Infinity,found=-1;
   const stencilHit=stencilCaps.pick(raycaster.ray,index=>canPick(index,hasSolid));
   if(stencilHit){nearest=stencilHit.distance;found=stencilHit.index;if(pivotOut)raycaster.ray.at(nearest,pivotOut);}
   // The visible cut face is its own mesh; raycasting only the source surface
   // misses the interior of a section and cannot identify the tissue there.
   for(const mesh of capMeshes){
    const index=mesh.userData.partIndex as number,cutIndex=mesh.userData.cutIndex as number;
    if(!mesh.visible||!canPick(index,hasSolid))continue;
    const hit=raycaster.intersectObject(mesh,false).find(candidate=>clipPlanes.every((plane,i)=>i===cutIndex||plane.distanceToPoint(candidate.point)>=-.00001));
    if(hit&&hit.distance<nearest){nearest=hit.distance;found=index;pivotOut?.copy(hit.point);}
   }
   pickers.forEach((mesh,i)=>{
    if(!mesh||!canPick(i,hasSolid))return;
    if(enabledSections(latest.current).length&&pairedSkinPart(atlas.parts[i]))return;
    if(!rotationPartIndices.has(i)){
     worldBox.copy(bounds[i]).translate(mesh.position);
     if(!raycaster.ray.intersectBox(worldBox,hitPoint))return;
    }
    const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));
    if(hit&&hit.distance<nearest){nearest=hit.distance;found=i;pivotOut?.copy(hit.point);}
   });
   if(found<0&&amount>.45&&!enabledSections(latest.current).length)found=findTarget(clientX-rect.left,clientY-rect.top,radius,i=>canPick(i,hasSolid));
   return found;
  };
  const pointerPivot=(clientX:number,clientY:number)=>{
   const point=new T.Vector3();pickAt(clientX,clientY,16,point);
   if(Number.isFinite(point.x))return point;
   const rect=renderer.domElement.getBoundingClientRect();
   pointer.set((clientX-rect.left)/rect.width*2-1,-(clientY-rect.top)/rect.height*2+1);
   raycaster.setFromCamera(pointer,camera);
   const plane=new T.Plane().setFromNormalAndCoplanarPoint(camera.getWorldDirection(new T.Vector3()),controls.target);
   return raycaster.ray.intersectPlane(plane,point)??controls.target.clone();
  };
  const down=(e:PointerEvent)=>{if(remoteRef.current)return;cancelCameraTween();const toggling=e.ctrlKey||e.metaKey,rotating=!toggling&&!sectionToolRef.current&&amount>.04&&!latest.current.isolate&&rotationPartIndices.size>0&&e.button===0,orbitCandidate=!toggling&&!rotating&&!sectionToolRef.current&&(amount<=.04||latest.current.isolate)&&e.button===0&&e.pointerType!=='touch';pressPoints.set(e.pointerId,{x:e.clientX,y:e.clientY,lastX:e.clientX,lastY:e.clientY,rotating,orbitCandidate,orbiting:false});hover.hidden=true;tap.down(e.pointerId,e.clientX,e.clientY,e.pointerType==='touch'?16:10);if(rotating||toggling||orbitCandidate){e.preventDefault();e.stopImmediatePropagation();renderer.domElement.setPointerCapture(e.pointerId);}};
  const move=(e:PointerEvent)=>{
   if(remoteRef.current)return;
   tap.move(e.pointerId,e.clientX,e.clientY);
   const press=pressPoints.get(e.pointerId);
   if(press?.rotating){
    e.preventDefault();e.stopImmediatePropagation();
    const dx=e.clientX-press.lastX;press.lastX=e.clientX;
    if(dx){selectedRotation.premultiply(new T.Quaternion().setFromAxisAngle(new T.Vector3(0,1,0),dx*.009)).normalize();syncSelectedRotation();}
    hover.hidden=true;renderer.domElement.style.cursor='grabbing';return;
   }
   if(press?.orbitCandidate){
    if(!press.orbiting&&Math.hypot(e.clientX-press.x,e.clientY-press.y)>8){press.pivot=pointerPivot(press.x,press.y);press.orbiting=true;}
    if(press.orbiting&&press.pivot){
     e.preventDefault();e.stopImmediatePropagation();
     const dx=e.clientX-press.lastX,dy=e.clientY-press.lastY;press.lastX=e.clientX;press.lastY=e.clientY;
     if(dx||dy){const rect=renderer.domElement.getBoundingClientRect();pointer.set((press.x-rect.left)/rect.width*2-1,-(press.y-rect.top)/rect.height*2+1);orbitAroundPointer(camera,controls.target,press.pivot,pointer,dx,dy,rect.height,controls.minPolarAngle,controls.maxPolarAngle);controls.update();dirty=true;}
     hover.hidden=true;renderer.domElement.style.cursor='grabbing';return;
    }
   }
   if(e.buttons||!ready||sectionToolRef.current&&!enabledSections(latest.current).length||e.pointerType==='touch'){hover.hidden=true;renderer.domElement.style.cursor=e.buttons?'grabbing':'default';return;}
   const rect=el.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top,index=pickAt(e.clientX,e.clientY);
   hover.hidden=index<0;renderer.domElement.style.cursor='default';
   if(index>=0){const fullName=structureName(atlas.parts[index].name),display=displayLaterality(fullName);hoverLabel.textContent=display.label;hover.className=`part-hover ${lateralityClass(display.side)}`;hover.title=fullName;hover.setAttribute('aria-label',fullName);hover.style.left=`${Math.max(8,Math.min(x+14,el.clientWidth-260))}px`;hover.style.top=`${Math.max(8,Math.min(y+18,el.clientHeight-55))}px`;}
  };
  const cancel=(e:PointerEvent)=>{pressPoints.delete(e.pointerId);tap.cancel(e.pointerId);hover.hidden=true;renderer.domElement.style.cursor='default';};
  const leave=()=>{hover.hidden=true;renderer.domElement.style.cursor='default';};
  const up=(e:PointerEvent)=>{
   if(remoteRef.current)return;
   const press=pressPoints.get(e.pointerId);pressPoints.delete(e.pointerId);if(press?.rotating||press?.orbitCandidate){e.preventDefault();e.stopImmediatePropagation();}if(press?.orbiting){tap.cancel(e.pointerId);cameraCallback.current?.(cameraValues());return;}const validTap=tap.up(e.pointerId,e.clientX,e.clientY);if(!validTap||!ready)return;const rect=renderer.domElement.getBoundingClientRect();pointer.set((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1);raycaster.setFromCamera(pointer,camera);
   // Covering-tissue scans shorten raycaster.far; each click must trace the
   // entire visible anatomy, including structures beyond the last scan point.
   raycaster.near=0;raycaster.far=Infinity;
   if(sectionToolRef.current&&!enabledSections(latest.current).length&&!e.shiftKey){let nearest=Infinity,point:T.Vector3|undefined;pickers.forEach((mesh,i)=>{if(!mesh||!sectionPartVisible(atlas.parts[i],latest.current))return;const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));if(hit&&hit.distance<nearest){nearest=hit.distance;point=hit.point;}});if(point){const section=latest.current.section??{enabled:false,axis:'axial',position:.38,flip:true};sectionPositionRef.current?.(sectionFraction(sectionBounds,section,point.toArray() as Point3));}return;}
   const found=pickAt(e.clientX,e.clientY,e.pointerType==='touch'?24:16);
   if(found>=0){hover.hidden=true;renderer.domElement.style.cursor='default';const id=atlas.parts[found].id;if(e.shiftKey)hidePart.current(id);else select.current(id,e.ctrlKey||e.metaKey);}
   else if(sectionToolRef.current){let nearest=Infinity,point:T.Vector3|undefined;pickers.forEach((mesh,i)=>{if(!mesh||!sectionPartVisible(atlas.parts[i],latest.current))return;const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));if(hit&&hit.distance<nearest){nearest=hit.distance;point=hit.point;}});if(point){const section=latest.current.section??{enabled:false,axis:'axial',position:.38,flip:true};sectionPositionRef.current?.(sectionFraction(sectionBounds,section,point.toArray() as Point3));}}
  };
  const preventContextMenu=(e:MouseEvent)=>{if(e.ctrlKey)e.preventDefault();};
  renderer.domElement.addEventListener('contextmenu',preventContextMenu);
  renderer.domElement.addEventListener('pointerdown',down,true);renderer.domElement.addEventListener('pointermove',move);renderer.domElement.addEventListener('pointerup',up);renderer.domElement.addEventListener('pointercancel',cancel);renderer.domElement.addEventListener('pointerleave',leave);
  renderer.domElement.addEventListener('wheel',cancelCameraTween);
  const xrDirection=new T.Vector3(0,0,-1),xrOrigin=new T.Vector3(),xrQuaternion=new T.Quaternion(),xrAim=new T.Vector3(),xrControllerQuaternion=new T.Quaternion(),xrOrbitQuaternion=new T.Quaternion(),xrOrbitEuler=new T.Euler(0,0,0,'YXZ');
  const setLightweight=(enabled:boolean)=>{lightweightActive=enabled;desktopAnatomy.visible=!enabled;rim.visible=!enabled;if(vrAnatomy)vrAnatomy.mesh.visible=enabled;
   if(enabled){clearTimeout(capTimer);capTimer=0;capGeneration++;capQueue.length=0;pendingCapSignatures.clear();capKey='vr';guideFill.visible=false;guideBorder.visible=false;}
   lastState=null;dirty=true;
  };
  function prepareVr(){
   if(vrAnatomy)return Promise.resolve();if(vrPreparation)return vrPreparation;
   xrButton.disabled=true;xrButton.textContent='Preparing VR…';xrArButton.disabled=true;xrArButton.textContent='Preparing AR…';
   vrPreparation=(async()=>{
    const {loadVrAnatomy}=await import('./vr-anatomy');
    const prepared=await loadVrAnatomy(atlas,vrModelUrl,{rotationState:rotationTexture,rotation:rotationUniform,pivot:pivotUniform,planes:clipPlanes},abort.signal);
    if(disposed){prepared.dispose();return;}
    prepared.sync(latest.current,data,pickers,contextUniform.value);
    // Compile and upload before entering the headset's frame loop.
    const warmup=new T.Scene();warmup.add(hemisphere.clone(),key.clone(),prepared.mesh);prepared.mesh.visible=true;
    try{await renderer.compileAsync(warmup,camera);if(disposed){prepared.dispose();return;}
     const target=new T.WebGLRenderTarget(64,64),previous=renderer.getRenderTarget();target.texture.colorSpace=renderer.outputColorSpace;
     try{renderer.setRenderTarget(target);renderer.render(warmup,camera);}finally{renderer.setRenderTarget(previous);target.dispose();}
    }catch(error){prepared.dispose();throw error;}
    prepared.mesh.visible=false;scene.add(prepared.mesh);vrAnatomy=prepared;
    renderer.domElement.dataset.xrSourceTriangles=String(prepared.sourceTriangles);renderer.domElement.dataset.xrModelTriangles=String(prepared.triangles);
    xrButton.textContent='Enter VR';xrButton.disabled=false;xrArButton.textContent='Enter AR';xrArButton.disabled=false;if(vrPreview)setLightweight(true);
   })().catch(error=>{vrPreparation=null;if(!disposed){xrButton.textContent='Retry VR preparation';xrButton.disabled=false;xrArButton.textContent='Retry AR preparation';xrArButton.disabled=false;onError(error instanceof Error?error.message:'Could not prepare the lightweight VR model.');}});
   return vrPreparation;
  };
  const xrPose=(frame:XRFrame,source:XRInputSource)=>{const space=renderer.xr.getReferenceSpace(),pose=space&&frame.getPose(source.targetRaySpace,space);if(!pose)return false;const p=pose.transform.position,q=pose.transform.orientation;xrOrigin.set(p.x,p.y,p.z).applyMatrix4(xrRig.matrixWorld);xrQuaternion.set(q.x,q.y,q.z,q.w).premultiply(xrRig.quaternion);xrDirection.set(0,0,-1).applyQuaternion(xrQuaternion);return true;};
  const xrPick=(frame:XRFrame,source:XRInputSource)=>{if(!ready||questVr&&xrActiveMode==='immersive-vr'&&!xrAnchorMatrix||!xrPose(frame,source))return -1;raycaster.set(xrOrigin,xrDirection);raycaster.near=0;raycaster.far=xrPresentationCamera.position.distanceTo(controls.target)*2+modelBounds.getSize(new T.Vector3()).length();
   const hasSolid=solidContext();let closest=Infinity,found=-1;
   (lightweightActive&&vrAnatomy?vrAnatomy.pickers:pickers).forEach((mesh,i)=>{if(!mesh||!canPick(i,hasSolid)||enabledSections(latest.current).length&&pairedSkinPart(atlas.parts[i]))return;
    worldBox.copy(mesh.geometry.boundingBox!).applyMatrix4(mesh.matrixWorld);if(!raycaster.ray.intersectsBox(worldBox))return;
    const hit=raycaster.intersectObject(mesh,false).find(candidate=>keptBySections(candidate.point));if(hit&&hit.distance<closest){closest=hit.distance;found=i;}
   });
   return found;
  };
  const xrSelect=(event:XRInputSourceEvent)=>{if(remoteRef.current||questVr&&xrActiveMode==='immersive-vr')return;const found=xrPick(event.frame,event.inputSource);if(found>=0){select.current(atlas.parts[found].id);updateXrInfo(structureName(atlas.parts[found].name));}};
  const xrSelectStart=(event:XRInputSourceEvent)=>{if(remoteRef.current||!questVr||xrActiveMode!=='immersive-vr')return;
   if(event.inputSource.handedness==='right'){xrOrbitHeld=true;xrLastLeftAim=null;return;}
   if(event.inputSource.handedness!=='left'||xrOrbitHeld)return;
   const found=xrPick(event.frame,event.inputSource);if(found<0)return;
   const part=atlas.parts[found];hidePart.current(part.id,true);updateXrInfo(`Peeled ${structureName(part.name)}`);
  };
  const xrSelectEnd=(event:XRInputSourceEvent)=>{if(event.inputSource.handedness==='right'){xrOrbitHeld=false;xrLastLeftAim=null;}};
  const xrEnded=()=>{xrSession?.removeEventListener('select',xrSelect);xrSession?.removeEventListener('selectstart',xrSelectStart);xrSession?.removeEventListener('selectend',xrSelectEnd);xrSession=null;xrActiveMode=null;xrBaseSpace=null;xrAnchorMatrix=null;xrOrbitHeld=false;xrLastLeftAim=null;xrOrbit={yaw:0,pitch:0};xrRays.forEach(ray=>ray.visible=false);xrInfo.visible=false;labelLayer.style.display='';setLightweight(vrPreview);xrButton.textContent='Enter VR';xrButton.hidden=false;xrArButton.hidden=false;delete (studio as HTMLElement).dataset.xrMode;renderer.setAnimationLoop(null);renderer.xr.enabled=false;controls.enabled=true;applyTheme();markLayoutDirty();dirty=true;if(xrLoopActive){xrLoopActive=false;if(!disposed)frame=requestAnimationFrame(animate);}};
  renderer.xr.addEventListener('sessionend',xrEnded);
  const enterXr=async(mode:'immersive-vr'|'immersive-ar')=>{if(!navigator.xr||xrSession||!ready)return;if((questVr||mode==='immersive-ar')&&!vrAnatomy){await prepareVr();return;}xrButton.disabled=xrArButton.disabled=true;try{
   renderer.xr.enabled=true;renderer.xr.setReferenceSpaceType(mode==='immersive-ar'?'local':'local-floor');renderer.xr.setFramebufferScaleFactor(.8);
   const session=await navigator.xr.requestSession(mode,mode==='immersive-ar'?{requiredFeatures:['local'],optionalFeatures:['dom-overlay'],domOverlay:{root:studio}}:{requiredFeatures:['local-floor']});xrSession=session;xrActiveMode=mode;
   session.addEventListener('select',xrSelect);
   session.addEventListener('selectstart',xrSelectStart);session.addEventListener('selectend',xrSelectEnd);
   if(questVr||mode==='immersive-ar')setLightweight(true);
   await renderer.xr.setSession(session);
   renderer.xr.setFoveation(questVr?.5:1);
   xrBaseSpace=renderer.xr.getReferenceSpace();xrZoom=1;xrAnchorMatrix=null;xrOrbitHeld=false;xrLastLeftAim=null;xrOrbit={yaw:0,pitch:0};(studio as HTMLElement).dataset.xrMode=mode==='immersive-ar'?'ar':'vr';applyTheme();
   controls.enabled=false;labelLayer.style.display='none';hover.hidden=true;xrInfo.visible=true;xrInfoText='';updateXrInfo(latest.current.selected.length?structureName(atlas.parts[partIndices.get(latest.current.selected[0])??0].name):'Human Atlas');
   xrButton.hidden=xrArButton.hidden=true;cancelAnimationFrame(frame);xrLoopActive=true;renderer.setAnimationLoop(animate);
  }catch(error){renderer.xr.enabled=false;setLightweight(vrPreview);xrButton.hidden=xrArButton.hidden=false;delete (studio as HTMLElement).dataset.xrMode;applyTheme();onError(error instanceof Error?`Could not enter XR: ${error.message}`:'Could not enter XR.');if(xrSession)await xrSession.end().catch(()=>{});xrActiveMode=null;}finally{xrButton.disabled=xrArButton.disabled=false;}};
  xrButton.onclick=()=>{void enterXr('immersive-vr');};xrArButton.onclick=()=>{void enterXr('immersive-ar');};
  const overlaySelect=(event:Event)=>{if(event.target instanceof Element&&event.target.closest('button,input,select,.advanced-panel,.detail-sheet,.download-panel,.layers-panel,.section-panel'))event.preventDefault();};studio.addEventListener('beforexrselect',overlaySelect);
  if(navigator.xr){navigator.xr.isSessionSupported('immersive-vr').then(supported=>{if(supported&&!disposed){xrSupported=true;el.appendChild(xrButton);if(questVr){xrButton.disabled=true;xrButton.textContent='Preparing VR…';if(ready)void prepareVr();}}}).catch(()=>{});navigator.xr.isSessionSupported('immersive-ar').then(supported=>{if(supported&&!disposed){xrArSupported=true;el.appendChild(xrArButton);xrArButton.disabled=true;xrArButton.textContent='Preparing AR…';if(ready)void prepareVr();}}).catch(()=>{});}
  const updateXrPlacement=(xrFrame:XRFrame)=>{
   const tracked=xrBaseSpace&&xrFrame.getViewerPose(xrBaseSpace);if(!tracked)return;
   if(!xrAnchorMatrix)xrAnchorMatrix=new T.Matrix4().fromArray(tracked.transform.matrix);
   const projection=tracked.views[0]?.projectionMatrix;if(!projection)return;
   const xrArea=xrSession?.domOverlayState?viewArea():immersiveArea;
   xrPresentationCamera.aspect=Math.abs(projection[5]/projection[0]);xrPresentationCamera.fov=T.MathUtils.radToDeg(2*Math.atan(1/Math.abs(projection[5])));
   const values=cameraValues();applyFrame(xrPresentationCamera,new T.Vector3(),values,viewBounds,xrArea);
   const framed=decodeFrame(values,viewBounds,xrPresentationCamera.aspect,xrArea,xrPresentationCamera.fov),distance=framed.position.distanceTo(framed.target)*xrZoom;
   const direction=framed.position.sub(framed.target).normalize();xrPresentationCamera.position.copy(framed.target).addScaledVector(direction,distance);
   if(questVr&&xrActiveMode==='immersive-vr'&&(xrOrbit.yaw||xrOrbit.pitch)){
    xrOrbitQuaternion.setFromEuler(xrOrbitEuler.set(xrOrbit.pitch,xrOrbit.yaw,0,'YXZ'));
    xrPresentationCamera.position.sub(framed.target).applyQuaternion(xrOrbitQuaternion).add(framed.target);xrPresentationCamera.quaternion.premultiply(xrOrbitQuaternion);direction.applyQuaternion(xrOrbitQuaternion);
   }
   const right=new T.Vector3(1,0,0).applyQuaternion(xrPresentationCamera.quaternion),up=new T.Vector3(0,1,0).applyQuaternion(xrPresentationCamera.quaternion),tan=Math.tan(T.MathUtils.degToRad(xrPresentationCamera.fov/2)),[x,y]=screenAnchor(xrArea,values.slice(6,8));
   // Translate the presentation rig instead of changing the headset projection.
   xrPresentationCamera.position.addScaledVector(right,(.5-x)*2*distance*tan*xrPresentationCamera.aspect).addScaledVector(up,(y-.5)*2*distance*tan);
   xrPresentationCamera.updateMatrix();xrRig.matrix.multiplyMatrices(xrPresentationCamera.matrix,xrAnchorMatrix.clone().invert());xrRig.matrix.decompose(xrRig.position,xrRig.quaternion,xrRig.scale);xrRig.updateMatrixWorld(true);
   xrRenderCamera.near=Math.max(.00001,Math.min(.05,distance/500));xrRenderCamera.far=Math.max(100,distance*10);
   const hudDistance=distance*.45;xrInfo.position.copy(xrPresentationCamera.position).addScaledVector(direction,-hudDistance).addScaledVector(up,-.84*hudDistance*tan);xrInfo.scale.set(1.3*hudDistance*tan*xrPresentationCamera.aspect,.28*hudDistance*tan,1);
  };
  const sectionWheel=(event:WheelEvent)=>{const section=latest.current.section;if(remoteRef.current||!sectionToolRef.current||!section?.enabled||event.ctrlKey||event.metaKey)return;event.preventDefault();event.stopImmediatePropagation();sectionPositionRef.current?.(Math.max(0,Math.min(1,section.position+Math.sign(event.deltaY)*Math.min(.04,Math.abs(event.deltaY)/2500))));};
  renderer.domElement.addEventListener('wheel',sectionWheel,{capture:true,passive:false});
  const drawLabels=(s:SceneState)=>{
   labelLayer.replaceChildren();if(s.labels===false)return;const w=el.clientWidth,h=el.clientHeight;labelLayer.setAttribute('viewBox',`0 0 ${w} ${h}`);
   if(amount>.04&&explosionLayout&&!s.isolate){
    const progress=amount*explosionLayout.steps,from=Math.min(explosionLayout.steps,Math.floor(progress)),to=Math.min(explosionLayout.steps,from+1),mix=smoothStep(progress-from);
    const depth=Math.min(4,Math.max(1,Math.floor(progress+.05))),first=new Map(explosionLayout.stages[from].groups.map(group=>[group.id,group]));
    const candidates=explosionLayout.stages[to].groups.filter(group=>group.depth<=depth),names=new Set(candidates.map(group=>group.id)),parents=new Set<string>();
    for(const group of candidates){let parent=group.id.lastIndexOf('/');while(parent>0){const id=group.id.slice(0,parent);if(names.has(id))parents.add(id);parent=id.lastIndexOf('/');}}
    const area=currentLabelArea,occupied:{left:number;right:number;top:number;bottom:number}[]=[];
    for(const group of candidates.filter(group=>!parents.has(group.id)).sort((a,b)=>b.count-a.count).slice(0,w<768?50:180)){
     const previous=first.get(group.id)??group;
     const x=T.MathUtils.lerp(previous.x,group.x,mix),y=T.MathUtils.lerp(previous.y,group.y,mix),z=T.MathUtils.lerp(previous.z,group.z,mix),height=T.MathUtils.lerp(previous.height,group.height,mix);
     projected.set(x,y+height/2,z).project(camera);if(projected.z< -1||projected.z>1)continue;
     const screenX=(projected.x+1)*w/2,screenY=(1-projected.y)*h/2-10;
     const display=displayLaterality(group.name),name=display.label.length>32?`${display.label.slice(0,31)}…`:display.label,labelWidth=Math.min(220,name.length*6.2+14);
     const box={left:screenX-labelWidth/2,right:screenX+labelWidth/2,top:screenY-17,bottom:screenY+2};
     if(box.left<area.left||box.right>area.right||box.top<area.top||box.bottom>area.bottom||occupied.some(other=>box.left<other.right+4&&box.right>other.left-4&&box.top<other.bottom+3&&box.bottom>other.top-3))continue;
     occupied.push(box);
     const element=document.createElementNS('http://www.w3.org/2000/svg','g');element.setAttribute('class',`explosion-group-label ${lateralityClass(display.side)}`);
     const background=document.createElementNS('http://www.w3.org/2000/svg','rect');background.setAttribute('x',String(box.left));background.setAttribute('y',String(box.top));background.setAttribute('width',String(labelWidth));background.setAttribute('height','19');background.setAttribute('rx','5');element.appendChild(background);
     const label=document.createElementNS('http://www.w3.org/2000/svg','text');label.setAttribute('x',String(screenX));label.setAttribute('y',String(screenY-3));label.setAttribute('text-anchor','middle');label.textContent=name;element.appendChild(label);const title=document.createElementNS('http://www.w3.org/2000/svg','title');title.textContent=group.name;element.appendChild(title);labelLayer.appendChild(element);
    }
   }
   if(s.selectionGroup)return;
   const anchors:LabelAnchor[]=[],named=new Set<string>();
   for(const id of s.selected){const i=partIndices.get(id)??-1,mesh=pickers[i];if(i<0||!mesh||data[i*4+3]<.5)continue;
    const name=structureName(atlas.parts[i].name);if(named.has(name))continue;
    let anchor=labelAnchors.get(i);if(!anchor){const g=mesh.geometry,positions=g.getAttribute('position'),indices=g.index!;let best=Infinity;const center=centers[i];for(let j=0;j<indices.count;j+=Math.max(1,Math.floor(indices.count/1800))*3){const point=new T.Vector3();for(let k=0;k<3;k++)point.add(new T.Vector3().fromBufferAttribute(positions,indices.getX(Math.min(j+k,indices.count-1))));point.multiplyScalar(1/3);if(!keptBySections(point.clone().add(mesh.position)))continue;const d=point.distanceToSquared(center);if(d<best){best=d;anchor=point;}}if(!anchor)continue;labelAnchors.set(i,anchor);}
    const world=anchor.clone().applyMatrix4(mesh.matrixWorld);if(!keptBySections(world))continue;const point=world.clone().project(camera);
    if(point.z< -1||point.z>1||Math.abs(point.x)>1||Math.abs(point.y)>1)continue;
    const display=displayLaterality(name);named.add(name);anchors.push({id,x:(point.x+1)*w/2,y:(1-point.y)*h/2,text:display.label,fullText:name,laterality:display.side});
   }
   const labelBounds=currentLabelArea;
   const measure=(text:string)=>{const context=labelMeasure.getContext('2d');if(!context)return text.length*6.8;context.font='12px Arial';return context.measureText(text).width;};
   const placements=layoutAnatomyLabels(anchors,labelBounds,measure);
   const leaders=document.createElementNS('http://www.w3.org/2000/svg','g');leaders.setAttribute('class','selected-anatomy-label');labelLayer.appendChild(leaders);
   const pills=document.createElementNS('http://www.w3.org/2000/svg','g');pills.setAttribute('class','selected-anatomy-label');labelLayer.appendChild(pills);
   const make=(parent:SVGGElement,name:string,attrs:Record<string,string>)=>{const element=document.createElementNS('http://www.w3.org/2000/svg',name);Object.entries(attrs).forEach(([key,value])=>element.setAttribute(key,value));parent.appendChild(element);return element;};
   for(const placed of placements){
    make(leaders,'line',{class:lateralityClass(placed.laterality??null),x1:String(placed.x),y1:String(placed.y),x2:String(placed.leaderX),y2:String(placed.leaderY)});
    make(leaders,'circle',{class:lateralityClass(placed.laterality??null),cx:String(placed.x),cy:String(placed.y),r:'3'});
   }
   for(const placed of placements){
    make(pills,'rect',{class:lateralityClass(placed.laterality??null),x:String(placed.left),y:String(placed.top),width:String(placed.right-placed.left),height:String(placed.bottom-placed.top),rx:'5'});
    const label=make(pills,'text',{class:lateralityClass(placed.laterality??null),x:String(placed.left+8),y:String(placed.leaderY+3.5)});label.textContent=placed.label;
    const title=make(pills,'title',{});title.textContent=placed.fullText??placed.text;
   }
  };
  const clock=new T.Clock();let lastExtent=-1,lastCamera:number[]|undefined,lastFocus=0,lastRotationKey='',lastSectionOrientation='',lastSectionPoint:T.Vector3|undefined,selectionStencil=false,lastRemotePose:ScenePose|null=null,panelsAnimating=false;
  const animate=(now=performance.now())=>{
   if(disposed)return;const presenting=renderer.xr.isPresenting;if(!presenting)frame=requestAnimationFrame(animate);const dt=Math.min(clock.getDelta(),.05),requested=latest.current;if(!presenting)attachPreview();
   const following=!!remoteRef.current,reducedMotion=window.matchMedia('(prefers-reduced-motion: reduce)').matches;
   let depthMoving=false;
   if(requested.depth===undefined){depthAmount=undefined;depthVelocity=0;}
   else{
    depthAmount??=ready?depthPosition(lastState??requested):requested.depth;
    depthMoving=Math.abs(depthAmount-requested.depth)>.00001||Math.abs(depthVelocity)>.0001;
    if(depthMoving&&!reducedMotion){const motion=dampMotion(depthAmount,depthVelocity,requested.depth,Math.max(depthAmount,requested.depth)*DEPTH_LAYERS.length<=2?.32:.22,dt);depthAmount=T.MathUtils.clamp(motion.value,0,1);depthVelocity=motion.velocity;}
    else{depthAmount=requested.depth;depthVelocity=0;depthMoving=false;}
   }
   const s=depthMoving&&depthAmount!==undefined?setDepthPosition(requested,depthAmount,atlas.parts):requested;
   const transition=transitionRef.current,savedTransition=!!transition?.id&&transition.id!==lastTransitionId&&!!s.camera&&validCamera(s.camera);
   const selectionCleared=!!lastState?.selected.length&&!s.selected.length;
   const newCamera=!!s.camera&&s.camera!==lastCamera&&!selectionCleared;
   const actionTransition=s.view!==lastView||s.reset!==lastReset||(s.region??'all')!==lastRegion||!!s.focus&&s.focus!==lastFocus||s.isolate&&s.selected.join(',')+':'+s.reset!==lastIsolate||!s.isolate&&!!lastIsolate&&!selectionCleared||newCamera;
   const startPose=ready&&!following&&!reducedMotion&&lastView!==''&&(savedTransition||actionTransition)?captureCameraPose(camera,controls.target):null;
   if(savedTransition)lastTransitionId=transition!.id;
   if(reducedMotion||!transition&&cameraTween?.source==='saved')cancelCameraTween();
   controls.enabled=!following&&!presenting;controls.enableDamping=!following&&!presenting&&!cameraTween;
   if(!presenting){const animating=studio.ownerDocument.getAnimations().some(animation=>{const target=(animation.effect as KeyframeEffect|null)?.target;return target instanceof Element&&(target.closest(panelSelectors)||target.querySelector(panelSelectors));});if(animating||panelsAnimating)layoutDirty=true;panelsAnimating=animating;}
   if(layoutDirty){
    layoutDirty=false;const areas=readLayoutAreas();currentArea=areas.view;currentLabelArea=areas.labels;
    const aspect=Math.max(1,el.clientWidth)/Math.max(1,el.clientHeight);if(camera.aspect!==aspect){const view=camera.view;camera.aspect=aspect;if(view?.enabled)camera.setViewOffset(aspect,1,view.offsetX/view.fullWidth*aspect,view.offsetY/view.fullHeight,aspect,1);}
    camera.updateProjectionMatrix();if(!presenting){syncCaptureResolution();renderer.setSize(el.clientWidth,el.clientHeight);}
    // Follow the initial detail-panel entrance before handing camera placement
    // to the user. Later panel changes preserve their chosen camera.
    if(initialSelectionFitPending&&!s.camera)lastFocus=-1;
    if(!areaInitialized)fit(s.view,Math.min(1,amount*6));
    areaInitialized=true;dirty=true;
    renderer.domElement.dataset.viewArea=JSON.stringify(currentArea);renderer.domElement.dataset.viewAnchor=JSON.stringify(screenAnchor(currentArea,preferredAnchor(currentArea)));
   }
   if(ready&&!panelsAnimating)initialSelectionFitPending=false;
   if(presenting){const xrFrame=renderer.xr.getFrame(),sources=xrSession?.inputSources;const left=[...(sources??[])].find(source=>source.handedness==='left');const stick=left?.gamepad?.axes[3]??0;
    if(!following&&Math.abs(stick)>.15)xrZoom=T.MathUtils.clamp(xrZoom*Math.exp(stick*dt*1.25),.25,4);
    if(questVr&&xrActiveMode==='immersive-vr'&&xrOrbitHeld&&!following&&left&&xrFrame){const space=renderer.xr.getReferenceSpace(),pose=space&&xrFrame.getPose(left.targetRaySpace,space);if(pose){const {x,y,z,w}=pose.transform.orientation,aim=xrAim.set(0,0,-1).applyQuaternion(xrControllerQuaternion.set(x,y,z,w));const current={yaw:Math.atan2(aim.x,-aim.z),pitch:Math.asin(T.MathUtils.clamp(aim.y,-1,1))};xrOrbit=advanceControllerOrbit(xrLastLeftAim,current,xrOrbit);xrLastLeftAim=current;}}
    else xrLastLeftAim=null;
    const selected=s.selected[0],index=selected===undefined?undefined:partIndices.get(selected);if(index!==undefined)updateXrInfo(structureName(atlas.parts[index].name));}
   const nextSectionActive=enabledSections(s).length?1:0;if(sectionActiveUniform.value!==nextSectionActive){sectionActiveUniform.value=nextSectionActive;dirty=true;}
   const nextContext=nextSectionActive?1:s.selected.length?(s.contextOpacity??1):1;
   if(contextUniform.value!==nextContext){contextUniform.value=nextContext;dirty=true;}
   const nextStencil=s.selected.length>0&&nextContext<.999;
   if(nextStencil!==selectionStencil){selectionStencil=nextStencil;
    for(const material of mats.values()){material.stencilWrite=nextStencil;material.stencilRef=1;material.stencilFunc=T.AlwaysStencilFunc;material.stencilZPass=T.ReplaceStencilOp;}
    for(const material of ghostMats.values()){material.stencilWrite=nextStencil;material.stencilRef=1;material.stencilFunc=T.NotEqualStencilFunc;material.stencilZPass=T.KeepStencilOp;}
    dirty=true;
   }
   const changed=lastState!==s;

   if(changed){previewDirty=true;labelAnchors.clear();const box=sectionBox(s),section=s.section??{enabled:false,axis:'axial',position:.38,flip:true};
    clipPlanes.forEach(plane=>{plane.constant=10000;});
    for(const cut of enabledSections(s)){const equation=sectionPlaneEquation(sectionBounds,cut.section),plane=clipPlanes[cut.index];plane.normal.fromArray(equation.normal);plane.constant=equation.constant;}
    const point=sectionPoint(sectionBounds,section),normal=sectionNormal(section),quaternion=new T.Quaternion().setFromUnitVectors(new T.Vector3(0,0,1),new T.Vector3().fromArray(normal)),size=box.getSize(new T.Vector3());
    const span=Math.max(section.axis==='axial'?size.x:size.y,section.axis==='sagittal'?size.z:size.x,size.z,.2)*1.12;
    guideFill.visible=!lightweightActive&&sectionToolRef.current&&!enabledSections(s).length;guideBorder.visible=guideFill.visible;guideFill.position.fromArray(point);guideBorder.position.fromArray(point);guideFill.quaternion.copy(quaternion);guideBorder.quaternion.copy(quaternion);guideFill.scale.set(span,span,1);guideBorder.scale.set(span,span,1);
    dirty=true;}
   const explosionRequested=s.explode!==lastExplode;
   if(explosionRequested&&!savedTransition&&!newCamera)cancelCameraTween();
   if(explosionRequested)implosionFillActive=s.explode<lastExplode&&!s.isolate;
   lastExplode=s.explode;
   const moving=Math.abs(amount-s.explode)>.00001||Math.abs(explosionVelocity)>.0001;
   if(moving&&!reducedMotion){const motion=dampMotion(amount,explosionVelocity,s.explode,Math.max(amount,s.explode)*(explosionLayout?.steps??1)<=2?.32:.22,dt);amount=T.MathUtils.clamp(motion.value,0,1);explosionVelocity=motion.velocity;dirty=true;}
   else if(amount!==s.explode||explosionVelocity){amount=s.explode;explosionVelocity=0;lastExtent=-1;dirty=true;}
   if(changed||moving||lastExtent<0){
    const selection=new Set(s.selected),activeCut=enabledSections(s).length>0;
    const visibleParts=atlas.parts.filter(p=>activeCut?sectionPartVisible(p,s):partVisible(p,s));
    const nextLayoutKey=visibleParts.map(p=>p.id).join(',')+':'+hierarchyRef.current+':'+(hierarchyRef.current==='guest'?guestKeyRef.current:'');
    const explosionChanged=!!lastState&&s.explode!==lastState.explode,layoutChanged=nextLayoutKey!==layoutKey||explosionChanged&&camera.aspect!==explosionAspect;
    if(layoutChanged){explosionLayout=createHierarchicalExplosionLayout(atlas,anatomyEntries,new Set(visibleParts.map(part=>part.id)),hierarchyRef.current,camera.aspect,guestNodesRef.current);layoutKey=nextLayoutKey;explosionAspect=camera.aspect;stepsCallback.current?.(explosionLayout.steps);}
    const stages=explosionLayout!.stages,progress=amount*explosionLayout!.steps,from=Math.min(stages.length-1,Math.floor(progress)),to=Math.min(stages.length-1,from+1),blend=smoothStep(progress-from);
    packingWidth=T.MathUtils.lerp(stages[from].width,stages[to].width,blend);packingHeight=T.MathUtils.lerp(stages[from].height,stages[to].height,blend);
    if(explosionRequested||layoutChanged&&amount>.001&&!selectionCleared){
     explosionFitPending=s.explode>0&&!s.isolate&&!implosionFillActive&&!s.camera;
     if(explosionFitPending){const box=nextExplosionBounds(explosionLayout!,atlas.parts,new Set(visibleParts.map(part=>part.id)),amount,s.explode);explosionFitBounds.min.fromArray(box.min);explosionFitBounds.max.fromArray(box.max);}
    }
    const clusterStageIndex=Math.min(stages.length-1,Math.floor(progress+.05)),clusterStage=stages[clusterStageIndex],selectedClusters=new Set<number>();
    if(clusterStageIndex>0)for(const id of s.selected){const i=partIndices.get(id),cluster=i===undefined?-1:clusterStage.clusterIds[i];if(cluster>=0)selectedClusters.add(cluster);}
    const rotationKey=s.selected.join(',')+':'+s.reset+':'+clusterStageIndex+':'+[...selectedClusters].join(',');
    if(layoutChanged||rotationKey!==lastRotationKey||s.explode<=.04)selectedRotation.identity();lastRotationKey=rotationKey;
    rotationPartIndices.clear();rotationData.fill(0);
    for(const cluster of selectedClusters)for(const i of clusterStage.clusters[cluster]){rotationPartIndices.add(i);rotationData[i*4]=255;}

    viewBounds.makeEmpty();explosionVisibleBounds.makeEmpty();
    atlas.parts.forEach((p,i)=>{
     const c=centers[i],j=i*3,dx=T.MathUtils.lerp(stages[from].positions[j],stages[to].positions[j],blend)-c.x,dy=T.MathUtils.lerp(stages[from].positions[j+1],stages[to].positions[j+1],blend)-c.y,dz=T.MathUtils.lerp(stages[from].positions[j+2],stages[to].positions[j+2],blend)-c.z;
     const selected=selection.has(p.id),checked=sectionPartVisible(p,s);
     data.set([dx,dy,dz,(activeCut?checked:partVisible(p,s))&&(!isSkinPart(p)||surfaceReady||selected||activeCut)?1:0],i*4);selectedData[i*4]=selected?255:0;selectedData[i*4+1]=Math.round(partLayerOpacity(p,s)*255);previewData[i*4]=(hasSurfacePreview?isSkinPart(p):checked)?255:0;
     if(activeCut?checked:partVisible(p,s)){explosionVisibleBounds.expandByPoint(framingPoint.set(bounds[i].min.x+dx,bounds[i].min.y+dy,bounds[i].min.z+dz));explosionVisibleBounds.expandByPoint(framingPoint.set(bounds[i].max.x+dx,bounds[i].max.y+dy,bounds[i].max.z+dz));}
     if(selection.size?selected:activeCut?checked:partVisible(p,s)){viewBounds.expandByPoint(framingPoint.set(bounds[i].min.x+dx,bounds[i].min.y+dy,bounds[i].min.z+dz));viewBounds.expandByPoint(framingPoint.set(bounds[i].max.x+dx,bounds[i].max.y+dy,bounds[i].max.z+dz));}
     const mesh=pickers[i];if(mesh){mesh.quaternion.identity();mesh.position.set(dx,dy,dz);mesh.updateMatrix();mesh.updateMatrixWorld(true);}
    });if(viewBounds.isEmpty())viewBounds.copy(modelBounds);syncSelectedRotation();partTexture.needsUpdate=true;selectionTexture.needsUpdate=true;previewTexture.needsUpdate=true;rotationTexture.needsUpdate=true;lastState=s;lastExtent=amount;dirty=true;if(changed)queueCaps(s);
   }
   if(s.view!==lastView||s.reset!==lastReset||s.region!==lastRegion){if(s.focus)lastFocus=-1;fit(amount>.04?'front':s.view,Math.min(1,amount*6));lastView=s.view;lastReset=s.reset;lastRegion=s.region??'all';}
   const isolateKey=s.isolate?s.selected.join(',')+':'+s.reset:'';
   const focusChanged=(s.focus??0)!==lastFocus;lastFocus=s.focus??0;
   if(isolateKey!==lastIsolate||(s.isolate&&moving)||focusChanged&&!!s.focus){
    // URL selection is available before its mesh finishes loading. Frame from
    // catalogue bounds immediately instead of consuming focus on an empty box.
    if(s.isolate||focusChanged&&!!s.focus){const box=new T.Box3();atlas.parts.forEach((p,i)=>{if(s.selected.includes(p.id))box.union(worldBox.copy(bounds[i]).applyMatrix4(pickers[i]?.matrixWorld??new T.Matrix4()));});
     if(!box.isEmpty()){const center=box.getCenter(new T.Vector3()),area=viewArea(),direction=directionFor(s.view).normalize();setFrameOffset(camera,area);const distance=Math.max(.006,frameDistance(box,direction,camera.up,camera.aspect,area,camera.fov)*1.35);controls.maxDistance=Math.max(.12,distance*1.6);controls.target.copy(center);camera.position.copy(center).addScaledVector(direction,distance);controls.update();dirty=true;}
    }else if(lastIsolate&&!selectionCleared){camera.clearViewOffset();fit(s.view,Math.min(1,amount*6));}
    lastIsolate=isolateKey;
   }
   const selectionKey=s.selected.join(',')+':'+s.reset;if(selectionKey!==lastSelection||moving&&s.selected.length>0){const pivot=selectionPivot();if(pivot){controls.target.copy(pivot);controls.update();}lastSelection=selectionKey;dirty=true;}
   const canOrbit=amount<=.04||s.isolate;controls.enableRotate=canOrbit;controls.mouseButtons.LEFT=canOrbit?T.MOUSE.ROTATE:T.MOUSE.PAN;controls.touches.ONE=canOrbit?T.TOUCH.ROTATE:T.TOUCH.PAN;controls.autoRotate=!following&&!presenting&&s.rotate&&!s.isolate&&canOrbit&&!cameraTween;controls.autoRotateSpeed=.65;if(!presenting)controls.update();if(controls.autoRotate)dirty=true;
   if(!s.camera)lastCamera=undefined;
   if(s.camera&&s.camera!==lastCamera){applyFrame(camera,controls.target,s.camera,viewBounds,viewArea());controls.maxDistance=Math.max(controls.maxDistance,camera.position.distanceTo(controls.target)*1.05);controls.update();lastCamera=s.camera;dirty=true;}
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
     const area=viewArea(),tan=Math.tan(T.MathUtils.degToRad(camera.fov/2));
     const distance=Math.max(.12,spanX/Math.max(.08,area.right-area.left)/camera.aspect,spanY/Math.max(.08,area.bottom-area.top))/(2*tan)*(referenceModel ? .9 : 1.12);
     setFrameOffset(camera,area);
     controls.maxDistance=Math.max(controls.maxDistance,distance*1.5);
     controls.target.copy(point);camera.position.copy(point).addScaledVector(direction,distance);controls.update();dirty=true;
    }
    else if(lastSectionPoint&&orientation===lastSectionOrientation&&!point.equals(lastSectionPoint)){const delta=point.clone().sub(lastSectionPoint);controls.target.add(delta);camera.position.add(delta);controls.update();dirty=true;}
    lastSectionPoint=point;
   }else if(lastSectionOrientation){camera.up.set(0,1,0);if(!s.camera)fit(s.view,Math.min(1,amount*6));lastSectionPoint=undefined;dirty=true;}
   lastSectionOrientation=orientation;
   if(startPose){
    let destination=captureCameraPose(camera,controls.target);
    if(s.camera&&validCamera(s.camera)&&explosionLayout){
     const finalParts=atlas.parts.filter(part=>enabledSections(requested).length?sectionPartVisible(part,requested):partVisible(part,requested)),ids=new Set(requested.selected.length?requested.selected:finalParts.map(part=>part.id));
     const finalLayout=s===requested?explosionLayout:createHierarchicalExplosionLayout(atlas,anatomyEntries,new Set(finalParts.map(part=>part.id)),hierarchyRef.current,camera.aspect,guestNodesRef.current);
     const box=explosionBoundsAt(finalLayout,atlas.parts,ids,requested.explode),finalBounds=new T.Box3(new T.Vector3().fromArray(box.min),new T.Vector3().fromArray(box.max));
     destination=cameraPoseForFrame(s.camera,finalBounds.isEmpty()?modelBounds:finalBounds,camera,viewArea());
    }
    cameraTween={from:startPose,to:destination,start:now,duration:1100,source:savedTransition?'saved':'action'};controls.enableDamping=false;
   }
   if(cameraTween){
    const {from,to,start,duration}=cameraTween,progress=Math.min(1,(now-start)/duration);
    interpolateCameraPose(camera,controls.target,from,to,smoothStep(progress));
    controls.maxDistance=Math.max(controls.maxDistance,camera.position.distanceTo(controls.target)*1.05);controls.update();dirty=true;
    if(progress===1)cancelCameraTween();
   }
   const remote=remoteRef.current?.current;
   if(following&&remote&&remote.model===remoteModelRef.current&&remote!==lastRemotePose){
    applyFrame(camera,controls.target,remote.camera,viewBounds,viewArea());if(remote.camera.length!==12)camera.up.fromArray(remote.up);controls.maxDistance=Math.max(controls.maxDistance,camera.position.distanceTo(controls.target)*1.05);controls.update();
    selectedRotation.fromArray(remote.rotation).normalize();syncSelectedRotation();lastRemotePose=remote;dirty=true;
   }
   const zoomStart=camera.position.clone(),zoomStartDistance=camera.position.distanceTo(controls.target);
   let fittedZoom=false;
   if(explosionFitPending&&!cameraTween){
    for(const i of rotationPartIndices){const mesh=pickers[i];if(mesh&&data[i*4+3]>.5)explosionFitBounds.union(worldBox.copy(bounds[i]).applyMatrix4(mesh.matrixWorld));}
    fittedZoom=zoomOutToFit(camera,explosionFitBounds,viewArea());
   }
   if(implosionFillActive&&!cameraTween&&!s.camera){
    const box=explosionVisibleBounds.clone();
    for(const i of rotationPartIndices){const mesh=pickers[i];if(mesh&&data[i*4+3]>.5)box.union(worldBox.copy(bounds[i]).applyMatrix4(mesh.matrixWorld));}
    let fitted=zoomOutToFit(camera,box,viewArea());
    fitted=zoomInToFill(camera,box,viewArea(),Math.max(0,camera.position.distanceTo(controls.target)-controls.minDistance))||fitted;
    fittedZoom=fitted||fittedZoom;
   }
   if((explosionFitPending||implosionFillActive)&&!cameraTween&&!s.camera){
    const desiredDistance=camera.position.distanceTo(controls.target),motion=dampMotion(Math.log(Math.max(.000001,zoomStartDistance)),zoomVelocity,Math.log(Math.max(.000001,desiredDistance)),.16,dt);
    const settled=Math.abs(motion.value-Math.log(Math.max(.000001,desiredDistance)))<.00001&&Math.abs(motion.velocity)<.0001;
    zoomVelocity=reducedMotion||settled?0:motion.velocity;
    if(!reducedMotion&&!settled){camera.position.copy(controls.target).addScaledVector(zoomStart.sub(controls.target).normalize(),Math.exp(motion.value));camera.updateMatrixWorld(true);dirty=true;}
    // Fit the moving geometry as well as the upcoming mark; gradual zoom must
    // never allow the current parts to cross a panel or viewport edge.
    fittedZoom=zoomOutToFit(camera,explosionVisibleBounds,viewArea())||fittedZoom;
    if(fittedZoom||zoomVelocity){controls.maxDistance=Math.max(controls.maxDistance,camera.position.distanceTo(controls.target)*1.25);controls.update();dirty=true;}
    if(!moving&&settled){explosionFitPending=false;implosionFillActive=false;}
   }else if(s.camera){explosionFitPending=false;implosionFillActive=false;zoomVelocity=0;}
   // Keep depth precision near the anatomy as the camera moves. A fixed tiny near
   // plane causes close skin and muscle surfaces to fight at whole-body distances.
   const near=Math.max(.00001,Math.min(.05,camera.position.distanceTo(controls.target)/500));
   if(!presenting&&Math.abs(camera.near-near)>near*.01){camera.near=near;camera.updateProjectionMatrix();dirty=true;}
   if(lightweightActive&&vrAnatomy&&dirty){const solidQuest=questVr&&xrActiveMode==='immersive-vr';vrAnatomy.sync(s,data,pickers,solidQuest?1:contextUniform.value,solidQuest);renderer.domElement.dataset.xrVisibleTriangles=String(vrAnatomy.mesh.geometry.drawRange.count/3);}
   if(presenting){const xrFrame=renderer.xr.getFrame();updateXrPlacement(xrFrame);xrRays.forEach((ray,i)=>{const source=xrSession?.inputSources[i];ray.visible=!!source&&xrPose(xrFrame,source);if(ray.visible){ray.position.copy(xrOrigin);ray.quaternion.copy(xrQuaternion);}});}
   if(dirty||presenting){
    renderer.render(scene,presenting?xrRenderCamera:camera);if(lightweightActive&&dirty){renderer.domElement.dataset.xrRenderCalls=String(renderer.info.render.calls);renderer.domElement.dataset.xrRenderedTriangles=String(renderer.info.render.triangles);}
    let values:number[],viewUp:number[];
    if(presenting){const tracked=renderer.xr.getCamera();tracked.matrixWorld.decompose(xrHeadCamera.position,xrHeadCamera.quaternion,xrHeadCamera.scale);xrHeadCamera.up.set(0,1,0).applyQuaternion(xrHeadCamera.quaternion);xrHeadCamera.fov=xrPresentationCamera.fov;xrHeadCamera.aspect=xrPresentationCamera.aspect;const target=xrHeadCamera.position.clone().addScaledVector(new T.Vector3(0,0,-1).applyQuaternion(xrHeadCamera.quaternion),xrPresentationCamera.position.distanceTo(controls.target));values=captureFrame(xrHeadCamera,target,viewBounds,xrSession?.domOverlayState?viewArea():immersiveArea);viewUp=xrHeadCamera.up.toArray();}
    else{values=cameraValues();viewUp=camera.up.toArray();drawLabels(s);}
    cameraCallback.current?.(values);if(!following)poseCallback.current?.({camera:values,up:viewUp,rotation:selectedRotation.toArray()});renderer.domElement.dataset.cameraFrame=JSON.stringify(values);targets=[];
    if(!presenting&&amount>.45){
     const hasSolid=solidContext();
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
    if(!pageFrameListeners.size)rememberPageFrame();for(const listener of pageFrameListeners)listener({canvas:renderer.domElement,labels:labelLayer,transitioning:!!cameraTween||moving||depthMoving||explosionFitPending||implosionFillActive},now);dirty=false;
   }
   if(!presenting&&previewCanvas&&(previewDirty||changed))drawPreview(s);
   if(!presenting&&ready){const key=[s.selected.join(','),s.visible.join(','),s.region,s.skinOpacity,s.depth,s.explode,...selectedRotation.toArray().map(v=>v.toFixed(2)),...camera.position.toArray().map(v=>v.toFixed(2)),...controls.target.toArray().map(v=>v.toFixed(2))].join(':');if((key!==lastCoveringKey||s.hidden!==lastCoveringHidden||s.depthHidden!==lastCoveringDepth)&&performance.now()-lastCoveringScan>250){scanCovering(s);lastCoveringKey=key;lastCoveringHidden=s.hidden;lastCoveringDepth=s.depthHidden;lastCoveringScan=performance.now();}}

  };animate();
  const capturePageFrame:SceneFrameCapture=Object.assign((size?:CaptureSize)=>{captureSize=size;const resolutionChanged=syncCaptureResolution();if(!pageCaptureCanvas||resolutionChanged){pageCaptureCanvas??=document.createElement('canvas');renderer.render(scene,renderer.xr.isPresenting?xrRenderCamera:camera);rememberPageFrame();}return {canvas:pageCaptureCanvas,labels:labelLayer,transitioning:!!cameraTween||Math.abs(amount-latest.current.explode)>.00001||depthAmount!==undefined&&Math.abs(depthAmount-(latest.current.depth??depthAmount))>.00001||explosionFitPending||implosionFillActive};},{release:()=>{pageFrameListeners.clear();pageCaptureCanvas=null;captureSize=undefined;syncCaptureResolution();},invalidate:()=>{dirty=true;},subscribe:(listener:Parameters<SceneFrameCapture['subscribe']>[0],size:CaptureSize)=>{captureSize=size;syncCaptureResolution();pageCaptureCanvas=null;pageFrameListeners.add(listener);dirty=true;return()=>{pageFrameListeners.delete(listener);};}});
  onFrameCapture?.(capturePageFrame);
  onCapture?.(async(options={labels:false,background:true})=>{
   if(!ready||disposed)throw new Error('The 3D view is still loading.');
   const previousAlpha=renderer.getClearAlpha();renderer.setClearAlpha(options.background?previousAlpha:0);
   const output=document.createElement('canvas');output.width=renderer.domElement.width;output.height=renderer.domElement.height;
   const context=output.getContext('2d');
   try{if(!context)throw new Error('The PNG could not be created.');renderer.render(scene,renderer.xr.isPresenting?xrRenderCamera:camera);context.drawImage(renderer.domElement,0,0);}
   finally{renderer.setClearAlpha(previousAlpha);dirty=true;}
   if(options.labels&&labelLayer.childElementCount){
    const clone=labelLayer.cloneNode(true) as SVGSVGElement,source=[labelLayer,...labelLayer.querySelectorAll('*')],target=[clone,...clone.querySelectorAll('*')];
    source.forEach((element,index)=>{const style=getComputedStyle(element);for(const property of ['fill','stroke','stroke-width','font-family','font-size','font-weight','opacity','paint-order','stroke-linejoin'])target[index].setAttribute(property,style.getPropertyValue(property));});
    clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('width',String(el.clientWidth));clone.setAttribute('height',String(el.clientHeight));
    const url=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)],{type:'image/svg+xml'}));
    try{const image=new Image();image.src=url;await image.decode();context.drawImage(image,0,0,output.width,output.height);}finally{URL.revokeObjectURL(url);}
   }
   return new Promise<Blob>((resolve,reject)=>output.toBlob(blob=>blob?resolve(blob):reject(new Error('The PNG could not be created.')),'image/png'));
  });
  const contextLost=(e:Event)=>{e.preventDefault();onError('The 3D session was paused by your device. Reload to continue.');};renderer.domElement.addEventListener('webglcontextlost',contextLost);
  return()=>{onCapture?.(null);onFrameCapture?.(null);disposed=true;abort.abort();vrAnatomy?.mesh.removeFromParent();vrAnatomy?.dispose();cancelAnimationFrame(frame);renderer.setAnimationLoop(null);renderer.xr.removeEventListener('sessionend',xrEnded);if(xrSession){xrSession.removeEventListener('select',xrSelect);xrSession.removeEventListener('selectstart',xrSelectStart);xrSession.removeEventListener('selectend',xrSelectEnd);void xrSession.end();}xrButton.remove();xrArButton.remove();delete (studio as HTMLElement).dataset.xrMode;xrRays.forEach(ray=>{ray.geometry.dispose();(ray.material as T.Material).dispose();});xrInfoTexture.dispose();(xrInfo.material as T.Material).dispose();clearTimeout(capTimer);clearCaps();stencilCaps.dispose();observer.disconnect();layoutObserver.disconnect();studio.removeEventListener('beforexrselect',overlaySelect);theme.removeEventListener('change',applyTheme);renderer.domElement.removeEventListener('wheel',cancelCameraTween);renderer.domElement.removeEventListener('wheel',sectionWheel,true);renderer.domElement.removeEventListener('contextmenu',preventContextMenu);previewCanvas?.removeEventListener('pointerdown',previewDown);previewCanvas?.removeEventListener('pointermove',previewMove);previewCanvas?.removeEventListener('pointerup',previewUp);previewCanvas?.removeEventListener('pointercancel',previewUp);controls.dispose();geometries.forEach(g=>g.dispose());materials.forEach(m=>m.dispose());textures.forEach(t=>t.dispose());textureUrls.forEach(url=>URL.revokeObjectURL(url));guideGeometry.dispose();guideBorder.geometry.dispose();guideFill.material.dispose();previewPlanes.forEach(({geometry,fill,edge})=>{geometry.dispose();edge.geometry.dispose();fill.material.dispose();edge.material.dispose();});previewTarget.dispose();scene.traverse(o=>{if(o instanceof T.Mesh&&!geometries.includes(o.geometry)&&o!==guideFill){o.geometry.dispose();const ms=Array.isArray(o.material)?o.material:[o.material];ms.forEach(m=>m.dispose());}});env.dispose();partTexture.dispose();selectionTexture.dispose();previewTexture.dispose();rotationTexture.dispose();hover.remove();labelLayer.remove();renderer.dispose();renderer.domElement.remove();};
 },[atlas]);
 return <div className="scene" ref={host}/>;
}
