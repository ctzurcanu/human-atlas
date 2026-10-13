import type {PerspectiveCamera,Vector3} from 'three';
import {captureExplosionView,restoreExplosionView,type ExplosionView} from './explosion-view';

export type SceneCameraIntent={model:string;navigation:number;view:string;reset:number;focus:number;camera?:number[];transition:number};
export type RawSceneCamera=ExplosionView&{aspect:number;near:number;far:number;filmGauge:number;filmOffset:number};
export type SceneCameraRestart={intent:SceneCameraIntent;pose:RawSceneCamera;minDistance:number;maxDistance:number;explosion:{amount:number;velocity:number;requested:number;start:ExplosionView|null;fitting:boolean}};

// The studio survives a scene effect restart or temporary catalogue loading.
// A fresh document/studio has no checkpoint. Keep this small cache itself
// across edits to this helper as well as ordinary React Fast Refresh updates.
const hot=(import.meta as ImportMeta&{hot?:{data:Record<string,unknown>}}).hot;
const checkpoints=(hot?.data.sceneCameraCheckpoints as WeakMap<object,SceneCameraRestart>|undefined)??new WeakMap<object,SceneCameraRestart>();
if(hot)hot.data.sceneCameraCheckpoints=checkpoints;

export function sceneCameraIntent(model:string,navigation:number,state:{view:string;reset:number;focus?:number;camera?:number[]},transition=0):SceneCameraIntent{
 return {model,navigation,view:state.view,reset:state.reset,focus:state.focus??0,camera:state.camera,transition};
}
export function captureRawSceneCamera(camera:PerspectiveCamera,target:Vector3):RawSceneCamera{
 return {...captureExplosionView(camera,target),aspect:camera.aspect,near:camera.near,far:camera.far,filmGauge:camera.filmGauge,filmOffset:camera.filmOffset};
}
export function restoreRawSceneCamera(camera:PerspectiveCamera,target:Vector3,pose:RawSceneCamera):void{
 camera.aspect=pose.aspect;camera.near=pose.near;camera.far=pose.far;camera.filmGauge=pose.filmGauge;camera.filmOffset=pose.filmOffset;
 restoreExplosionView(camera,target,pose);
}
export function rememberSceneCamera(owner:object,restart:SceneCameraRestart):void{checkpoints.set(owner,restart);}
export function resumeSceneCamera(owner:object,intent:SceneCameraIntent):SceneCameraRestart|undefined{
 const saved=checkpoints.get(owner),previous=saved?.intent;
 // Clearing a saved frame/focus is a display change. A nonempty new frame,
 // positive focus request, view/reset change or URL navigation is explicit.
 if(!saved||!previous||intent.model!==previous.model||intent.navigation!==previous.navigation
  ||intent.view!==previous.view||intent.reset!==previous.reset
  ||!!intent.focus&&intent.focus!==previous.focus||!!intent.camera&&intent.camera!==previous.camera
  ||!!intent.transition&&intent.transition!==previous.transition){checkpoints.delete(owner);return;}
 return saved;
}
