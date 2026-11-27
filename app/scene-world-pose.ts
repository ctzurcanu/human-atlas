import type {PerspectiveCamera,Vector3} from 'three';
import type {ScenePose,SceneWorldPose} from './connect-state';
import {validCamera} from '../shared/camera-frame.mjs';
import {fitViewToAspect} from './view-framing';

const object=(value:unknown):value is Record<string,unknown>=>!!value&&typeof value==='object'&&!Array.isArray(value);
export function validSceneWorldPose(value:unknown):value is SceneWorldPose{
 if(!object(value))return false;
 for(const [key,size]of [['position',3],['target',3],['up',3],['quaternion',4]]as const){const list=value[key];if(!Array.isArray(list)||list.length!==size||!list.every(item=>typeof item==='number'&&Number.isFinite(item)&&Math.abs(item)<10000))return false;}
 for(const key of ['fov','zoom','near','far'])if(typeof value[key]!=='number'||!Number.isFinite(value[key])||Number(value[key])<=0)return false;
 if(Number(value.fov)>=180||Number(value.far)<=Number(value.near)||Math.hypot(...value.up as number[])<1e-12||Math.hypot(...value.quaternion as number[])<1e-12||Math.hypot(...(value.position as number[]).map((n,i)=>n-(value.target as number[])[i]))<1e-12)return false;
 for(const key of ['aspect','filmGauge'])if(value[key]!==undefined&&(typeof value[key]!=='number'||!Number.isFinite(value[key])||Number(value[key])<=0))return false;
 if(value.filmOffset!==undefined&&(typeof value.filmOffset!=='number'||!Number.isFinite(value.filmOffset)))return false;
 if(value.view!==undefined&&value.view!==null){const v=value.view;if(!object(v)||typeof v.enabled!=='boolean'||!['fullWidth','fullHeight','offsetX','offsetY','width','height'].every(key=>typeof v[key]==='number'&&Number.isFinite(v[key])&&Math.abs(Number(v[key]))<1e7)||['fullWidth','fullHeight','width','height'].some(key=>Number(v[key])<=0))return false;}
 return true;
}
export function captureSceneWorldPose(camera:PerspectiveCamera,target:Vector3):SceneWorldPose{
 return {position:camera.position.toArray(),target:target.toArray(),up:camera.up.toArray(),quaternion:camera.quaternion.toArray(),fov:camera.fov,zoom:camera.zoom,near:camera.near,far:camera.far,aspect:camera.aspect,filmGauge:camera.filmGauge,filmOffset:camera.filmOffset,view:camera.view?{...camera.view}:null};
}
export function restoreSceneWorldPose(camera:PerspectiveCamera,target:Vector3,pose:SceneWorldPose){
 camera.position.fromArray(pose.position);target.fromArray(pose.target);camera.up.fromArray(pose.up);camera.quaternion.fromArray(pose.quaternion).normalize();camera.fov=pose.fov;camera.zoom=pose.zoom;camera.near=pose.near;camera.far=pose.far;
 // The aspect belongs to the canvas, not to the pose (a pose from another window size would stretch the body).
 if(pose.filmGauge!==undefined)camera.filmGauge=pose.filmGauge;if(pose.filmOffset!==undefined)camera.filmOffset=pose.filmOffset;
 if(pose.view!==undefined)camera.view=pose.view?{...pose.view}:null;fitViewToAspect(camera);camera.updateMatrixWorld(true);
}
export function writeWorldView(url:URL,model:string,pose:ScenePose|null){
 if(pose?.model===model&&pose.world&&validSceneWorldPose(pose.world)&&validCamera(pose.camera))url.searchParams.set('world',JSON.stringify(pose));
}
export function readWorldView(parameters:URLSearchParams,model:string):ScenePose|null{
 const raw=parameters.get('world');if(!raw||raw.length>6000)return null;
 try{const pose:unknown=JSON.parse(raw);if(!object(pose)||pose.model!==model||!validCamera(pose.camera)||!Array.isArray(pose.up)||pose.up.length!==3||!pose.up.every(Number.isFinite)||Math.hypot(...pose.up)<1e-12||!Array.isArray(pose.rotation)||pose.rotation.length!==4||!pose.rotation.every(Number.isFinite)||Math.hypot(...pose.rotation)<1e-12||!validSceneWorldPose(pose.world))return null;return pose as ScenePose;}catch{return null;}
}
