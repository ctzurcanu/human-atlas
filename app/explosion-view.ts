import {Box3,PerspectiveCamera,Vector3} from 'three';
import {zoomOutToFit,type ViewArea} from './view-framing';

/** A lossless view checkpoint, independent of anatomy bounds and URL framing. */
export function captureExplosionView(camera:PerspectiveCamera,target:Vector3){
 return {position:camera.position.clone(),target:target.clone(),up:camera.up.clone(),quaternion:camera.quaternion.clone(),zoom:camera.zoom,fov:camera.fov,view:camera.view?{...camera.view}:null};
}
export type ExplosionView=ReturnType<typeof captureExplosionView>;
export function restoreExplosionView(camera:PerspectiveCamera,target:Vector3,view:ExplosionView){
 camera.position.copy(view.position);target.copy(view.target);camera.up.copy(view.up);camera.quaternion.copy(view.quaternion);camera.zoom=view.zoom;camera.fov=view.fov;camera.view=view.view?{...view.view}:null;camera.updateProjectionMatrix();camera.updateMatrixWorld(true);
}
/** Fit only the current moving extent, starting from the original user view.
 * The assembled endpoint always restores that view, even if it was cropped. */
export function fitExplosionView(camera:PerspectiveCamera,target:Vector3,view:ExplosionView,bounds:Box3,area:ViewArea,amount:number){
 restoreExplosionView(camera,target,view);
 if(amount>0)zoomOutToFit(camera,bounds,area);
}
