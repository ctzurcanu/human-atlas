import {PerspectiveCamera,Quaternion,Spherical,Vector2,Vector3} from 'three';

/** Rotate around the viewing target without an upright axis or polar stops. */
export function freeOrbit(camera:PerspectiveCamera,target:Vector3,dx:number,dy:number,viewportHeight:number){
 const distance=Math.hypot(dx,dy);
 if(!distance||viewportHeight<=0)return;
 const axis=new Vector3(-dy,-dx,0).normalize().applyQuaternion(camera.quaternion);
 const rotation=new Quaternion().setFromAxisAngle(axis,2*Math.PI*distance/viewportHeight);
 camera.position.sub(target).applyQuaternion(rotation).add(target);
 camera.up.applyQuaternion(rotation).normalize();
 camera.lookAt(target);
 camera.updateMatrixWorld(true);
}

/** Free orbit about an anatomical point, with no re-aiming jump at drag start. */
export function freeOrbitAroundPivot(camera:PerspectiveCamera,target:Vector3,pivot:Vector3,dx:number,dy:number,viewportHeight:number){
 const distance=Math.hypot(dx,dy);if(!distance||viewportHeight<=0)return;
 const axis=new Vector3(-dy,-dx,0).normalize().applyQuaternion(camera.quaternion);
 const rotation=new Quaternion().setFromAxisAngle(axis,2*Math.PI*distance/viewportHeight);
 camera.position.sub(pivot).applyQuaternion(rotation).add(pivot);
 target.sub(pivot).applyQuaternion(rotation).add(pivot);
 camera.up.applyQuaternion(rotation).normalize();camera.lookAt(target);camera.updateMatrixWorld(true);
}

/** Two-finger screen-space pan and pinch, preserving a freely rotated view. */
export function panPinch(camera:PerspectiveCamera,target:Vector3,before:[Vector2,Vector2],after:[Vector2,Vector2],viewportHeight:number,minDistance:number,maxDistance:number){
 if(viewportHeight<=0)return;
 const offset=camera.position.clone().sub(target),distance=offset.length();
 if(!distance)return;
 const delta=after[0].clone().add(after[1]).sub(before[0]).sub(before[1]).multiplyScalar(.5);
 const unitsPerPixel=2*distance*Math.tan(camera.fov*Math.PI/360)/viewportHeight;
 const shift=new Vector3(-delta.x,delta.y,0).applyQuaternion(camera.quaternion).multiplyScalar(unitsPerPixel);
 const oldSpan=before[0].distanceTo(before[1]),newSpan=after[0].distanceTo(after[1]);
 const nextDistance=oldSpan>=1&&newSpan>=1?Math.max(minDistance,Math.min(maxDistance,distance*oldSpan/newSpan)):distance;
 target.add(shift);
 camera.position.copy(target).addScaledVector(offset.normalize(),nextDistance);
 camera.updateMatrixWorld(true);
}

/** Orbit with a fixed screen-space pivot while preserving the camera's up axis. */
export function orbitAroundPointer(camera:PerspectiveCamera,target:Vector3,pivot:Vector3,ndc:Vector2,dx:number,dy:number,viewportHeight:number,minPolar=.001,maxPolar=Math.PI-.001){
 const up=new Vector3().copy(camera.up).normalize();
 const toY=new Quaternion().setFromUnitVectors(up,new Vector3(0,1,0));
 const offset=new Vector3().subVectors(camera.position,target).applyQuaternion(toY);
 const spherical=new Spherical().setFromVector3(offset);
 spherical.theta-=2*Math.PI*dx/viewportHeight;
 spherical.phi=Math.max(minPolar,Math.min(maxPolar,spherical.phi-2*Math.PI*dy/viewportHeight));
 spherical.makeSafe();
 offset.setFromSpherical(spherical).applyQuaternion(toY.invert());
 camera.position.copy(target).add(offset);
 camera.lookAt(target);
 camera.updateMatrixWorld(true);

 // OrbitControls normally turns around its target. Translate that orbit in
 // the view plane so the surface point under the initial pointer stays put.
 const forward=camera.getWorldDirection(new Vector3());
 const ray=new Vector3(ndc.x,ndc.y,.5).unproject(camera).sub(camera.position).normalize();
 const depth=new Vector3().subVectors(pivot,camera.position).dot(forward);
 const projected=ray.dot(forward);
 if(depth>0&&projected>.001){
  const pointOnRay=camera.position.clone().addScaledVector(ray,depth/projected);
  const shift=pivot.clone().sub(pointOnRay);
  camera.position.add(shift);
  target.add(shift);
  camera.lookAt(target);
  camera.updateMatrixWorld(true);
 }
}
