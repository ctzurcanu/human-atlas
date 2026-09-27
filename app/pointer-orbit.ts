import {PerspectiveCamera,Quaternion,Spherical,Vector2,Vector3} from 'three';

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
