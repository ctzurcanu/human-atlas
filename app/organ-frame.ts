import {Box3,Matrix4,Quaternion,Vector3} from 'three';
import type {Part} from './anatomy';
type Instanced=Pick<Part,'geometryOrigin'|'geometryInstance'>;
const identity=new Quaternion(),poses=new WeakMap<object,Quaternion>();
export function organOrientation(part:Pick<Part,'geometryInstance'>){
 if(!part.geometryInstance)return identity;
 const known=poses.get(part);if(known)return known;
 const {rotation,scale}=part.geometryInstance;
 if(rotation.length!==4||!rotation.every(Number.isFinite)||!Number.isFinite(scale)||scale<=0||Math.abs(rotation.reduce((s,v)=>s+v*v,0)-1)>1e-8)throw new Error('Invalid anatomy instance placement.');
 const q=new Quaternion().fromArray(rotation).normalize();poses.set(part,q);return q;
}
export const organInstanceScale=(part:Pick<Part,'geometryInstance'>)=>part.geometryInstance?.scale??1;
export function organAnatomicalPoint(part:Instanced,point:Vector3){return point.clone().multiplyScalar(organInstanceScale(part)).applyQuaternion(organOrientation(part)).add(organOrigin(part));}
export function organGeometryPoint(part:Instanced,point:Vector3){return point.clone().sub(organOrigin(part)).applyQuaternion(organOrientation(part).clone().invert()).divideScalar(organInstanceScale(part));}
/** Coordinates stay near their organ; body bounds remain in the atlas frame. */
export function organOrigin(part:Pick<Part,'geometryOrigin'>){return new Vector3().fromArray(part.geometryOrigin??[0,0,0]);}
export function organBounds(part:Pick<Part,'bounds'|'geometryOrigin'|'geometryInstance'>,matrix?:Matrix4){
 const box=new Box3(new Vector3(...part.bounds[0]),new Vector3(...part.bounds[1]));
 if(matrix){
  if(part.geometryInstance){const assembled=new Matrix4().compose(organOrigin(part),organOrientation(part),new Vector3().setScalar(organInstanceScale(part)));box.applyMatrix4(matrix.clone().multiply(assembled.invert()));}
  else box.translate(organOrigin(part).negate()).applyMatrix4(matrix);
 }
 return box;
}
export function organPose(part:Instanced,offset:Vector3,rotation=new Quaternion(),pivot=new Vector3()){
 const position=organOrigin(part).add(offset).sub(pivot).applyQuaternion(rotation).add(pivot);
 return new Matrix4().compose(position,rotation.clone().multiply(organOrientation(part)),new Vector3().setScalar(organInstanceScale(part)));
}
