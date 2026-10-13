import {Box3,Matrix4,Quaternion,Vector3} from 'three';
import type {Part} from './anatomy';
/** Coordinates stay near their organ; body bounds remain in the atlas frame. */
export function organOrigin(part:Pick<Part,'geometryOrigin'>){return new Vector3().fromArray(part.geometryOrigin??[0,0,0]);}
export function organBounds(part:Pick<Part,'bounds'|'geometryOrigin'>,matrix?:Matrix4){
 const box=new Box3(new Vector3(...part.bounds[0]),new Vector3(...part.bounds[1]));
 if(matrix)box.translate(organOrigin(part).negate()).applyMatrix4(matrix);
 return box;
}
export function organPose(part:Pick<Part,'geometryOrigin'>,offset:Vector3,rotation=new Quaternion(),pivot=new Vector3()){
 const position=organOrigin(part).add(offset).sub(pivot).applyQuaternion(rotation).add(pivot);
 return new Matrix4().compose(position,rotation,new Vector3(1,1,1));
}
