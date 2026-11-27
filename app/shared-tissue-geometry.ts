import {BufferGeometry,InstancedBufferGeometry,BufferAttribute,InstancedBufferAttribute} from 'three';
import type {Part} from './anatomy';
import {partIndexBuffer} from './part-index-buffer';

/** An individually drawn bank view needs one selection ID, not one per
 * vertex of every other cell in its shared position bank. A single instance
 * supplies the constant attribute without changing indices or coordinates.
 */
export function anatomyPartGeometry(part:Pick<Part,'geometryOrigin'|'vertexCount'>,index:number,total:number):BufferGeometry{
 const geometry=part.geometryOrigin?new InstancedBufferGeometry():new BufferGeometry();
 if(geometry instanceof InstancedBufferGeometry){
  geometry.instanceCount=1;
  geometry.setAttribute('partIndex',new InstancedBufferAttribute(partIndexBuffer(index,1,total),1));
 }else geometry.setAttribute('partIndex',new BufferAttribute(partIndexBuffer(index,part.vertexCount,total),1));
 return geometry;
}

/** Display and picking views retain the very same attribute/index objects. */
export function anatomyGeometryView(source:BufferGeometry):BufferGeometry{
 const geometry=source instanceof InstancedBufferGeometry?new InstancedBufferGeometry():new BufferGeometry();
 if(geometry instanceof InstancedBufferGeometry)geometry.instanceCount=(source as InstancedBufferGeometry).instanceCount;
 for(const [name,attribute]of Object.entries(source.attributes))geometry.setAttribute(name,attribute);
 geometry.setIndex(source.getIndex());geometry.boundingBox=source.boundingBox;geometry.boundingSphere=source.boundingSphere;
 return geometry;
}
