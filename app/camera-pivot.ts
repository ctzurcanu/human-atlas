import {Box3,Vector3} from 'three';
import type {Part} from './anatomy';

/** Include the current exploded translation when framing a selection. */
export function selectionCenter(parts:Part[],selectedIds:string[],translations:Float32Array){
 const selected=new Set(selectedIds),box=new Box3();
 parts.forEach((part,i)=>{if(selected.has(part.id))box.union(new Box3(new Vector3().fromArray(part.bounds[0]),new Vector3().fromArray(part.bounds[1])).translate(new Vector3(translations[i*4],translations[i*4+1],translations[i*4+2])));});
 return box.isEmpty()?null:box.getCenter(new Vector3());
}
