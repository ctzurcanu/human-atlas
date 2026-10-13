import {Box3,Sphere,Vector3} from 'three';
import type {Part} from './anatomy';
/** A shared position bank includes other tissue interfaces. The generated
 * cell's stored world bounds cover its indexed vertices, so its conservative
 * picking sphere needs no scan of those unrelated interfaces.
 * Ordinary/donor anatomy retains BufferGeometry's usual calculation.
 */
export function bowelWallBoundingSphere(part:Part):Sphere|null{
 if(!part.bowelWallSource||!part.pairedTissue||!part.geometryOrigin||!part.tissueDrawGroups||part.sourceOffset)return null;
 const [lo,hi]=part.bounds,origin=part.geometryOrigin;
 if(lo.length!==3||hi.length!==3||origin.length!==3||[...lo,...hi,...origin].some(v=>!Number.isFinite(v))||lo.some((v,d)=>v>hi[d]))return null;
 const box=new Box3(new Vector3(...lo).sub(new Vector3(...origin)),new Vector3(...hi).sub(new Vector3(...origin)));
 const sphere=box.getBoundingSphere(new Sphere());
 // Float32 carrier coordinates and world-origin subtraction need a small
 // conservative allowance at a corner, especially on thin crest partitions.
 sphere.radius+=Math.max(1e-9,sphere.radius*1e-7);
 return sphere;
}
