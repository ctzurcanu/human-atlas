import * as T from 'three';
const ray=new T.Raycaster(),direction=new T.Vector3(.371,.529,.763).normalize();
/** Parity against the actual closed tissue boundary. Empty cap support quads
 * are not anatomy; a procedural source field may differ from a retained LOD.
 * Call with a DoubleSide material so both entering and leaving faces count.
 */
export function insideClosedTissue(point:T.Vector3,mesh:T.Mesh){
 ray.set(point,direction);const hits=ray.intersectObject(mesh,false);let crossings=0,last=-Infinity;
 for(const hit of hits){if(hit.distance<=1e-12||hit.distance-last<=1e-12)continue;crossings++;last=hit.distance;}
 return crossings%2===1;
}
