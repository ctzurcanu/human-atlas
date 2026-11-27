import * as T from 'three';
import type {LiverReferenceGroup} from './liver-completion-controls';
export function createLiverBoundaryOverlay(parent:T.Group,planes:T.Plane[]){
 const material=new T.LineBasicMaterial({color:0x9ae8ee,depthWrite:false,transparent:true,opacity:1,clippingPlanes:planes});
 const lines=new T.LineSegments(new T.BufferGeometry(),material);lines.visible=false;lines.matrixAutoUpdate=false;lines.renderOrder=4;lines.frustumCulled=false;parent.add(lines);
 let key='';
 return {update(group:LiverReferenceGroup|undefined,host:T.Mesh|undefined,origin:readonly number[]=[0,0,0]){
  lines.visible=!!group&&!!host;if(!group||!host)return;
  if(group.conceptId!==key){key=group.conceptId;const points:number[]=[];
   for(const curve of group.curves??[])for(let i=0;i<curve.witnesses.length-(curve.closed?0:1);i++)for(const witness of [curve.witnesses[i],curve.witnesses[(i+1)%curve.witnesses.length]])points.push(...witness.point.map((x,k)=>x-origin[k]));
   lines.geometry.dispose();lines.geometry=new T.BufferGeometry().setAttribute('position',new T.Float32BufferAttribute(points,3));
  }
  lines.matrix.copy(host.matrixWorld);lines.matrixWorldNeedsUpdate=true;
 },dispose(){lines.geometry.dispose();material.dispose();lines.removeFromParent();}};
}
