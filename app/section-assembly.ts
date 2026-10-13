import * as T from 'three';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';
import type {StencilCapCandidate} from './stencil-caps';
import {sectionTissue} from './section-tissue';
import {gastricRegionRim} from './gastric-region-rim';

/** Join visible tissue partitions once; moving the cut never rebuilds the mesh. */
export function createSectionAssemblyCache(){
 const cache=new Map<string,{geometry:T.BufferGeometry;part:StencilCapCandidate['part'];mode:'wall'|'solid';closedWall:boolean}>();
 const collect=(candidates:StencilCapCandidate[])=>{
  const groups=new Map<string,StencilCapCandidate[]>();
  for(const item of candidates){if(item.part.pairedTissue||!item.part.sectionAssembly)continue;
   // The female regions have independently verified paired walls and annular
   // tissue rims. Count separately moved regions in their own source frames.
   // Other assemblies still require their whole shared closure.
   const transform=item.part.sectionAssembly==='TA98:female-stomach-wall'?`:${item.matrix.elements.join(',')}`:'';
   const key=`${item.cutIndex}:${item.part.sectionAssembly}${transform}`;const group=groups.get(key)??[];group.push(item);groups.set(key,group);}
  const caps:StencilCapCandidate[]=candidates.filter(item=>item.part.pairedTissue).map(item=>({...item,mode:'solid',closedWall:true})),covered=new Set<string>(caps.map(item=>`${item.cutIndex}:${item.index}`)),used=new Set<string>();
  for(const group of groups.values()){
   const first=group[0],role=sectionTissue(first.part);
   // Different transforms or tissue owners must never share a stencil fill.
   if(group.some(item=>sectionTissue(item.part)!==role||item.mode!==first.mode||item.wallWidth!==first.wallWidth||item.matrix.elements.some((value,i)=>Math.abs(value-first.matrix.elements[i])>1e-9)))continue;
   const key=group.map(item=>`${item.index}:${item.geometry.uuid}`).join(',');let entry=cache.get(key);
   if(!entry){
    const fragments=group.map(({geometry})=>{const result=new T.BufferGeometry();result.setAttribute('position',geometry.getAttribute('position'));result.setIndex(geometry.getIndex());return result;});
    let geometry=mergeGeometries(fragments);for(const fragment of fragments)fragment.dispose();if(!geometry)continue;
    if(role==='gastric-wall'&&first.part.sectionAssembly==='TA98:female-stomach-wall'&&!closedWall(geometry)){
     const rim=gastricRegionRim(geometry);
     if(rim){const indexed=rim.clone();indexed.deleteAttribute('normal');indexed.setIndex(Array.from({length:rim.getAttribute('position').count},(_,i)=>i));const closed=mergeGeometries([geometry,indexed]);indexed.dispose();rim.dispose();if(closed){geometry.dispose();geometry=closed;}}
    }
    // A closed paired wall already contains tissue volume and an open lumen.
    // Insetting it again would subtract real wall tissue and leave only slivers.
    const pairedWall=role==='gastric-wall'&&['TA98:female-stomach-wall','TA98:gastric-tract-wall'].includes(first.part.sectionAssembly??'')&&closedWall(geometry),mode=pairedWall?'solid':first.mode??'solid';
    const bounds=new T.Box3();for(const {part} of group)bounds.union(new T.Box3(new T.Vector3(...part.bounds[0]),new T.Vector3(...part.bounds[1])));
    geometry.computeBoundingBox();geometry.computeBoundingSphere();
    entry={geometry,mode,closedWall:pairedWall,part:{...first.part,bounds:[bounds.min.toArray(),bounds.max.toArray()] as typeof first.part.bounds}};cache.set(key,entry);
   }
   used.add(key);caps.push({...first,...entry});for(const item of group)covered.add(`${item.cutIndex}:${item.index}`);
  }
  for(const [key,entry] of cache){if(!used.has(key)){entry.geometry.dispose();cache.delete(key);}}
  return {caps,covered};
 };
 return {collect,dispose:()=>{for(const {geometry} of cache.values())geometry.dispose();cache.clear();}};
}

function closedWall(geometry:T.BufferGeometry){
 const position=geometry.getAttribute('position'),indices=geometry.index;if(!indices)return false;
 const edges=new Map<string,{count:number;winding:number}>(),vertices=Array.from({length:position.count},(_,i)=>[position.getX(i),position.getY(i),position.getZ(i)].map(v=>Math.round(v*1e7)).join(','));
 for(let i=0;i<indices.count;i+=3){const face=[0,1,2].map(k=>vertices[indices.getX(i+k)]);
  for(let k=0;k<3;k++){const a=face[k],b=face[(k+1)%3],key=[a,b].sort().join('|'),entry=edges.get(key)??{count:0,winding:0};entry.count++;entry.winding+=a<b?1:-1;edges.set(key,entry);}
 }
 return edges.size>0&&[...edges.values()].every(edge=>edge.count===2&&edge.winding===0);
}
