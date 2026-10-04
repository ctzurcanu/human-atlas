import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

// The combined source is not a partition into TA98 bones. Connections refer
// to the named anatomical route, not certified surface contact or joint tissue.
// NASEM NBK557988 (TMJ); anatomical study PMID18769090 (C0–C1).
export function isNativeFemaleCranialSource(part:Part):boolean{
 return !part.suppressed&&part.system==='skeletal'&&part.conceptId==='LOCAL:female:skull_cranium'
  &&/^LOCAL:female:skull_cranium:(999|1000|1001)$/.test(part.id)&&part.name==='Skull cranium';
}
export function cranialSourceLimitation(part:Part):string|undefined{
 if(isNativeFemaleCranialSource(part))return 'Combined cranial source placed under cranium, without separate bone identities. Its three material pieces and seven connected components are not anatomical subdivisions. The original source retains 388 detected triangle intersection pairs, nonmanifold edges and degenerate faces. Bone layers, sutures, foraminal passages and joint interfaces remain unverified; local repair trials are not installed.';
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function cranialRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const previous=cache.get(parts);if(previous)return previous;
 const graph=new Map<string,ResolvedRelation[]>();
 const targets=[
  {id:'LOCAL:female:mandible:960',concept:'LOCAL:female:mandible',name:'Mandible',note:'Mandibular condyles articulate with the temporal bones through the temporomandibular joints. Here the temporal regions remain within the combined cranial source; separate articular discs, cartilage, capsules and physical joint contacts are not verified.'},
  {id:'LOCAL:female:atlas_vertebra_c1:1006',concept:'LOCAL:female:atlas_vertebra_c1',name:'Atlas vertebra c1',note:'Occipital condyles articulate with the superior articular facets of the atlas at the atlanto-occipital joints. Here the occipital region remains within the combined cranial source; separate joint tissues, foraminal clearance and physical articular contacts are not verified.'},
 ];
 for(const source of parts.filter(isNativeFemaleCranialSource))for(const rule of targets){
  const target=parts.find(p=>!p.suppressed&&p.system==='skeletal'&&p.id===rule.id&&p.conceptId===rule.concept&&p.name===rule.name&&nativePartFrame(p)===nativePartFrame(source));
  if(!target)continue;
  for(const [a,b] of [[source,target],[target,source]]){
   const links=graph.get(a.id)??[];links.push({kind:'articulates',target:b,note:rule.note});graph.set(a.id,links);
  }
 }
 cache.set(parts,graph);return graph;
}
