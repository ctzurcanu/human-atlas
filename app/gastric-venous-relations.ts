import type {Part} from './anatomy';
import type {RelationKind,ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const native=(p:Part,id:string,concept:string,source:string,name:string,system:Part['system'])=>
 !p.suppressed&&!p.sourceOffset&&p.id===id&&p.conceptId===concept&&p.sourceId===source&&p.name===name
 &&p.system===system&&nativePartFrame(p)==='LOCAL:female:';
export const isNativeFemaleGastricVein=(p:Part)=>native(p,'LOCAL:female:l_gastric_vein:1793','LOCAL:female:l_gastric_vein','l_gastric_vein','Gastric vein (left)','venous');
const gastricArtery=(p:Part)=>native(p,'LOCAL:female:l_gastric_artery:1720','LOCAL:female:l_gastric_artery','l_gastric_artery','Gastric artery (left)','arterial');
const portalVein=(p:Part)=>native(p,'LOCAL:female:portal_vein:1795','LOCAL:female:portal_vein','portal_vein','Portal vein','venous');

/** A reviewed native artery–vein context is distinct from arterial branching.
 * Do not let this exception admit matching donor names or another source frame. */
export function isReviewedGastricVenousRelation(kind:RelationKind,a:Part,b:Part):boolean{
 return kind==='adjacent'&&(isNativeFemaleGastricVein(a)&&gastricArtery(b)||isNativeFemaleGastricVein(b)&&gastricArtery(a));
}

const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** IFAA A12.3.12.015; Huang et al., PMID23793790: gastric venous routes vary.
 * These named native associations certify neither a complete vessel nor its
 * blood lumen, wall coats or an arterial–venous physical junction. */
export function gastricVenousRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const previous=cache.get(parts);if(previous)return previous;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(a:Part,b:Part,kind:RelationKind,note:string)=>{const rows=graph.get(a.id)??[];rows.push({kind,target:b,note});graph.set(a.id,rows);};
 for(const vein of parts.filter(isNativeFemaleGastricVein)){
  for(const portal of parts.filter(portalVein)){
   const note='Named gastric venous outflow to the native portal source. Shared exterior rims do not establish a patent blood lumen, vessel wall coats or complete tributary anatomy.';
   add(vein,portal,'after',note);add(portal,vein,'before',note);
  }
  for(const artery of parts.filter(gastricArtery)){
   const note='Native left gastric arterial and venous context along the lesser curvature. Their courses vary; this association does not establish a physical artery–vein junction.';
   add(vein,artery,'adjacent',note);add(artery,vein,'adjacent',note);
  }
 }
 cache.set(parts,graph);return graph;
}
