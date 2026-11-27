import {liverCompletionRegistry} from './liver-completion-controls';
import type {Atlas} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {liverFeatureRegistry} from './liver-surface-features';
import {liverLandmarkRegistry} from './liver-landmark';
type Reference={choice:HierarchyChoice;host:HierarchyChoice};
const cache=new WeakMap<Atlas,Map<string,Reference>>();
export function liverReferenceChoices(atlas:Atlas){
 let result=cache.get(atlas);if(result)return result;
 result=new Map();const concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const groups=[...(liverFeatureRegistry(atlas)?.termGroups??[]),...(liverLandmarkRegistry(atlas)?.marker.termGroups??[]),...(liverCompletionRegistry(atlas)?.termGroups.filter(g=>g.curves)??[])];
 for(const g of groups){const c=concepts.get(g.conceptId),p=atlas.parts.find(p=>p.id===g.partIds[0]);
  if(c&&p)result.set(g.term,{choice:c,host:{id:p.conceptId,name:p.name,elements:[p.id]}});
 }
 cache.set(atlas,result);return result;
}
export function liverReferenceActive(reference:Reference,choice:HierarchyChoice|undefined){
 return !!choice&&(choice.id===reference.choice.id||(choice.ta98Term??choice.terminology?.ta98)===reference.choice.ta98Term)
  &&choice.elements.length===reference.choice.elements.length&&choice.elements.every(id=>reference.choice.elements.includes(id));
}
/** A reference checkbox switches selection ink. Clearing it selects the same
 * host, preserving physical visibility, sections and every camera field. */
export function toggleLiverReferenceChoice(reference:Reference,choice:HierarchyChoice|undefined){
 return liverReferenceActive(reference,choice)?reference.host:reference.choice;
}
