import type {Atlas,Part,Concept} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';
import {mergeRelationshipNotes} from './relationship-summary';
import {spinalCordRegistry,type SpinalBoundaryReference} from './spinal-cord-controls';
import inventory from './spinal-cord-system-inventory';

const parents=new Map<string,string>(inventory.terms.map(row=>[row.id,row.parent]));
const descendant=(term:string,parent:string)=>{
 const seen=new Set<string>();
 while(parents.has(term)&&!seen.has(term)){seen.add(term);term=parents.get(term)!;if(term===parent)return true;}
 return false;
};

/** Mixed system/body branches must also count names, rather than thousands
 * of rendering fragments, whenever they include the shared spinal bank. */
export function spinalCordDisplayCount(parts:readonly Part[]):number|undefined{
 if(!parts.some(p=>{
  const boundary=(p as Part&{spinalBoundaryReference?:SpinalBoundaryReference}).spinalBoundaryReference;
  return /^TA98CORD:(male|female):/.test(p.id)&&!p.suppressed&&
   (/^spinal-cord:(male|female)$/.test(p.tissueWallAssembly??'')||boundary?.kind==='lumen-boundary'&&boundary.physicalOwnerRetained===true);
 }))return;
 return new Set(parts.filter(p=>!p.suppressed).map(p=>p.conceptId)).size;
}

/** A physical fragment never becomes a separate component of the whole cord.
 * Keep each installed term once, and retain complete scoped tissue ownership.
 */
export function spinalCordComponentChoices(atlas:Atlas,choice:HierarchyChoice|undefined):HierarchyChoice[]|undefined{
 const r=spinalCordRegistry(atlas);if(!r||!choice||choice.id==='selection-set'||!choice.elements.length)return;
 const selected=new Set(choice.elements),parts=new Map([...r.parts,...r.referenceParts].map(p=>[p.id,p]));
 const referenceOwners=new Map(r.referenceParts.map(p=>[p.id,(p as Part&{spinalBoundaryReference:{ownerPartId:string}}).spinalBoundaryReference.ownerPartId]));
 if(choice.elements.some(id=>!parts.has(id)))return;
 const concepts=new Map(atlas.concepts.map(c=>[c.id,c])),groups=r.marker.termGroups;
 const segment=r.marker.segments?.find(s=>s.conceptId===choice.id||groups.some(g=>`${g.conceptId}@${s.name}`===choice.id));
 const term=choice.ta98Term??choice.terminology?.ta98;
 let root=groups.find(g=>g.conceptId===choice.id||segment&&`${g.conceptId}@${segment.name}`===choice.id)
  ??groups.find(g=>g.term===term);
 const inside=(id:string)=>selected.has(id)||selected.has(referenceOwners.get(id)??'');
 const scoped=(ids:string[])=>segment?ids.filter(inside):ids;
 const make=(id:string,ids:string[]):HierarchyChoice|undefined=>{
  const c=concepts.get(id);if(!c||!ids.length)return;
  return {...c,id:segment?`${id}@${segment.name}`:id,name:segment?`${c.name} · ${segment.name}`:c.name,elements:ids};
 };
 const stored=concepts.get(choice.id),exactStored=!!stored&&stored.elements.length===selected.size&&stored.elements.every(id=>selected.has(id));
 if(!root&&!segment&&!exactStored)return;
 if(root&&(!scoped(root.partIds).length||scoped(root.partIds).some(id=>!inside(id)))){
  // A sided nucleus is a valid atomic selection, not its bilateral parent.
  if(!exactStored)return;
  root=undefined;
 }
 const result:HierarchyChoice[]=[];
 if(segment&&choice.id===segment.conceptId){
  for(const term of ['A14.1.02.201','A14.1.02.020','A14.1.02.019','A14.1.01.309']){
   const g=groups.find(g=>g.term===term),c=g&&make(g.conceptId,scoped(g.partIds));if(c)result.push(c);
  }
 }else if(root){
  const regionSegments=r.marker.segments?.filter(s=>s.parentTerm===root.term&&!segment)??[];
  if(regionSegments.length){
   for(const s of regionSegments){const c=concepts.get(s.conceptId);if(c&&s.partIds.every(id=>selected.has(id)))result.push({...c,elements:[...s.partIds]});}
  }else{
   for(const g of groups){
    if(g.term===root.term)continue;
    const ids=scoped(g.partIds);
    if(!ids.length||ids.some(id=>!inside(id))||!descendant(g.term,root.term))continue;
    const c=make(g.conceptId,ids);if(c)result.push(c);
   }
  }
 }
 // A bounded residual groups fragments by their existing anatomical concept.
 // It covers material outside named descendants, without one row per segment.
 const covered=new Set(result.flatMap(c=>c.elements)),remaining=new Map<string,string[]>();
 for(const id of choice.elements){
  if(covered.has(id))continue;const p=parts.get(id)!;
  if(p.conceptId===choice.id||segment&&`${p.conceptId}@${segment.name}`===choice.id)continue;
  const ids=remaining.get(p.conceptId)??[];ids.push(id);remaining.set(p.conceptId,ids);
 }
 for(const [id,ids]of remaining){const c=make(id,ids);if(c&&!(c.name===choice.name&&ids.length===choice.elements.length))result.push(c);}
 return result;
}

/** Relations to a tissue name select its participating pieces together. */
export function spinalCordRelations(relations:ResolvedRelation[],concepts:Concept[]):ResolvedRelation[]{
 const byId=new Map(concepts.map(c=>[c.id,c])),groups=new Map<string,ResolvedRelation[]>();
 for(const relation of relations){
  const target=relation.target,spinal=spinalCordDisplayCount([target])!==undefined;
  const key=spinal&&!relation.targetConcept?`${relation.kind}:${target.conceptId}`:`${relation.kind}:${target.id}`;
  const group=groups.get(key)??[];group.push(relation);groups.set(key,group);
 }
 return [...groups.values()].flatMap(group=>{
  const first=group[0],concept=byId.get(first.target.conceptId);
  if(!concept||first.targetConcept||spinalCordDisplayCount([first.target])===undefined)return group;
  const targetConcept={...concept,elements:[...new Set(group.map(r=>r.target.id))]};
  return mergeRelationshipNotes(group.map(r=>({...r,target:first.target,targetConcept})));
 });
}
