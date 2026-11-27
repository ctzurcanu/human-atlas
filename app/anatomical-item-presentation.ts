import type {Part,Concept} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';
import {structureLabelName} from './laterality';
import {mergeRelationshipNotes} from './relationship-summary';

/** Tessellation, chunks and regional slices are not additional anatomy names. */
export function namedAnatomyMembers(parts:Part[],concepts:Concept[]):HierarchyChoice[]{
 const byId=new Map(concepts.map(c=>[c.id,c])),groups=new Map<string,Part[]>();
 for(const part of parts){if(part.suppressed)continue;const members=groups.get(part.conceptId)??[];members.push(part);groups.set(part.conceptId,members);}
 return [...groups].map(([id,members])=>({...byId.get(id),id,name:byId.get(id)?.name??structureLabelName(members[0].name),elements:members.map(p=>p.id)}));
}

/** Keep the participating subset when several fragments share a relation. */
export function namedAnatomyRelations(relations:ResolvedRelation[],concepts:Concept[]):ResolvedRelation[]{
 const byId=new Map(concepts.map(c=>[c.id,c])),groups=new Map<string,ResolvedRelation[]>();
 for(const relation of relations){
  const key=`${relation.kind}:${relation.targetConcept?.id??relation.target.conceptId}`;
  const members=groups.get(key)??[];members.push(relation);groups.set(key,members);
 }
 return [...groups.values()].flatMap(group=>{
  const first=group[0],concept=first.targetConcept??byId.get(first.target.conceptId);
  if(!concept)return group;
  const targetConcept={...concept,elements:[...new Set(group.flatMap(r=>r.targetConcept?.elements??[r.target.id]))]};
  return mergeRelationshipNotes(group.map(r=>({...r,target:first.target,targetConcept})));
 });
}
