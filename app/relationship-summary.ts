import type {ResolvedRelation} from './anatomical-relations';

/** A shared target can have different anatomical roles for selected pieces. */
export function mergeRelationshipNotes(relations:Iterable<ResolvedRelation>):ResolvedRelation[]{
 const merged=new Map<string,{relation:ResolvedRelation;notes:Set<string>;routes:Map<string,{via:string[];modeled:boolean}>}>();
 for(const relation of relations){
  const key=`${relation.kind}:${relation.target.id}`;
  let entry=merged.get(key);
  if(!entry){entry={relation:{...relation},notes:new Set(),routes:new Map()};merged.set(key,entry);}
  if(relation.note)entry.notes.add(relation.note);
  for(const route of relation.viaRoutes??(relation.via?.length?[{via:relation.via,modeled:!!relation.viaModeled}]:[]))entry.routes.set(JSON.stringify(route),{via:[...route.via],modeled:route.modeled});
 }
 return [...merged.values()].map(({relation,notes,routes})=>{
  const result=notes.size?{...relation,note:[...notes].join('\n')}:relation;
  if(routes.size>1){delete result.via;delete result.viaModeled;result.viaRoutes=[...routes.values()];}
  else if(routes.size===1){const route=[...routes.values()][0];delete result.viaRoutes;result.via=route.via;if(route.modeled)result.viaModeled=true;else delete result.viaModeled;}
  return result;
 });
}
