import type {Atlas,Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';
import type {WallPeelStep} from './bowel-wall-peeling';
import {nativePartFrame} from './ta98-modeled-parent';

const sourceId='ZA:Hepatic portal vein';
const assemblyId='portal:male:native-trunk';
export const vascularCoatKeys=['externa','media','intima'] as const;
export const vascularCoatTerms=['A12.0.00.017','A12.0.00.019','A12.0.00.018'] as const;
const ids=vascularCoatKeys.map(key=>`VEINWALL:male:portal:${key}`);
type Registry={source:Part;members:Part[];ids:string[];reference:Part[]};
const registries=new WeakMap<Part[],Registry|null>();

/** Only the explicit replacement of this registered native trunk is eligible.
 * Named tributaries, donors and generic vascular terms do not establish coats. */
export function vascularWallRegistry(parts:Part[]):Registry|null{
 if(registries.has(parts))return registries.get(parts)!;
 const source=parts.find(p=>p.id===sourceId),members=ids.map(id=>parts.find(p=>p.id===id)),first=members[0];
 const valid=source?.suppressed&&source.conceptId===sourceId&&source.sourceId==='Hepatic portal vein'
  &&source.system==='venous'&&!source.sourceOffset&&source.replacementParts?.length===3&&ids.every(id=>source.replacementParts!.includes(id))
  &&members.every((p,k)=>p&&!p.suppressed&&p.system==='venous'&&p.conceptId===p.id&&p.sourceId===source.sourceId
   &&p.pairedTissue&&p.sectionSolid&&!p.sourceOffset&&p.tissueWallAssembly===assemblyId&&p.tissueWallLayer===vascularCoatKeys[k]
   &&nativePartFrame(p)==='ZA:'&&p.normalSigns!==undefined&&p.tissueChart!==undefined&&p.tissueDrawGroups?.length
   &&p.geometryOrigin?.length===3&&p.geometryOrigin.every(Number.isFinite)&&p.chunk===first!.chunk&&p.positions===first!.positions
   &&p.normals===first!.normals&&p.vertexCount===first!.vertexCount&&JSON.stringify(p.geometryOrigin)===JSON.stringify(first!.geometryOrigin))
  &&parts.filter(p=>!p.suppressed&&p.tissueWallAssembly===assemblyId).length===3;
 const result=valid?{source:source!,members:members as Part[],ids,reference:parts.map(p=>p===source?{...p,suppressed:false}:p)}:null;
 registries.set(parts,result);return result;
}
const sameMembers=(elements:readonly string[],r:Registry)=>elements.length===r.ids.length&&new Set(elements).size===r.ids.length&&r.ids.every(id=>elements.includes(id));
export function vascularWallIntegration(atlas:Atlas):Registry|null{
 const r=vascularWallRegistry(atlas.parts);
 if(!r||atlas.sex!=='male')return null;
 const whole=atlas.concepts.find(c=>c.id===r.source.conceptId);
 if(whole?.ta98Term!=='A12.3.12.001'||whole.ta98Kind!=='exact'||!sameMembers(whole.elements,r)
  ||r.members.some((p,k)=>{const c=atlas.concepts.find(c=>c.id===p.conceptId);return c?.partOf!==whole.id||c.ta98Term!==vascularCoatTerms[k]||c.ta98Kind!=='exact'||c.elements.length!==1||c.elements[0]!==p.id;}))return null;
 return atlas.vascularWallIntegration?.some(marker=>marker.sex==='male'&&marker.sourcePartId===r.source.id
  &&marker.wholeConceptId===r.source.conceptId&&marker.assemblyId===assemblyId&&sameMembers(marker.partIds,r))?r:null;
}
const componentRegistries=new WeakMap<Map<string,Part>,Registry|null>();
export function vascularWallComponentIds(choice:HierarchyChoice,parts:Map<string,Part>):string[]|undefined{
 if(!componentRegistries.has(parts))componentRegistries.set(parts,vascularWallRegistry([...parts.values()]));
 const r=componentRegistries.get(parts);
 return r&&choice.id!=='selection-set'&&sameMembers(choice.elements,r)?[...r.ids]:undefined;
}
export function vascularWallPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const r=vascularWallIntegration(atlas);
 return r&&choice.id!=='selection-set'&&sameMembers(choice.elements,r)
  ?r.members.map((p,k)=>({name:`Tunica ${vascularCoatKeys[k]} of hepatic portal vein`,ids:[p.id]})):null;
}
type Resolver=(source:Part,parts:Part[])=>ResolvedRelation[];
/** Flow belongs to the whole registered trunk, not every histological coat. */
export function vascularWallWholeRelations(concept:Concept,selected:Part[],parts:Part[],resolve:Resolver):ResolvedRelation[]|null{
 const r=vascularWallRegistry(parts);
 if(!r||concept.id==='selection-set'||!sameMembers(concept.elements,r)||!sameMembers(selected.map(p=>p.id),r))return null;
 return resolve(r.reference.find(p=>p.id===r.source.id)!,r.reference).filter(row=>['before','after'].includes(row.kind)&&!r.ids.includes(row.target.id));
}
const graphs=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function vascularWallRelations(parts:Part[],resolve:Resolver):Map<string,ResolvedRelation[]>{
 const prior=graphs.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),r=vascularWallRegistry(parts);
 // The semantic reference has an active original source, so recursive source
 // resolution cannot identify it as a replacement registry.
 if(r){
  const whole:Concept={id:r.source.conceptId,name:r.source.name,elements:r.ids,ta98Term:'A12.3.12.001',ta98Kind:'exact'};
  const target={...r.members[0],name:r.source.name,conceptId:whole.id};
  for(const row of vascularWallWholeRelations(whole,r.members,parts,resolve)??[]){
   const kind=row.kind==='before'?'after':'before';
   const rows=graph.get(row.target.id)??[];rows.push({kind,target,targetConcept:whole,note:row.note});graph.set(row.target.id,rows);
  }
  for(const p of r.members)for(const group of p.tissueDrawGroups??[]){
   if(!group.neighbor)continue;const neighbor=r.members.find(q=>q.id===group.neighbor);
   if(!neighbor||!(neighbor.tissueDrawGroups??[]).some(g=>g.neighbor===p.id&&g.count===group.count))continue;
   const rows=graph.get(p.id)??[];
   if(!rows.some(row=>row.target.id===neighbor.id))rows.push({kind:'adjacent',target:neighbor,note:'Opposed shared interface between finite vascular coats. Thickness and histological proportions remain provisional.'});
   graph.set(p.id,rows);
  }
 }
 graphs.set(parts,graph);return graph;
}
