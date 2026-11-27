import type {Atlas,Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';
import type {WallPeelStep} from './bowel-wall-peeling';
import {nativePartFrame} from './ta98-modeled-parent';

export const vascularCoatKeys=['externa','media','intima'] as const;
export const vascularCoatTerms=['A12.0.00.017','A12.0.00.019','A12.0.00.018'] as const;
const configurations=[
 {sourceId:'ZA:Hepatic portal vein',sourceName:'Hepatic portal vein',system:'venous',assemblyId:'portal:male:native-trunk',wholeTerm:'A12.3.12.001',prefix:'VEINWALL:male:portal',conceptSuffix:'',chart:'tissueChart'},
 {sourceId:'ZA:Thoracic aorta',sourceName:'Thoracic aorta',system:'arterial',assemblyId:'thoracic-aorta:male:native-trunk',wholeTerm:'A12.2.11.001',prefix:'TA98WALL:male:thoracic-aorta',conceptSuffix:':concept',chart:'uvs'},
] as const;
type Configuration=typeof configurations[number];
type Registry={source:Part;members:Part[];ids:string[];reference:Part[];configuration:Configuration};
const registries=new WeakMap<Part[],Registry[]>();

/** Only explicitly registered native replacements establish finite coats. */
export function vascularWallRegistries(parts:Part[]):Registry[]{
 const cached=registries.get(parts);if(cached)return cached;
 const result:Registry[]=[];
 for(const configuration of configurations){
  const {sourceId,sourceName,system,assemblyId,conceptSuffix,chart}=configuration;
  const ids=vascularCoatKeys.map(key=>`${configuration.prefix}:${key}`);
  const source=parts.find(p=>p.id===sourceId),members=ids.map(id=>parts.find(p=>p.id===id)),first=members[0];
  const valid=source?.suppressed&&source.conceptId===sourceId&&source.sourceId===sourceName
   &&source.system===system&&!source.sourceOffset&&source.replacementParts?.length===3&&ids.every(id=>source.replacementParts!.includes(id))
   &&members.every((p,k)=>p&&!p.suppressed&&p.system===system&&p.conceptId===p.id+conceptSuffix&&p.sourceId===sourceName
    &&p.pairedTissue&&p.sectionSolid&&!p.sourceOffset&&p.tissueWallAssembly===assemblyId&&p.tissueWallLayer===vascularCoatKeys[k]
    &&nativePartFrame(p)==='ZA:'&&p.normalSigns!==undefined&&p[chart]!==undefined&&p.tissueDrawGroups?.length
    &&p.geometryOrigin?.length===3&&p.geometryOrigin.every(Number.isFinite)&&p.chunk===first!.chunk&&p.positions===first!.positions
    &&p.normals===first!.normals&&p[chart]===first![chart]&&p.vertexCount===first!.vertexCount&&JSON.stringify(p.geometryOrigin)===JSON.stringify(first!.geometryOrigin))
   &&parts.filter(p=>!p.suppressed&&p.tissueWallAssembly===assemblyId).length===3;
  if(valid)result.push({source:source!,members:members as Part[],ids,configuration,reference:parts.map(p=>p===source?{...p,suppressed:false}:p)});
 }
 registries.set(parts,result);return result;
}
/** Compatibility accessor; choice-sensitive callers use the complete registry. */
export function vascularWallRegistry(parts:Part[]):Registry|null{return vascularWallRegistries(parts)[0]??null;}
const sameMembers=(elements:readonly string[],r:Registry)=>elements.length===r.ids.length&&new Set(elements).size===r.ids.length&&r.ids.every(id=>elements.includes(id));
export function vascularWallIntegrations(atlas:Atlas):Registry[]{
 if(atlas.sex!=='male')return [];
 return vascularWallRegistries(atlas.parts).filter(r=>{
  const whole=atlas.concepts.find(c=>c.id===r.source.conceptId);
  if(whole?.ta98Term!==r.configuration.wholeTerm||whole.ta98Kind!=='exact'||!sameMembers(whole.elements,r)
   ||r.members.some((p,k)=>{const c=atlas.concepts.find(c=>c.id===p.conceptId);return c?.partOf!==whole.id||c.ta98Term!==vascularCoatTerms[k]||c.ta98Kind!=='exact'||c.elements.length!==1||c.elements[0]!==p.id;}))return false;
  return atlas.vascularWallIntegration?.some(marker=>marker.sex==='male'&&marker.sourcePartId===r.source.id
   &&marker.wholeConceptId===r.source.conceptId&&marker.assemblyId===r.configuration.assemblyId&&sameMembers(marker.partIds,r))??false;
 });
}
export function vascularWallIntegration(atlas:Atlas):Registry|null{return vascularWallIntegrations(atlas)[0]??null;}
const componentRegistries=new WeakMap<Map<string,Part>,Registry[]>();
export function vascularWallComponentIds(choice:HierarchyChoice,parts:Map<string,Part>):string[]|undefined{
 if(!componentRegistries.has(parts))componentRegistries.set(parts,vascularWallRegistries([...parts.values()]));
 const r=choice.id!=='selection-set'&&componentRegistries.get(parts)!.find(r=>sameMembers(choice.elements,r));
 return r?[...r.ids]:undefined;
}
export function vascularWallPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const r=choice.id!=='selection-set'&&vascularWallIntegrations(atlas).find(r=>sameMembers(choice.elements,r));
 return r?r.members.map((p,k)=>({name:`Tunica ${vascularCoatKeys[k]} of ${r.source.name.toLowerCase()}`,ids:[p.id]})):null;
}
type Resolver=(source:Part,parts:Part[])=>ResolvedRelation[];
/** Flow belongs to the whole registered trunk, not every histological coat. */
export function vascularWallWholeRelations(concept:Concept,selected:Part[],parts:Part[],resolve:Resolver):ResolvedRelation[]|null{
 const r=vascularWallRegistries(parts).find(r=>sameMembers(concept.elements,r));
 if(!r||concept.id==='selection-set'||!sameMembers(concept.elements,r)||!sameMembers(selected.map(p=>p.id),r))return null;
 return resolve(r.reference.find(p=>p.id===r.source.id)!,r.reference).filter(row=>['before','after'].includes(row.kind)&&!r.ids.includes(row.target.id));
}
const graphs=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function vascularWallRelations(parts:Part[],resolve:Resolver):Map<string,ResolvedRelation[]>{
 const prior=graphs.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 // The semantic reference has an active original source, so recursive source
 // resolution cannot identify it as a replacement registry.
 for(const r of vascularWallRegistries(parts)){
  const whole:Concept={id:r.source.conceptId,name:r.source.name,elements:r.ids,ta98Term:r.configuration.wholeTerm,ta98Kind:'exact'};
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
