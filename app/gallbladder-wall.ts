import type {Atlas,Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';
import type {WallPeelStep} from './bowel-wall-peeling';
import {gallbladderRelations} from './gallbladder-relations';
import {biliaryJunctionInterface} from './biliary-junction';
import {digestiveDuctRelations} from './digestive-duct-relations';

const femaleSourceId='LOCAL:female:gallbladder:1030';
const maleSourceId='ZA:Gallbladder';
const keys=['covering','subserosa','muscular','mucosa'];
const femaleIds=keys.map(key=>'LOCAL:female:GALLWALL:'+key);
const maleIds=keys.map(key=>'TA98PART:ZA:Gallbladder:GALLWALL:'+key);
const names=['External covering (allocation unresolved)','Tela subserosa','Tunica muscularis','Tunica mucosa'];
type Registry={sex:'male'|'female';ids:string[];source:Part;members:Part[];reference:Part[]};
const registries=new WeakMap<Part[],Registry|null>();
function registry(parts:Part[]):Registry|null{
 if(registries.has(parts))return registries.get(parts)!;
 let result:Registry|null=null;
 for(const sex of ['female','male'] as const){
  const sourceId=sex==='female'?femaleSourceId:maleSourceId;
  const source=parts.find(p=>p.id===sourceId);
  const regional=source?.replacementParts?.length===16;
  const regionCodes=['A05.8.02.002','A05.8.02.004','A05.8.02.003','A05.8.02.005'];
  const ids=regional?keys.flatMap(key=>regionCodes.map(code=>(sex==='male'?'TA98PART:ZA:Gallbladder:GALLWALL:':'LOCAL:female:GALLWALL:')+key+':'+code)):sex==='female'?femaleIds:maleIds;
  const members=ids.map(id=>parts.find(p=>p.id===id)),first=members[0];
  const valid=source?.suppressed&&source.name==='Gallbladder'&&source.sourceId===(sex==='female'?'gallbladder':'Gallbladder')&&!source.sourceOffset
   &&source.replacementParts?.length===ids.length&&ids.every(id=>source.replacementParts!.includes(id))
   &&members.every((p,k)=>p&&!p.suppressed&&p.system==='digestive'&&p.pairedTissue&&p.sectionSolid&&!p.sourceOffset
    &&p.bowelWallSource===sourceId&&p.tissueWallAssembly===`gallbladder:${sex}`&&p.tissueWallLayer===keys[regional?Math.floor(k/4):k]
    &&p.geometryOrigin?.length===3&&p.geometryOrigin.every(Number.isFinite)
    &&p.chunk===first!.chunk&&p.positions===first!.positions&&p.normals===first!.normals
    &&JSON.stringify(p.geometryOrigin)===JSON.stringify(first!.geometryOrigin)&&p.tissueDrawGroups?.length);
  if(valid){result={sex,ids,source:source!,members:members as Part[],reference:parts.map(p=>p===source?{...p,suppressed:false}:p)};break;}
 }
 registries.set(parts,result);return result;
}
export function gallbladderWholeSource(selected:Part[],parts:Part[]):Part|undefined{
 const r=registry(parts),selectedIds=new Set(selected.map(p=>p.id));
 return r&&selectedIds.size===r.ids.length&&r.members.every(p=>selectedIds.has(p.id))?r.source:undefined;
}
/** Whole-organ associations use the original registered source identity.
 * The tissue allocation does not establish a local supply to every coat. */
export function gallbladderWholeRelations(concept:Concept,selected:Part[],parts:Part[]):ResolvedRelation[]|null{
 const r=registry(parts);if(!r||!gallbladderWholeSource(selected,parts))return null;
 const chosen=new Set(concept.elements);if(chosen.size!==r.ids.length||!r.ids.every(id=>chosen.has(id)))return null;
 return [...(gallbladderRelations(r.reference).get(r.source.id)??[]),...(digestiveDuctRelations(r.reference).get(r.source.id)??[])];
}
export function gallbladderWallPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const marker=atlas.gallbladderWallPreview,r=registry(atlas.parts);
 if(!marker||!r||marker.sourcePartId!==r.source.id||marker.sex!==r.sex||marker.partIds.length!==r.ids.length||!r.ids.every(id=>marker.partIds.includes(id))||choice.id==='selection-set')return null;
 const selected=new Set(choice.elements);
 if(selected.size&&[...selected].every(id=>r.ids.includes(id))){
  const steps=keys.map((key,k)=>({name:names[k]+' of gallbladder',ids:r.members.filter(p=>p.tissueWallLayer===key&&selected.has(p.id)).map(p=>p.id)})).filter(step=>step.ids.length);
  const complete=selected.size===r.ids.length;
  const oneCoat=steps.length===1;
  const regional=r.source.replacementParts?.length===16&&selected.size===4&&steps.length===4&&steps.every(step=>step.ids.length===1)&&new Set([...selected].map(id=>id.split(':').at(-1))).size===1;
  if(complete||oneCoat||regional)return steps;
 }
 return null;
}
const graphs=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function gallbladderWallInterface(a:Part,b:Part):boolean{
 const femaleCoat=(p:Part)=>(femaleIds.includes(p.id)||p.id===`LOCAL:female:GALLWALL:${p.tissueWallLayer}:A05.8.02.005`)&&p.bowelWallSource===femaleSourceId&&p.tissueWallAssembly==='gallbladder:female'&&keys.includes(p.tissueWallLayer??'');
 const femaleDuct=(p:Part)=>p.id==='LOCAL:female:bile_duct:1031'&&p.sourceId==='bile_duct'&&p.conceptId==='LOCAL:female:bile_duct'||p.id==='LOCAL:female:cystic_duct:authored'&&p.sourceId==='cystic_duct'&&p.conceptId==='LOCAL:female:cystic_duct'&&p.name==='Cystic duct'&&p.tissueWallAssembly==='cystic:female'&&p.sectionAssembly==='gallbladder-cystic:female:folds';
 const maleCoat=(p:Part)=>p.id===`TA98PART:ZA:Gallbladder:GALLWALL:${p.tissueWallLayer}:A05.8.02.005`&&keys.includes(p.tissueWallLayer??'')&&p.bowelWallSource===maleSourceId&&p.tissueWallAssembly==='gallbladder:male';
 const maleDuct=(p:Part)=>p.id==='FJ3080'&&p.sourceId==='FJ3080'&&p.conceptId==='FMA14539'&&p.tissueWallAssembly==='cystic:male'&&p.sectionAssembly==='gallbladder-cystic:male:folds';
 const matches=(x:Part,y:Part)=>femaleCoat(x)&&femaleDuct(y)||maleCoat(x)&&maleDuct(y);
 if(!(matches(a,b)||matches(b,a))||[a,b].some(p=>p.suppressed||p.system!=='digestive'||!p.pairedTissue||!p.sectionSolid||p.sourceOffset)
  ||a.chunk!==b.chunk||a.positions!==b.positions||a.normals!==b.normals||JSON.stringify(a.geometryOrigin)!==JSON.stringify(b.geometryOrigin))return false;
 const count=(p:Part,id:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id).reduce((sum,g)=>sum+g.count,0);
 return count(a,b.id)>0&&count(a,b.id)===count(b,a.id);
}
/** Actual reciprocal tissue draw contacts, including the shared male neck.
 * The unpartitioned female compound duct does not imply cystic branch identity. */
export function gallbladderWallRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=graphs.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),r=registry(parts),byId=new Map(parts.map(p=>[p.id,p]));
 if(r)for(const p of r.members){
  const seen=new Set<string>();
  for(const group of p.tissueDrawGroups??[]){
   if(!group.neighbor||seen.has(group.neighbor))continue;
   const target=byId.get(group.neighbor);
   if(!target||target.suppressed||!target.pairedTissue||!target.sectionSolid||target.sourceOffset
    ||target.chunk!==p.chunk||target.positions!==p.positions||target.normals!==p.normals
    ||JSON.stringify(target.geometryOrigin)!==JSON.stringify(p.geometryOrigin))continue;
   const reciprocal=(target.tissueDrawGroups??[]).filter(g=>g.neighbor===p.id).reduce((sum,g)=>sum+g.count,0);
   if(reciprocal!==group.count)continue;
   seen.add(target.id);
   const duct=gallbladderWallInterface(p,target);
   if(!duct&&!r.ids.includes(target.id))continue;
   const kind=duct?'continuous':'adjacent';
   const confluence=duct&&r.sex==='male'&&parts.some(other=>biliaryJunctionInterface(target,other));
   const note=duct?r.sex==='male'?(confluence?'Literal shared gallbladder–cystic neck and open cystic/common hepatic/bile confluence, with reciprocal tissue-only rims. Clinical calibre and upstream/downstream organ interfaces remain provisional.':'Literal shared gallbladder–cystic neck, with an open assembled lumen and reciprocal tissue-only rims. Terminal hepatic/bile junction and clinical calibre remain unverified.'):target.id==='LOCAL:female:cystic_duct:authored'?'Literal shared gallbladder–cystic neck with an open lumen and reciprocal tissue-only rims. The cystic wall joins the retained compound trunk; clinical calibre and finer duct regions remain provisional.':'Shared finite-wall neck interface with the compound native duct. The distinct cystic branch and clinical wall allocation remain unverified.':'Opposed shared interface of adjoining gallbladder tissue intervals. Thickness and local histological proportions remain provisional.';
   for(const [from,to] of [[p,target],[target,p]]){
    const rows=graph.get(from.id)??[];
    if(!rows.some(row=>row.kind===kind&&row.target.id===to.id))rows.push({kind,target:to,note});
    graph.set(from.id,rows);
   }
  }
 }
 if(r){
  const targetConcept:Concept={id:r.source.conceptId,name:'Gallbladder',elements:r.source.replacementParts!.filter(id=>r.ids.includes(id)),ta98Term:'A05.8.02.001',ta98Kind:'exact',coverageLimitation:r.members[0].coverageLimitation};
  const target={...r.members[0],name:'Gallbladder',conceptId:targetConcept.id};
  const reverse={arterial:'supplies',venous:'drains',innervation:'innervates',lymphaticDrainage:'lymphaticTributaries'} as const;
  for(const row of gallbladderRelations(r.reference).get(r.source.id)??[]){
   const kind=reverse[row.kind as keyof typeof reverse];if(!kind)continue;
   const rows=graph.get(row.target.id)??[];
   rows.push({kind,target,targetConcept,note:row.note+' This reference selects the complete registered gallbladder wall; it does not establish a local attachment or supply to every individual coat.',via:row.via});
   graph.set(row.target.id,rows);
  }
 }
 graphs.set(parts,graph);return graph;
}
