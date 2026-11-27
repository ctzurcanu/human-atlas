import type {Atlas,Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';

/** This marker describes stored tissues, not a claim of anatomical acceptance. */
export interface ThyroidWallAssembly {
 schema:'ta98-thyroid-wall-v1';
 sex:'male'|'female';
 wholeConcept:string;
 regionConcepts:string[];
 capsuleConcept:string;
 parenchymaConcept:string;
 stromaConcept?:string;
 lobulesConcept?:string;
 capsuleParts:string[];
 parenchymaParts:string[];
 stromaParts?:string[];
 /** Alias union of genuinely modeled parenchymal domains; no extra coat. */
 lobuleParts?:string[];
 /** Local extraction territories; these do not establish microscopic lobules. */
 illustrativeDomainConcept?:string;
 illustrativeDomainParts?:string[];
 /** Shared finite vessel walls are interface partners, not endocrine coats. */
 vascularParts?:string[];
 sourceReview:{status:'provisional'|'reviewed'};
}
export interface ThyroidWallRegistry {
 marker:ThyroidWallAssembly;
 whole:Concept;
 regions:Concept[];
 capsule:Concept;
 parenchyma:Concept;
 stroma?:Concept;
 lobules?:Concept;
 domains?:Concept;
 members:Part[];
}
const registries=new WeakMap<Atlas,ThyroidWallRegistry|null>();
const sameMembers=(a:readonly string[],b:readonly string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));
const ids=(value:unknown):value is string[]=>Array.isArray(value)&&value.length>0&&value.every(id=>typeof id==='string'&&id.length>0)&&new Set(value).size===value.length;

/** Only a registered shared-bank finite assembly activates tissue controls.
 * Provisional boundaries remain provisional; metadata cannot validate shape.
 */
export function thyroidWallRegistry(atlas:Atlas):ThyroidWallRegistry|null {
 if(registries.has(atlas))return registries.get(atlas)!;
 const marker=atlas.thyroidWallAssembly;
 const reject=()=>{registries.set(atlas,null);return null;};
 if(!marker||marker.schema!=='ta98-thyroid-wall-v1'||!['male','female'].includes(marker.sex)||marker.sex!==atlas.sex
  ||!['provisional','reviewed'].includes(marker.sourceReview?.status)
  ||!ids(marker.capsuleParts)||!ids(marker.parenchymaParts)||!ids(marker.regionConcepts))return reject();
 const roles:[string,readonly string[]][]=[['capsule',marker.capsuleParts],['parenchyma',marker.parenchymaParts]];
 if(marker.stromaParts?.length){
  if(!ids(marker.stromaParts)||!marker.stromaConcept)return reject();
  roles.push(['stroma',marker.stromaParts]);
 }else if(marker.stromaConcept)return reject();
 if(marker.lobuleParts?.length){
  if(!ids(marker.lobuleParts)||!marker.lobulesConcept||marker.lobuleParts.some(id=>!marker.parenchymaParts.includes(id)))return reject();
 }else if(marker.lobulesConcept)return reject();
 if(marker.illustrativeDomainParts?.length){
  if(!ids(marker.illustrativeDomainParts)||!marker.illustrativeDomainConcept
   ||marker.illustrativeDomainParts.some(id=>!marker.parenchymaParts.includes(id)))return reject();
 }else if(marker.illustrativeDomainConcept)return reject();
 const allIds=roles.flatMap(([,members])=>members);
 if(new Set(allIds).size!==allIds.length)return reject();
 if(marker.vascularParts&&(!ids(marker.vascularParts)||marker.vascularParts.some(id=>allIds.includes(id))))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const members=allIds.map(id=>byId.get(id)),first=members[0];
 if(!first||members.some(p=>!p))return reject();
 const whole=concepts.get(marker.wholeConcept),capsule=concepts.get(marker.capsuleConcept),parenchyma=concepts.get(marker.parenchymaConcept);
 const tissueConcept=(c:Concept|undefined,term:string,parts:readonly string[])=>!!c&&c.ta98Term===term&&c.ta98Kind==='exact'&&c.partOf===marker.wholeConcept&&sameMembers(c.elements,parts);
 if(!whole||whole.ta98Term!=='A11.3.00.001'||whole.ta98Kind!=='exact'||!sameMembers(whole.elements,allIds)
  ||!tissueConcept(capsule,'A11.3.00.006',marker.capsuleParts)||!tissueConcept(parenchyma,'A11.3.00.008',marker.parenchymaParts)
  ||marker.stromaParts?.length&&!tissueConcept(concepts.get(marker.stromaConcept!),'A11.3.00.007',marker.stromaParts)
  ||marker.lobuleParts?.length&&!tissueConcept(concepts.get(marker.lobulesConcept!),'A11.3.00.009',marker.lobuleParts))return reject();
 const domains=marker.illustrativeDomainConcept?concepts.get(marker.illustrativeDomainConcept):undefined;
 if(marker.illustrativeDomainParts?.length&&(!domains||domains.ta98Term||domains.partOf!==marker.parenchymaConcept
  ||!sameMembers(domains.elements,marker.illustrativeDomainParts)))return reject();
 const roleById=new Map(roles.flatMap(([role,members])=>members.map(id=>[id,role] as const)));
 const assembly=`thyroid:${marker.sex}`;
 const vascularIds=new Set(marker.vascularParts??[]);
 const vascular=[...vascularIds].map(id=>byId.get(id));
 if(members.some(p=>!p||p.suppressed||p.sourceOffset||p.system!=='endocrine'||!p.pairedTissue||!p.sectionSolid
  ||p.tissueWallAssembly!==assembly||p.tissueWallLayer!==roleById.get(p.id)||p.normalSigns===undefined||!p.tissueAppearance
  ||p.geometryOrigin?.length!==3||!p.geometryOrigin.every(Number.isFinite)
  ||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount
  ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin))
  ||atlas.parts.some(p=>!p.suppressed&&p.tissueWallAssembly===assembly&&!roleById.has(p.id)))return reject();
 if(vascular.some(p=>!p||p.suppressed||p.sourceOffset||!['arterial','venous'].includes(p.system)
  ||!p.pairedTissue||!p.sectionSolid||p.normalSigns===undefined
  ||p.tissueWallAssembly!==`thyroid-vascular:${marker.sex}`||p.tissueWallLayer!=='vascular-wall'
  ||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount
  ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return reject();
 const contactCount=(p:Part,id:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id).reduce((n,g)=>n+g.count,0);
 for(const p of [...members,...vascular] as Part[]){
  const groups=p.tissueDrawGroups;if(!groups?.length)return reject();
  let end=0;
  for(const group of groups){
   if(!Number.isInteger(group.start)||!Number.isInteger(group.count)||group.start!==end||group.count<=0||group.count%3!==0)return reject();
   end+=group.count;
   if(group.neighbor){const peer=byId.get(group.neighbor);if(!peer||peer.id===p.id
    ||!roleById.has(peer.id)&&!vascularIds.has(peer.id)||contactCount(p,peer.id)!==contactCount(peer,p.id))return reject();}
  }
  if(end!==p.indexCount)return reject();
  // Each capsule interval meets actual underlying gland tissue. Merely
  // naming two independent exteriors must not establish peeling layers.
  if(p.tissueWallLayer==='capsule'&&!groups.some(g=>g.neighbor&&['parenchyma','stroma'].includes(roleById.get(g.neighbor)??'')))return reject();
 }
 // Every registered vascular component must belong to the contact graph
 // reached from actual gland tissue; an unrelated closed vessel is ineligible.
 const reached=new Set(allIds),pending=[...allIds];
 for(let i=0;i<pending.length;i++)for(const g of byId.get(pending[i])!.tissueDrawGroups??[]){
  if(g.neighbor&&vascularIds.has(g.neighbor)&&!reached.has(g.neighbor)){reached.add(g.neighbor);pending.push(g.neighbor);}
 }
 if([...vascularIds].some(id=>!reached.has(id)))return reject();
 const regions=marker.regionConcepts.map(id=>concepts.get(id));
 if(regions.some(c=>!c||c.partOf!==whole.id||!['A11.3.00.002','A11.3.00.003','A11.3.00.004'].includes(c.ta98Term??'')
  ||!ids(c.elements)||c.elements.some(id=>!roleById.has(id))
  ||!c.elements.some(id=>roleById.get(id)==='capsule')||!c.elements.some(id=>roleById.get(id)==='parenchyma'))
  ||allIds.some(id=>regions.filter(c=>c?.elements.includes(id)).length!==1))return reject();
 const result={marker,whole,regions:regions as Concept[],capsule:capsule!,parenchyma:parenchyma!,
  ...(marker.stromaConcept?{stroma:concepts.get(marker.stromaConcept)!}:{}),
  ...(marker.lobulesConcept?{lobules:concepts.get(marker.lobulesConcept)!}:{}),
  ...(domains?{domains}:{}),members:members as Part[]};
 registries.set(atlas,result);return result;
}

/** Anatomical regions are membership unions, not concentric peel stages. */
export function thyroidWallPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null {
 const r=thyroidWallRegistry(atlas);if(!r||choice.id==='selection-set'||!choice.elements.length)return null;
 const complete=sameMembers(choice.elements,r.whole.elements);
 const region=r.regions.find(c=>sameMembers(choice.elements,c.elements));
 const term=choice.ta98Term??choice.terminology?.ta98??choice.id.replace(/^hierarchy:guest:ta98:/,'');
 const lobes=r.regions.filter(c=>c.ta98Term==='A11.3.00.002');
 const collectiveLobes=term==='A11.3.00.002'&&sameMembers(choice.elements,lobes.flatMap(c=>c.elements));
 const capsule=sameMembers(choice.elements,r.capsule.elements),parenchyma=sameMembers(choice.elements,r.parenchyma.elements);
 const internalUnion=[r.stroma,r.lobules,r.domains].some(c=>c&&sameMembers(choice.elements,c.elements));
 if(!complete&&!region&&!collectiveLobes&&!capsule&&!parenchyma&&!internalUnion)return null;
 const selected=new Set(choice.elements),inner=r.members.filter(p=>p.tissueWallLayer!=='capsule'&&selected.has(p.id));
 const steps=[{name:'Fibrous capsule',ids:r.marker.capsuleParts.filter(id=>selected.has(id))},
  {name:inner.some(p=>p.tissueWallLayer!=='parenchyma')?'Internal thyroid tissue':'Parenchyma',ids:inner.map(p=>p.id)}].filter(step=>step.ids.length);
 return steps.flatMap(s=>s.ids).length===selected.size?steps:null;
}
