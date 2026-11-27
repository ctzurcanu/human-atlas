import {catalogueSourceReference} from "./source-reference";
import type {Concept,Part} from './anatomy';
import type {RelationKind,ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const names=new Set(['Coeliac trunk','Celiac trunk','Left gastric artery','Gastric artery (left)','Common hepatic artery','Proper hepatic artery','Gastroduodenal artery','Splenic artery']);
const nativeLeaf=(p:Part)=>!!p.sourceId&&(p.id===`ZA:${p.sourceId}`||!!nativePartFrame(p)?.startsWith('LOCAL:')&&p.id.startsWith(`${nativePartFrame(p)}${p.sourceId}:`));
const splenicIds={trunk:'TA98PART:ZA:Splenic artery:trunk',branches:'TA98PART:ZA:Splenic artery:A12.2.12.051'} as const;
type SplenicRole=keyof typeof splenicIds;
type SplenicOwnership={schema?:string;role?:string};
const sourceUrl='/ATTRIBUTION.md#open-anatomy-source';
const same=(a:unknown,b:unknown)=>JSON.stringify(a)===JSON.stringify(b);
const nativePose=(p:Part)=>!p.suppressed&&!p.sourceOffset&&!p.geometryInstance&&(p.geometryScale??1)===1&&catalogueSourceReference(p.provenance?.url)===sourceUrl&&p.license==='CC BY-SA 4.0';
const splenicCandidate=(p:Part)=>p.id.startsWith('TA98PART:ZA:Splenic artery:')||(p as Part&{splenicBranchOwnership?:SplenicOwnership}).splenicBranchOwnership?.schema==='source-bound-hilar-tree/v1';
function ownedSplenicRole(p:Part):SplenicRole|undefined{
 const role=(Object.keys(splenicIds) as SplenicRole[]).find(role=>p.id===splenicIds[role]),ownership=(p as Part&{splenicBranchOwnership?:SplenicOwnership}).splenicBranchOwnership;
 if(!role||!nativePose(p)||nativePartFrame(p)!=='ZA:'||p.system!=='arterial'||p.sourceId!=='Splenic artery'||p.conceptId!==p.id
  ||p.pairedTissue!==true||p.sectionSolid!==true||p.normalSigns===undefined||p.tissueWallAssembly!=='splenic-artery:male:owned-branches'
  ||p.ta98Term!==(role==='trunk'?'A12.2.12.040':'A12.2.12.051')||!same(p.geometryOrigin,[0,1.25,0])
  ||ownership?.schema!=='source-bound-hilar-tree/v1'||ownership.role!==role)return;
 return role;
}
function sharedSplenicRim(trunk:Part,branches:Part):boolean{
 if(ownedSplenicRole(trunk)!=='trunk'||ownedSplenicRole(branches)!=='branches'||trunk.chunk!==branches.chunk||!same(trunk.geometryOrigin,branches.geometryOrigin))return false;
 const rim=(a:Part,b:Part)=>(a.tissueDrawGroups??[]).filter(g=>g.neighbor===b.id&&Number.isInteger(g.start)&&g.start>=0&&g.start%3===0&&Number.isInteger(g.count)&&g.count>0&&g.count%3===0&&g.start+g.count<=a.indexCount);
 const a=rim(trunk,branches),b=rim(branches,trunk);
 return a.length===1&&b.length===1&&a[0].count===b[0].count;
}
const spleenRoles=['serosa','capsule','pulp','trabeculae'] as const;
const spleenTerms=['A13.2.01.018','A13.2.01.002','A13.2.01.004','A13.2.01.003'];
function nativeMaleSpleenTissue(p:Part):boolean{
 const role=spleenRoles.findIndex(role=>p.id===`TA98PART:male:Spleen:${role}`);
 return role>=0&&nativePose(p)&&p.system==='lymphatic'&&p.sourceId===`Spleen:finite:${spleenRoles[role]}`
  &&(p as Part&{spleenSourcePartId?:string}).spleenSourcePartId==='ZA:Spleen'&&p.ta98Term===spleenTerms[role]
  &&p.tissueWallAssembly==='splenic:male'&&p.pairedTissue===true&&p.sectionSolid===true&&p.normalSigns!==undefined&&same(p.geometryOrigin,[0,1,0]);
}
function ownedSplenicFamily(parts:Part[]):{trunk:Part;branches:Part}|undefined{
 const trunk=parts.find(p=>p.id===splenicIds.trunk),branches=parts.find(p=>p.id===splenicIds.branches),retired=parts.find(p=>p.id==='ZA:Splenic artery');
 if(!trunk||!branches||!sharedSplenicRim(trunk,branches)||retired&&(!retired.suppressed||retired.replacementParts?.length!==2||!Object.values(splenicIds).every(id=>retired.replacementParts?.includes(id))))return;
 return {trunk,branches};
}
function registeredMaleSpleen(parts:Part[]):Part[]|undefined{
 const tissues=spleenRoles.map(role=>parts.find(p=>p.id===`TA98PART:male:Spleen:${role}`));
 if(tissues.some(p=>!p||!nativeMaleSpleenTissue(p)))return;
 return tissues as Part[];
}
const nativeCoeliac=(p:Part)=>p.id==='ZA:Coeliac trunk'&&p.conceptId===p.id&&p.sourceId==='Coeliac trunk'&&p.system==='arterial'&&nativePose(p)&&nativePartFrame(p)==='ZA:';
function ownedSplenicPairCompatible(kind:RelationKind,a:Part,b:Part):boolean{
 if(a.suppressed||b.suppressed)return false;
 const ar=ownedSplenicRole(a),br=ownedSplenicRole(b);
 if(ar==='trunk'&&br==='branches')return kind==='after'&&sharedSplenicRim(a,b);
 if(ar==='branches'&&br==='trunk')return kind==='before'&&sharedSplenicRim(b,a);
 if(ar==='trunk'&&nativeCoeliac(b))return kind==='before';
 if(br==='trunk'&&nativeCoeliac(a))return kind==='after';
 const carrier=(p:Part)=>p.id==='TA98PART:male:Spleen:serosa'&&nativeMaleSpleenTissue(p);
 return ar==='branches'&&kind==='supplies'&&carrier(b)||br==='branches'&&kind==='arterial'&&carrier(a);
}
const splenicSupplyNote='Whole-spleen arterial reference from the source-bound distal hilar tree. Recipient tissue entry, biological perfusion and microscopic distribution remain unverified; no individual splenic tissue layer is assigned direct supply.';
/** Complete named recipient only: no arbitrary selection or donor aggregate. */
export function foregutWholeFunctionalRelations(concept:Concept,selected:Part[],parts:Part[]):ResolvedRelation[]|undefined{
 if(concept.id!=='ZA:Spleen')return;
 const tissues=registeredMaleSpleen(parts),family=ownedSplenicFamily(parts),ids=new Set(selected.map(p=>p.id));
 if(!family)return;
 if(!tissues||selected.length!==4||ids.size!==4||concept.elements.length!==4||!tissues.every(p=>selected.includes(p)&&ids.has(p.id)&&concept.elements.includes(p.id)))return [];
 return [{kind:'arterial',target:family.branches,note:splenicSupplyNote}];
}

/** A native arterial itinerary cannot jump into an independently placed donor.
 * This guard also applies to reciprocal generic routes in the final resolver. */
export function foregutArterialFrameCompatible(kind:RelationKind,a:Part,b:Part):boolean{
 if(['before','after','adjacent','arterial','supplies'].includes(kind)&&[a,b].some(splenicCandidate))return ownedSplenicPairCompatible(kind,a,b);
 if(!['before','after','adjacent'].includes(kind)||![a,b].some(p=>names.has(p.name)))return true;
 const frame=nativePartFrame(a),other=nativePartFrame(b);
 if(!frame&&!other)return true;
 return !a.suppressed&&!b.suppressed&&a.system==='arterial'&&b.system==='arterial'&&!!frame&&frame===other&&[a,b].every(p=>!names.has(p.name)||nativeLeaf(p));
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** NBK459241: usual foregut arterial branching, with anatomical variants.
 * Named sources and shared rims do not prove complete walls or lumen joins. */
export function foregutArterialRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),active=parts.filter(p=>!p.suppressed&&p.system==='arterial'&&nativePartFrame(p)&&nativeLeaf(p)&&names.has(p.name));
 for(const [parents,children] of [
  [['Coeliac trunk','Celiac trunk'],['Left gastric artery','Gastric artery (left)','Common hepatic artery','Splenic artery']],
  [['Common hepatic artery'],['Proper hepatic artery','Gastroduodenal artery']],
 ])for(const a of active.filter(p=>parents.includes(p.name)))for(const b of active.filter(p=>children.includes(p.name)&&nativePartFrame(p)===nativePartFrame(a))){
  const note='Usual named foregut arterial branching in this native model. Arterial variants occur; complete vessel walls, biological lumen continuity and source attachment anatomy remain unverified.';
  for(const [source,target,kind] of [[a,b,'after'],[b,a,'before']] as const){const rows=graph.get(source.id)??[];rows.push({kind,target,note});graph.set(source.id,rows);}
 }
 const family=ownedSplenicFamily(parts);
 if(family){
  const add=(source:Part,target:Part,kind:RelationKind,note:string,targetConcept?:Concept)=>{const rows=graph.get(source.id)??[];rows.push({kind,target,note,...(targetConcept?{targetConcept}:{})});graph.set(source.id,rows);};
  const interfaceNote='Source-bound owned splenic cells share the same reciprocal annular tissue interface. The source-station ownership boundary is provisional; recipient perfusion and other vascular junctions remain unverified.';
  add(family.trunk,family.branches,'after',interfaceNote);add(family.branches,family.trunk,'before',interfaceNote);
  const coeliac=parts.find(nativeCoeliac);
  if(coeliac){const note='Named native coeliac parent of the splenic arterial trunk. The coeliac-to-splenic communicating mesh junction remains unverified.';add(coeliac,family.trunk,'after',note);add(family.trunk,coeliac,'before',note);}
  const tissues=registeredMaleSpleen(parts);
  if(tissues){const whole:Concept={id:'ZA:Spleen',name:'Spleen',elements:tissues.map(p=>p.id)};add(family.branches,tissues[0],'supplies',splenicSupplyNote,whole);add(tissues[0],family.branches,'arterial',splenicSupplyNote);}
 }
 cache.set(parts,graph);return graph;
}
