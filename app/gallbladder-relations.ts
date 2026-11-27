import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const regions=new Map([
 ['A05.8.02.002','Fundus of gallbladder'],['A05.8.02.004','Body of gallbladder'],
 ['A05.8.02.003','Infundibulum of gallbladder'],['A05.8.02.005','Neck of gallbladder'],
]);
/** These are the stored whole organ and its four source exterior regions.
 * A histological coat, duct or independent donor does not inherit this scope. */
export function nativeGallbladderPart(p:Part):boolean{
 if(p.suppressed||p.system!=='digestive'||!nativePartFrame(p))return false;
 if(p.id==='ZA:Gallbladder')return p.sourceId==='Gallbladder'&&p.name==='Gallbladder';
 if(/^LOCAL:(?:male|female):gallbladder:\d+$/.test(p.id))return p.sourceId==='gallbladder'&&p.name==='Gallbladder';
 const code=p.id.replace('TA98PART:ZA:Gallbladder:','');
 return p.id.startsWith('TA98PART:ZA:Gallbladder:')&&p.sourceId==='Gallbladder'&&regions.get(code)===p.name;
}
const base=(p:Part)=>p.name.replace(/\s*\((?:left|right)\)\s*$/i,'').toLowerCase();
const nativeLeaf=(p:Part)=>!!p.sourceId&&(p.id===`ZA:${p.sourceId}`||!!nativePartFrame(p)?.startsWith('LOCAL:')&&p.id.startsWith(`${nativePartFrame(p)}${p.sourceId}:`));
/** This installed donor has a registered native hepatic junction. Its HRA
 * source name alone does not register other donors or certify perfusion. */
function authoredMaleCysticVessel(parts:Part[],role:'artery'|'vein'):Part|undefined{
 const artery=parts.find(p=>p.id===`TA98CYSTICBLOOD:male:${role}`);
 const parent=parts.find(p=>p.id===`TA98HEPATIC:male:${role==='artery'?'artery':'portal'}`);
 const system=role==='artery'?'arterial':'venous',assembly=`hepatic-connections:male:${role==='artery'?'artery':'portal'}`;
 const model=(p:Part)=>(p as Part&{hepaticSourceModel?:string}).hepaticSourceModel;
 if(!artery||!parent||artery.conceptId!==artery.id||artery.name!==`Cystic ${role}`||artery.system!==system
  ||artery.sourceId!==`HRA:VH_F_cystic_${role}`||artery.ta98Term!==(role==='artery'?'A12.2.12.031':'A12.3.12.012')||model(artery)!=='local-ta98'
  ||parent.system!==system||model(parent)!=='local-ta98'||nativePartFrame(parent)!=='ZA:'
  ||[artery,parent].some(p=>p.suppressed||!p.pairedTissue||!p.sectionSolid||p.normalSigns===undefined||p.sourceOffset
   ||p.tissueWallAssembly!==assembly)
  ||artery.chunk!==parent.chunk||artery.positions!==parent.positions||artery.normals!==parent.normals
  ||JSON.stringify(artery.geometryOrigin)!==JSON.stringify(parent.geometryOrigin)
  ||JSON.stringify(artery.geometryInstance)!==JSON.stringify(parent.geometryInstance)
  ||(artery.geometryScale??1)!==(parent.geometryScale??1))return;
 const rim=(p:Part,id:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id
  &&Number.isInteger(g.start)&&g.start>=0&&g.start%3===0&&Number.isInteger(g.count)&&g.count>0&&g.count%3===0
  &&g.start+g.count<=p.indexCount).reduce((sum,g)=>sum+g.count,0);
 const count=rim(artery,parent.id);
 return count>0&&count===rim(parent,artery.id)?artery:undefined;
}
function registeredMaleGallbladderCovering(parts:Part[]):Part[]{
 const source=parts.find(p=>p.id==='ZA:Gallbladder'),byId=new Map(parts.map(p=>[p.id,p]));
 const ids=['covering','subserosa','muscular','mucosa'].flatMap(layer=>[...regions.keys()].map(code=>`TA98PART:ZA:Gallbladder:GALLWALL:${layer}:${code}`));
 if(!source||source.name!=='Gallbladder'||source.sourceId!=='Gallbladder'||source.sourceOffset
  ||source.replacementParts?.length!==ids.length||!ids.every(id=>source.replacementParts!.includes(id)))return [];
 const members=ids.map(id=>byId.get(id)),first=members[0];
 if(!members.every(p=>p&&!p.suppressed&&p.system==='digestive'&&p.sourceId==='Gallbladder'&&nativePartFrame(p)==='ZA:'
  &&p.bowelWallSource===source.id&&p.tissueWallAssembly==='gallbladder:male'&&p.pairedTissue&&p.sectionSolid&&p.normalSigns!==undefined&&!p.sourceOffset
  &&p.id.startsWith(`TA98PART:ZA:Gallbladder:GALLWALL:${p.tissueWallLayer}:`)
  &&p.chunk===first!.chunk&&p.positions===first!.positions&&p.normals===first!.normals
  &&JSON.stringify(p.geometryOrigin)===JSON.stringify(first!.geometryOrigin)
  &&JSON.stringify(p.geometryInstance)===JSON.stringify(first!.geometryInstance)
  &&(p.geometryScale??1)===(first!.geometryScale??1)))return [];
 return (members as Part[]).filter(p=>p.tissueWallLayer==='covering');
}
const authoredCysticNote=(gall:Part)=>{
 const otherRegion=gall.id.startsWith('TA98PART:')&&!gall.id.endsWith(':A05.8.02.002');
 return (otherRegion?'Whole-gallbladder arterial reference only; no local entry into this exterior region is established. ':'')
  +'The installed cystic artery has a shared finite right-hepatic arterial junction and a verified macroscopic entry into the gallbladder fundus. Microscopic perfusion and exact wall thickness remain provisional; upstream coeliac/hepatic connections outside this local junction remain unverified.';
};
const authoredCysticVenousNote=(gall:Part)=>{
 const otherRegion=gall.id.startsWith('TA98PART:')&&!gall.id.endsWith(':A05.8.02.002');
 return (otherRegion?'Whole-gallbladder venous reference only; no local entry into this exterior region is established. ':'')
  +'The installed cystic vein has a shared finite right-portal venous junction and a verified macroscopic entry into the gallbladder fundus wall without breaching its inner surface. This is one source-specific venous route; microscopic drainage, other venules and exact wall thickness remain provisional.';
};
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Sources: NBK459288 (variable venules, cystic node); PMC4319126
 * (hepatic autonomic route); PMID26698600 (cystic arterial variants).
 * References associate existing upstream anatomy, never certify mesh joins. */
export function gallbladderRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=cache.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>();
 const leaves=parts.filter(p=>!p.suppressed&&nativeLeaf(p));
 const authored=authoredMaleCysticVessel(parts,'artery'),authoredVein=authoredMaleCysticVessel(parts,'vein');
 const add=(from:Part,target:Part,kind:'arterial'|'venous'|'innervation'|'lymphaticDrainage',reverse:'supplies'|'drains'|'innervates'|'lymphaticTributaries',note:string,via?:string[])=>{
  for(const [a,b,k] of [[from,target,kind],[target,from,reverse]] as const){
   const rows=graph.get(a.id)??[];rows.push({kind:k,target:b,note,...(via?.length?{via}: {})});graph.set(a.id,rows);
  }
 };
 for(const gall of parts.filter(nativeGallbladderPart)){
  const frame=nativePartFrame(gall),available=leaves.filter(p=>nativePartFrame(p)===frame);
  const named=(system:Part['system'],names:string[])=>available.filter(p=>p.system===system&&names.includes(base(p)));
  const scope=gall.id.startsWith('TA98PART:')?'Whole-gallbladder reference route; this exterior region has no verified local vascular or nerve attachment. ':'';
  const direct=authored&&frame==='ZA:'?[authored]:named('arterial',['cystic artery','cystic arteries']);
  const right=available.filter(p=>p.system==='arterial'&&['Right hepatic artery','Hepatic artery (right)'].includes(p.name));
  const proper=named('arterial',['proper hepatic artery']);
  const arteries=direct.length?direct:right.length?right:proper;
  for(const artery of arteries){
   if(artery===authored){
    add(gall,artery,'arterial','supplies',authoredCysticNote(gall));
   }else add(gall,artery,'arterial','supplies',scope+'Usual arterial route through a cystic branch, commonly from the right hepatic artery. Origins and branch patterns vary; the native cystic branch and physical lumen junction remain unverified.',direct.length?undefined:right.length?['Cystic artery']:['Right hepatic artery','Cystic artery']);
  }
  const cystic=authoredVein&&frame==='ZA:'?[authoredVein]:named('venous',['cystic vein','cystic veins']);
  if(cystic.length){
   for(const vein of cystic)add(gall,vein,'venous','drains',vein===authoredVein?authoredCysticVenousNote(gall):scope+'This source depicts a cystic venous route. Small venules and direct liver-bed drainage also occur; there is no assumed single universal cystic vein or verified source junction.');
  }else{
   for(const vein of named('venous',['hepatic portal vein','portal vein']))add(gall,vein,'venous','drains',scope+'Portal-tree association through variable cystic venous tributaries, particularly around the neck and cystic duct. This does not make the extrahepatic portal trunk the direct outlet of every region; local tributaries remain unverified.',['Cystic venous tributaries']);
   for(const vein of named('venous',['hepatic veins']))add(gall,vein,'venous','drains',scope+'Liver-bed venules can reach hepatic venous outflow through the liver sinusoids. These are downstream reference structures; short venules, sinusoids and local source junctions remain unverified.',['Short cystic venules','Hepatic sinusoids']);
  }
  for(const nerve of named('nervous',['vagus nerve','vagus nerve (x)']))add(gall,nerve,'innervation','innervates',scope+'Parasympathetic contribution through the hepatic plexus. This upstream vagal source does not represent the individual gallbladder nerve branches, sympathetic contribution or intramural fibers; local innervation remains unverified.',['Hepatic plexus']);
  for(const node of named('lymphatic',['cystic node','cystic nodes','cystic lymph node']))add(gall,node,'lymphaticDrainage','lymphaticTributaries',scope+'Cystic-node lymphatic association. Lymphatic channels, regional drainage variants and native tissue attachments remain unverified.');
 }
 if(authored||authoredVein)for(const covering of registeredMaleGallbladderCovering(parts)){
  if(authored)add(covering,authored,'arterial','supplies',authoredCysticNote(covering));
  if(authoredVein)add(covering,authoredVein,'venous','drains',authoredCysticVenousNote(covering));
 }
 cache.set(parts,graph);return graph;
}

/** Bypass unrelated native-bowel frame filters only for a pair already
 * registered by this source-specific blood-vessel graph, in either direction. */
export function isReviewedCysticBloodRelation(kind:ResolvedRelation['kind'],source:Part,target:Part,parts:Part[]):boolean{
 if(kind!=='arterial'&&kind!=='supplies'&&kind!=='venous'&&kind!=='drains')return false;
 const direct=kind==='arterial'||kind==='venous',role=kind==='arterial'||kind==='supplies'?'artery':'vein';
 const artery=direct?target:source;
 if(artery.id!==`TA98CYSTICBLOOD:male:${role}`||!parts.includes(artery)||direct&&!parts.includes(source))return false;
 return !!gallbladderRelations(parts).get(source.id)?.some(row=>row.kind===kind&&row.target.id===target.id);
}
