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
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Sources: NBK459288 (variable venules, cystic node); PMC4319126
 * (hepatic autonomic route); PMID26698600 (cystic arterial variants).
 * References associate existing upstream anatomy, never certify mesh joins. */
export function gallbladderRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=cache.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>();
 const leaves=parts.filter(p=>!p.suppressed&&nativeLeaf(p));
 const add=(from:Part,target:Part,kind:'arterial'|'venous'|'innervation'|'lymphaticDrainage',reverse:'supplies'|'drains'|'innervates'|'lymphaticTributaries',note:string,via?:string[])=>{
  for(const [a,b,k] of [[from,target,kind],[target,from,reverse]] as const){
   const rows=graph.get(a.id)??[];rows.push({kind:k,target:b,note,...(via?.length?{via}: {})});graph.set(a.id,rows);
  }
 };
 for(const gall of parts.filter(nativeGallbladderPart)){
  const frame=nativePartFrame(gall),available=leaves.filter(p=>nativePartFrame(p)===frame);
  const named=(system:Part['system'],names:string[])=>available.filter(p=>p.system===system&&names.includes(base(p)));
  const scope=gall.id.startsWith('TA98PART:')?'Whole-gallbladder reference route; this exterior region has no verified local vascular or nerve attachment. ':'';
  const direct=named('arterial',['cystic artery','cystic arteries']);
  const right=available.filter(p=>p.system==='arterial'&&['Right hepatic artery','Hepatic artery (right)'].includes(p.name));
  const proper=named('arterial',['proper hepatic artery']);
  const arteries=direct.length?direct:right.length?right:proper;
  for(const artery of arteries)add(gall,artery,'arterial','supplies',scope+'Usual arterial route through a cystic branch, commonly from the right hepatic artery. Origins and branch patterns vary; the native cystic branch and physical lumen junction remain unverified.',direct.length?undefined:right.length?['Cystic artery']:['Right hepatic artery','Cystic artery']);
  const cystic=named('venous',['cystic vein','cystic veins']);
  if(cystic.length){
   for(const vein of cystic)add(gall,vein,'venous','drains',scope+'This source depicts a cystic venous route. Small venules and direct liver-bed drainage also occur; there is no assumed single universal cystic vein or verified source junction.');
  }else{
   for(const vein of named('venous',['hepatic portal vein','portal vein']))add(gall,vein,'venous','drains',scope+'Portal-tree association through variable cystic venous tributaries, particularly around the neck and cystic duct. This does not make the extrahepatic portal trunk the direct outlet of every region; local tributaries remain unverified.',['Cystic venous tributaries']);
   for(const vein of named('venous',['hepatic veins']))add(gall,vein,'venous','drains',scope+'Liver-bed venules can reach hepatic venous outflow through the liver sinusoids. These are downstream reference structures; short venules, sinusoids and local source junctions remain unverified.',['Short cystic venules','Hepatic sinusoids']);
  }
  for(const nerve of named('nervous',['vagus nerve','vagus nerve (x)']))add(gall,nerve,'innervation','innervates',scope+'Parasympathetic contribution through the hepatic plexus. This upstream vagal source does not represent the individual gallbladder nerve branches, sympathetic contribution or intramural fibers; local innervation remains unverified.',['Hepatic plexus']);
  for(const node of named('lymphatic',['cystic node','cystic nodes','cystic lymph node']))add(gall,node,'lymphaticDrainage','lymphaticTributaries',scope+'Cystic-node lymphatic association. Lymphatic channels, regional drainage variants and native tissue attachments remain unverified.');
 }
 cache.set(parts,graph);return graph;
}
