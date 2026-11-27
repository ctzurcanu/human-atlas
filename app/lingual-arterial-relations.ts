import {catalogueSourceReference} from "./source-reference";
import type {Concept,Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';

type Side='left'|'right';
const sourceUrl='/ATTRIBUTION.md#open-anatomy-source';
const roles={lingual:'A12.2.05.015',suprahyoid:'A12.2.05.016','dorsal-lingual-0':'A12.2.05.017','dorsal-lingual-1':'A12.2.05.017',sublingual:'A12.2.05.018','deep-lingual':'A12.2.05.019'} as const;
type Role=keyof typeof roles;
const arteryId=(role:Role,side:Side)=>`TA98ARTERIAL:male:${role}:${side}:finite-wall`;
const unchangedPose=(p:Part)=>!p.suppressed&&!p.sourceOffset&&!p.geometryInstance&&(p.geometryScale??1)===1&&catalogueSourceReference(p.provenance?.url)===sourceUrl;
const same=(a:unknown,b:unknown)=>JSON.stringify(a)===JSON.stringify(b);
const origins={left:[.016594554721975487,1.623810152874436,.04152941869862918],right:[-.01985116613921688,1.6236474710160156,.043219807199565606]};
type LingualRevisionPart=Part&{lingualArterialWallRevision?:string;carotidFiniteWallRevision?:string;carotidFiniteWallBankSha256?:string};
const leftCommonRevision='left-four-cell-common-bank-v2';
const leftCommonBankSha256='e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36';
const rightCommonRevision='right-four-cell-common-bank-v1';
const rightCommonBankSha256='7b7bcbf93627134a3e240bb10716c202a1455e4779a4f429ae309a66145eb737';
function lingualBankFrame(p:Part,side:Side):'original'|'left-common'|'right-common'|undefined{
 const marked=p as LingualRevisionPart;
 if(marked.carotidFiniteWallRevision!==undefined||marked.carotidFiniteWallBankSha256!==undefined){
  const revision=side==='left'?leftCommonRevision:rightCommonRevision,bank=side==='left'?leftCommonBankSha256:rightCommonBankSha256;
  return marked.carotidFiniteWallRevision===revision&&marked.carotidFiniteWallBankSha256===bank
   &&same(p.geometryOrigin,[0,1.625,0])?(side==='left'?'left-common':'right-common'):undefined;
 }
 return same(p.geometryOrigin,origins[side])?'original':undefined;
}

/** A constructed region belongs only to its side's audited common bank. */
export function isInstalledCarotidBifurcationPart(p:Part):boolean{
 const match=/^TA98ARTERIAL:male:carotid-bifurcation:(left|right):finite-wall$/.exec(p.id);
 if(!match)return false;
 const side=match[1] as Side;
 return p.conceptId===p.id&&p.ta98Term==='A12.2.04.009'
  &&p.sourceId===`Constructed ${side} carotid bifurcation tissue region`
  &&unchangedPose(p)&&p.license==='CC BY-SA 4.0'&&p.system==='arterial'
  &&p.pairedTissue===true&&p.sectionSolid===true&&Number.isInteger(p.normalSigns)&&p.normalSigns!>=0&&p.uvs!==undefined
  &&p.sectionAssembly===`${side}-carotid-lingual-coupled-finite-wall:male`&&lingualBankFrame(p,side)===`${side}-common`;
}

/** Frame identity for these installed cells only; the complete-family query
 * below additionally requires the actual shared bank and reciprocal rims. */
export function isInstalledLingualArterialPart(p:Part):boolean{
 const match=/^TA98ARTERIAL:male:(lingual|suprahyoid|dorsal-lingual-[01]|sublingual|deep-lingual):(left|right):finite-wall$/.exec(p.id);
 if(!match)return false;
 const role=match[1] as Role,side=match[2] as Side;
 return unchangedPose(p)&&p.system==='arterial'&&p.pairedTissue===true&&p.sectionSolid===true&&p.normalSigns!==undefined
  &&p.sourceId==='independent-lingual-course/native-ECA-and-recipient-constrained'&&p.ta98Term===roles[role]
  &&p.conceptId===(role.startsWith('dorsal-lingual-')?`TA98ARTERIAL:male:dorsal-lingual:${side}:whole`:p.id)
  &&p.sectionAssembly===`${side}-lingual-arterial-family:male`&&lingualBankFrame(p,side)!==undefined
  &&(p as LingualRevisionPart).lingualArterialWallRevision===`male-${side}-lingual-finite-network-v${side==='left'?5:2}`;
}

/** Only these reviewed bank layouts establish the seven-cell lingual family.
 * A matching vessel name, TA98 code or donor ID does not register a junction. */
export function registeredLingualFamily(parts:Part[],side:Side):Map<Role,Part>|undefined{
 const byId=new Map(parts.map(p=>[p.id,p])),parent=byId.get(`ZA:External carotid artery.${side[0]}`);
 const branches=Object.keys(roles).map(role=>byId.get(arteryId(role as Role,side)));
 if(!parent||branches.some(p=>!p))return;
 const members=[parent,...branches as Part[]],revision=`male-${side}-lingual-finite-network-v${side==='left'?5:2}`;
 const frame=lingualBankFrame(parent,side);if(frame===undefined)return;
 if(parent.sourceId!==`External carotid artery.${side[0]}`||parent.conceptId!==parent.id
  ||members.some(p=>!unchangedPose(p)||p.system!=='arterial'||!p.pairedTissue||!p.sectionSolid||p.normalSigns===undefined
   ||p.sectionAssembly!==`${side}-lingual-arterial-family:male`
   ||(p as LingualRevisionPart).lingualArterialWallRevision!==revision
   ||lingualBankFrame(p,side)!==frame||p.chunk!==parent.chunk||p.positions!==parent.positions||p.normals!==parent.normals||p.uvs!==parent.uvs||p.vertexCount!==parent.vertexCount))return;
 const family=new Map<Role,Part>();
 for(const role of Object.keys(roles) as Role[]){
  const p=byId.get(arteryId(role,side))!;
  if(!isInstalledLingualArterialPart(p))return;
  family.set(role,p);
 }
 const rim=(p:Part,id:string)=>{
  const groups=(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id);
  return groups.length===1&&groups[0].count===72&&Number.isInteger(groups[0].start)&&groups[0].start>=0
   &&groups[0].start%3===0&&groups[0].start+72<=p.indexCount;
 };
 const main=family.get('lingual')!;
 for(const child of members.filter(p=>p!==main))if(!rim(main,child.id)||!rim(child,main.id))return;
 return family;
}

const tongueIds=['mucosa','aponeurosis','intrinsic:superior-longitudinal','intrinsic:inferior-longitudinal-left','intrinsic:inferior-longitudinal-right','intrinsic:transverse-vertical-shared','styloglossus:left','styloglossus:right'].map(role=>`TA98TONGUE:male:${role}`)
 .concat(['Genioglossus','Hyoglossus'].flatMap(name=>['l','r'].map(side=>`ZA:${name} muscle.${side}`)));
function registeredTongue(parts:Part[]):Part[]|undefined{
 const byId=new Map(parts.map(p=>[p.id,p])),members=tongueIds.map(id=>byId.get(id));
 if(members.some(p=>!p||!unchangedPose(p)))return;
 for(const p of members as Part[]){
  const styloid=p.id.startsWith('TA98TONGUE:male:styloglossus:');
  if(!same(p.geometryOrigin,[0,styloid?1.5:1.625,0]))return;
  if(p.id.startsWith('TA98TONGUE:')&&!p.sectionSolid)return;
 }
 return members as Part[];
}
function nativeRecipient(parts:Part[],name:string,side:Side):Part|undefined{
 const id=`ZA:${name}.${side[0]}`,p=parts.find(p=>p.id===id);
 return p&&p.sourceId===`${name}.${side[0]}`&&p.conceptId===id&&unchangedPose(p)
  &&p.system===(name==='Hyoglossus muscle'?'muscular':'digestive')
  &&same(p.geometryOrigin,name==='Hyoglossus muscle'?[0,1.625,0]:undefined)?p:undefined;
}
const wholeNote='Whole-tongue arterial reference. The installed branch family communicates with the external carotid artery, and its dorsal/deep terminal annuli enter the tongue parenchymal envelope. This does not establish entry into an individual mucosal or intrinsic-muscle layer; microscopic perfusion and illustrative wall thickness remain provisional.';
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function lingualArterialRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=cache.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),tongue=registeredTongue(parts);
 const add=(p:Part,row:ResolvedRelation)=>{const rows=graph.get(p.id)??[];rows.push(row);graph.set(p.id,rows);};
 for(const side of ['left','right'] as const){
  const family=registeredLingualFamily(parts,side);if(!family)continue;
  const main=family.get('lingual')!,parent=parts.find(p=>p.id===`ZA:External carotid artery.${side[0]}`)!;
  for(const [from,to] of [[parent,main],...[...family.values()].filter(p=>p!==main).map(p=>[main,p])] as [Part,Part][]){
   const note='Registered reciprocal annular interface in the installed finite lingual arterial family. Other native carotid branch mouths and microscopic perfusion remain unfinished.';
   add(from,{kind:'before',target:to,note});add(to,{kind:'after',target:from,note});
  }
  for(const [role,name] of [['suprahyoid','Hyoglossus muscle'],['sublingual','Sublingual gland']] as const){
   const recipient=nativeRecipient(parts,name,side);if(!recipient)continue;
   const artery=family.get(role)!,note=`The installed ${role} terminal annulus enters the ipsilateral ${name.toLowerCase()} tissue. This is a macroscopic recipient interface; microscopic perfusion and exact wall thickness remain provisional.`;
   add(recipient,{kind:'arterial',target:artery,note});add(artery,{kind:'supplies',target:recipient,note});
  }
  if(tongue){
   const concept:Concept={id:'TA98TONGUE:male:whole',name:'Tongue',elements:tongue.map(p=>p.id)};
   // Mucosa is only the selection carrier for the whole assembly. No direct
   // mucosal-coat supply row is introduced.
   for(const role of ['lingual','dorsal-lingual-0','dorsal-lingual-1','deep-lingual'] as const)
    add(family.get(role)!,{kind:'supplies',target:tongue[0],targetConcept:concept,note:wholeNote});
  }
 }
 cache.set(parts,graph);return graph;
}
export function lingualWholeFunctionalRelations(concept:Concept,selected:Part[],parts:Part[]):ResolvedRelation[]|undefined{
 if(concept.id!=='TA98TONGUE:male:whole')return;
 const tongue=registeredTongue(parts),ids=new Set(selected.map(p=>p.id));
 if(!tongue||ids.size!==tongueIds.length||selected.length!==ids.size||concept.elements.length!==ids.size
  ||!tongueIds.every(id=>ids.has(id)&&concept.elements.includes(id))||selected.some(p=>!parts.includes(p)))return [];
 const graph=lingualArterialRelations(parts),result:ResolvedRelation[]=[];
 for(const p of parts)for(const row of graph.get(p.id)??[])
  if(row.kind==='supplies'&&row.targetConcept?.id===concept.id)result.push({kind:'arterial',target:p,note:row.note});
 return result;
}
export function isReviewedLingualArterialRelation(kind:ResolvedRelation['kind'],source:Part,target:Part,parts:Part[]):boolean{
 if(!parts.includes(source)||!parts.includes(target))return false;
 return !!lingualArterialRelations(parts).get(source.id)?.some(row=>row.kind===kind&&row.target===target);
}
