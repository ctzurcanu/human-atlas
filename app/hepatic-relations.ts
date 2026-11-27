import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';
import {isInstalledLiverLobe} from './liver-controls';
import {hepaticConnectionRole,isHepaticNativeParent} from './hepatic-connection-controls';

export const isNativeWholeLiver=(p:Part)=>!p.suppressed&&p.system==='digestive'&&(
 p.id==='ZA:Liver:Organ'&&p.conceptId==='ZA:Liver'&&p.name==='Liver · Organ'&&p.sourceId==='Liver'
 ||p.id==='LOCAL:female:liver:1032'&&p.conceptId==='LOCAL:female:liver'&&p.name==='Liver'&&p.sourceId==='liver');
const isVein=(p:Part,name:string)=>!p.suppressed&&p.system==='venous'&&(p.name===name||name==='Hepatic veins'&&/^TA98HEPATIC:(male|female):hepatic_vein$/.test(p.id)&&hepaticConnectionRole(p)==='hepatic_vein')&&!!nativePartFrame(p);
export const isNativePortalVein=(p:Part)=>!p.suppressed&&p.system==='venous'&&!!nativePartFrame(p)&&(
 /^TA98HEPATIC:(male|female):portal$/.test(p.id)&&hepaticConnectionRole(p)==='portal'
 ||p.id==='ZA:Hepatic portal vein'&&p.conceptId==='ZA:Hepatic portal vein'&&p.sourceId==='Hepatic portal vein'&&p.name==='Hepatic portal vein'
 ||p.id==='LOCAL:female:portal_vein:1795'&&p.conceptId==='LOCAL:female:portal_vein'&&p.sourceId==='portal_vein'&&p.name==='Portal vein');
export function hepaticSourceLimitation(p:Part):string|undefined{
 if(isInstalledLiverLobe(p))return p.coverageLimitation;
 if(isNativeWholeLiver(p))return 'Whole native liver exterior. Lobular microcirculation, sinusoids, portal triads, segmental vascular territories and vessel openings are not established by this surface. Related blood-flow routes describe anatomy; they do not certify mesh junctions or lumen continuity.';
}
/** Keep the named right hepatic branch distinct from a search result for
 * the common hepatic artery. This source does not establish a cystic outlet. */
export function hepaticDescription(p:Part):{text:string;url:string}|undefined{
 if(p.suppressed||p.system!=='arterial'||!nativePartFrame(p))return;
 if(p.id==='LOCAL:female:r_hepatic_artery:1719'&&p.conceptId==='LOCAL:female:r_hepatic_artery'&&p.sourceId==='r_hepatic_artery'&&p.name==='Hepatic artery (right)')
  return {text:'The right hepatic artery supplies the liver and is the usual origin of the cystic artery supplying the gallbladder. Cystic branch origins vary. This named source is distinct from the common hepatic artery; its local cystic branch and physical vessel junction remain unverified.',url:'https://www.ncbi.nlm.nih.gov/sites/books/NBK459288/'};
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Liver Anatomy, PMC4038911: portal inflow mixes in sinusoids; hepatic
 * veins drain to the IVC. Never classify portal inflow as hepatic drainage. */
export function hepaticRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(a:Part,b:Part,forward:'after'|'venous',reverse:'before'|'drains',note:string,via?:string[])=>{
  for(const [from,to,kind] of [[a,b,forward],[b,a,reverse]] as const){const links=graph.get(from.id)??[];links.push({kind,target:to,note,...(via?{via}: {})});graph.set(from.id,links);}
 };
 for(const liver of parts.filter(p=>isNativeWholeLiver(p)||isInstalledLiverLobe(p))){
  const frame=nativePartFrame(liver);
  for(const portal of parts.filter(p=>isNativePortalVein(p)&&nativePartFrame(p)===frame))add(portal,liver,'after','before','Portal venous inflow enters the liver and mixes with arterial inflow in hepatic sinusoids. This is a physiological route to the whole organ, not a direct portal-to-hepatic-vein shunt or proof of source vessel openings.');
  for(const veins of parts.filter(p=>isVein(p,'Hepatic veins')&&nativePartFrame(p)===frame))add(liver,veins,'venous','drains','Hepatic veins carry liver blood to the systemic circulation after sinusoidal and central-vein collection. Separate microvascular collectors, segmental territories and patent mesh junctions remain unverified.',['Hepatic sinusoids','Central veins']);
 }
 for(const veins of parts.filter(p=>isVein(p,'Hepatic veins'))){
  const frame=nativePartFrame(veins);
  for(const cava of parts.filter(p=>!p.suppressed&&nativePartFrame(p)===frame&&(p.system==='venous'||p.system==='cardiac')&&['Inferior vena cava','Inferior vena cava (abdominal part)'].includes(p.name)))add(veins,cava,'after','before','Hepatic venous outflow enters the inferior vena cava. Named native blood-flow association; individual ostia, vessel walls and source lumen continuity remain unverified.');
 }
 const put=(a:Part,b:Part,kind:ResolvedRelation['kind'],note:string)=>{const links=graph.get(a.id)??[];if(!links.some(r=>r.kind===kind&&r.target.id===b.id))links.push({kind,target:b,note});graph.set(a.id,links);};
 const arterial=parts.find(p=>/^TA98HEPATIC:(male|female):artery$/.test(p.id)&&hepaticConnectionRole(p)==='artery');
 if(arterial)for(const liver of parts.filter(p=>isInstalledLiverLobe(p))){put(liver,arterial,'arterial','Macroscopic hepatic arterial branches. Arterioles, sinusoids and their parenchymal interfaces remain unmodelled.');put(arterial,liver,'supplies','Macroscopic arterial supply; no direct arterial-to-venous shunt is represented.');}
 for(const role of ['portal','biliary']as const){
  const root=parts.find(p=>hepaticConnectionRole(p)===role&&!p.id.includes(':region:'));
  const sex=root?.id.split(':')[1],parent=role==='portal'?`TA98HEPATIC:${sex}:portal-distal`:sex==='female'?'TA98HEPATIC:female:common-hepatic-parent':'TA98PART:FJ3079:A05.8.01.061';
  const lower=parts.find(p=>p.id===parent&&isHepaticNativeParent(p));
  if(root&&lower){put(root,lower,'continuous','Repaired communicating network meets the retained native downstream parent.');put(lower,root,'continuous','Repaired communicating network meets the retained native downstream parent.');}
 }
 const lower=parts.find(p=>p.id==='TA98PART:FJ3079:A05.8.01.061'&&isHepaticNativeParent(p));
 if(lower)for(const id of ['FJ3080','TA98PART:FJ3079:A05.8.02.013']){const peer=parts.find(p=>p.id===id&&!p.suppressed);if(peer){put(lower,peer,'continuous','Retained native cystic/bile interface.');put(peer,lower,'continuous','Retained native cystic/bile interface.');}}
 cache.set(parts,graph);return graph;
}
