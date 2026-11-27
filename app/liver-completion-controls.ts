import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {liverControlRegistry} from './liver-controls';
export type LiverWitness={partId:string;face:number;barycentric:[number,number,number];point:[number,number,number]};
export type LiverReferenceGroup={term:string;conceptId:string;partIds:string[];representation:string;anchor?:LiverWitness;curves?:{closed:boolean;witnesses:LiverWitness[]}[];detailModel?:string};
type Marker={schema:string;hostPartId:string;hostBankSha256:string;sourceCatalogueSha256:string;partIds:string[];functionalPartIds:string[];microscopicPartIds:string[];ligamentPartIds:string[];termGroups:LiverReferenceGroup[]};
const cache=new WeakMap<Atlas,Marker|null>();
export function liverCompletionRegistry(atlas:Atlas):Marker|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;},liver=liverControlRegistry(atlas);
 const m=(atlas as Atlas&{liverCompletionIntegration?:Marker}).liverCompletionIntegration;
 if(!m||!liver||m.schema!==`${atlas.sex}-liver-completion/v1`||m.hostPartId!=='TA98PART:ZA:Liver:capsule'||m.hostBankSha256!==liver.marker.bankSha256
  ||!Array.isArray(m.partIds)||m.partIds.length!==16||new Set(m.partIds).size!==16||m.functionalPartIds?.length!==8||m.microscopicPartIds?.length!==6||m.ligamentPartIds?.length!==2
  ||!Array.isArray(m.termGroups)||m.termGroups.length!==30||new Set(m.termGroups.map(g=>g.term)).size!==30)return reject();
 const parts=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 if(m.partIds.some(id=>!parts.has(id)||parts.get(id)!.suppressed||!parts.get(id)!.pairedTissue||!parts.get(id)!.sectionSolid))return reject();
 for(const g of m.termGroups){const c=concepts.get(g.conceptId);
  if(!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||!g.partIds.length||c.elements.length!==g.partIds.length||!c.elements.every(id=>g.partIds.includes(id))
   ||g.partIds.some(id=>!parts.has(id)||parts.get(id)!.suppressed))return reject();
  if(g.curves){if(g.partIds.length!==1||g.partIds[0]!==m.hostPartId||!g.anchor||!g.curves.length)return reject();
   for(const curve of g.curves){if(curve.witnesses.length<(curve.closed?3:2))return reject();
    for(const w of curve.witnesses)if(w.partId!==m.hostPartId||!Number.isSafeInteger(w.face)||w.face<0||w.face>=parts.get(m.hostPartId)!.indexCount/3||w.point?.length!==3||!w.point.every(Number.isFinite)||w.barycentric?.length!==3||w.barycentric.some(x=>!Number.isFinite(x)||x<0||x>1)||Math.abs(w.barycentric.reduce((s,x)=>s+x,0)-1)>1e-12)return reject();
   }
  }
 }
 cache.set(atlas,m);return m;
}
export function liverCompletionTermParts(atlas:Atlas){
 const r=liverCompletionRegistry(atlas),parts=new Map(atlas.parts.map(p=>[p.id,p]));return new Map(r?.termGroups.map(g=>[g.term,g.partIds.map(id=>parts.get(id)!)])??[]);
}
export function liverWholePresentationParts(atlas:Atlas,term:string,parts:Part[]){
 const m=liverCompletionRegistry(atlas);if(!m||term!=='A05.8.01.001')return parts;
 const alternatives=new Set([...m.functionalPartIds,...m.microscopicPartIds]);return parts.filter(p=>!alternatives.has(p.id));
}
export function selectedLiverBoundary(atlas:Atlas,choice:HierarchyChoice|undefined){
 const m=liverCompletionRegistry(atlas);if(!m||!choice||choice.elements.length!==1||choice.elements[0]!==m.hostPartId)return;
 return m.termGroups.find(g=>g.curves&&(g.conceptId===choice.id||g.term===(choice.ta98Term??choice.terminology?.ta98)));
}
export function liverCompletionLabelReference(atlas:Atlas,choice:HierarchyChoice){return selectedLiverBoundary(atlas,choice)?.anchor;}
export function liverMicroscopicDetailUrl(choice:HierarchyChoice){
 const names=['perivascular-sheath','lobule','arteries','portal-venules','central-vein','bile-ductules'];
 const term=choice.ta98Term??choice.terminology?.ta98??choice.id.match(/A05\.8\.01\.0(?:55|56|57|58|59|60)$/)?.[0];
 const n=Number(term?.split('.').at(-1));
 return n>=55&&n<=60&&choice.elements.length===1&&choice.elements[0]==='TA98LIVER:micro:'+names[n-55]
  ?`/?model=liver-unit&select=${encodeURIComponent('TA98LIVER:male:term:'+term)}&skin=0&context=0&isolate=1&focus=1&view=front`:null;
}
