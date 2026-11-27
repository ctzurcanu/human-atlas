import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';

export type TonguePoint=[number,number,number];
export type TongueAnchor={partId:string;face:number;barycentric:TonguePoint;point:TonguePoint};
export type TongueCurve={partId:string;points:TonguePoint[];closed:boolean;sourceBankSha256:string};
export type TongueTermGroup={key:string;name:string;term:string;conceptId:string;parent?:string;partIds:string[];kind:'exact';representation?:'surface-region'|'point-landmark'|'boundary-curve';surface?:{partId:string;faces:number[];sourceBankSha256:string};anchor?:TongueAnchor;curves?:TongueCurve[]};
export type TongueIntegration={schema:'tongue-local/v1';sex?:'male'|'female';wholeConceptId:string;partIds:string[];termGroups:TongueTermGroup[];bankSha256?:string;bankUrl?:string;carrierBankHashes?:Record<string,string>;legacyWholeAlias?:boolean;sourcePartId?:string};
type Registry={marker:TongueIntegration;parts:Part[];byId:Map<string,Part>};
const cache=new WeakMap<Atlas,Registry|null>();
const point=(p:unknown):p is TonguePoint=>Array.isArray(p)&&p.length===3&&p.every(Number.isFinite);
const sha=(s:unknown):s is string=>typeof s==='string'&&/^[a-f0-9]{64}$/.test(s);
const terms=new Set(['A04.1.05.001',...Array.from({length:25},(_,i)=>`A05.1.04.${String(i+1).padStart(3,'0')}`),...Array.from({length:10},(_,i)=>`A05.1.04.${101+i}`)]);
const equalIds=(a:string[],b:string[])=>a.length===b.length&&a.every(id=>b.includes(id));

/** Exact identities and anchored references are controls, not a claim that
 * anatomical boundaries or muscle attachments have been accepted. */
export function tongueRegistry(atlas:Atlas):Registry|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const marker=(atlas as Atlas&{finiteTongueIntegration?:TongueIntegration}).finiteTongueIntegration;
 if(!marker||marker.schema!=='tongue-local/v1'||!Array.isArray(marker.partIds)||!marker.partIds.length||new Set(marker.partIds).size!==marker.partIds.length||!Array.isArray(marker.termGroups)||!marker.termGroups.length||marker.termGroups.length>36)return reject();
 if(marker.bankSha256!==undefined&&!sha(marker.bankSha256)||marker.carrierBankHashes&&Object.values(marker.carrierBankHashes).some(hash=>!sha(hash)))return reject();
 if(marker.termGroups.some(g=>!g||typeof g!=='object'||!Array.isArray(g.partIds)))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const owned=new Set(marker.partIds),groups=new Map(marker.termGroups.map(g=>[g.key,g]));
 if(marker.partIds.some(id=>!byId.has(id)||byId.get(id)!.suppressed||byId.get(id)!.indexCount<3)||groups.size!==marker.termGroups.length||new Set(marker.termGroups.map(g=>g.term)).size!==marker.termGroups.length||new Set(marker.termGroups.map(g=>g.conceptId)).size!==marker.termGroups.length)return reject();
 const whole=marker.termGroups.find(g=>g.conceptId===marker.wholeConceptId);
 if(!whole||whole.term!=='A05.1.04.001'||whole.parent||!equalIds(whole.partIds,marker.partIds))return reject();
 const binding=(id:string,hash:string)=>{
  const p=byId.get(id),expected=marker.carrierBankHashes?.[id]??marker.bankSha256;
  return !!p&&sha(hash)&&sha(expected)&&hash===expected&&(!marker.bankUrl||!!marker.carrierBankHashes?.[id]||atlas.chunks[p.chunk]?.url===marker.bankUrl);
 };
 const anchorValid=(a:TongueAnchor|undefined,ids:string[])=>!!a&&ids.includes(a.partId)&&point(a.point)&&point(a.barycentric)&&a.barycentric.every(w=>w>=0&&w<=1)&&Math.abs(a.barycentric.reduce((n,w)=>n+w,0)-1)<1e-10&&Number.isSafeInteger(a.face)&&a.face>=0&&a.face<byId.get(a.partId)!.indexCount/3;
 for(const g of marker.termGroups){
  const c=concepts.get(g.conceptId);
  if(!g.key||!g.name||!terms.has(g.term)||g.kind!=='exact'||!Array.isArray(g.partIds)||!g.partIds.length||new Set(g.partIds).size!==g.partIds.length||g.partIds.some(id=>!owned.has(id))||!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||!equalIds(c.elements,g.partIds)||g.parent&&!groups.has(g.parent))return reject();
  const visited=new Set([g.key]);let parent=g.parent;
  while(parent){if(visited.has(parent))return reject();visited.add(parent);parent=groups.get(parent)?.parent;}
  if(g.anchor&&!anchorValid(g.anchor,g.partIds))return reject();
  if(g.representation&&!['surface-region','point-landmark','boundary-curve'].includes(g.representation))return reject();
  if(g.representation&&!anchorValid(g.anchor,g.partIds))return reject();
  if(g.representation&&(g.partIds.length!==1||g.partIds[0]!==g.anchor!.partId))return reject();
  if(g.representation==='surface-region'){
   const s=g.surface;
   if(!s||!g.partIds.includes(s.partId)||!binding(s.partId,s.sourceBankSha256)||!Array.isArray(s.faces)||!s.faces.length||new Set(s.faces).size!==s.faces.length||s.faces.some(face=>!Number.isSafeInteger(face)||face<0||face>=byId.get(s.partId)!.indexCount/3)||g.anchor!.partId!==s.partId||!s.faces.includes(g.anchor!.face))return reject();
  }
  if(g.representation==='boundary-curve'&&(!Array.isArray(g.curves)||!g.curves.length))return reject();
  if(g.curves?.some(curve=>!g.partIds.includes(curve.partId)||!binding(curve.partId,curve.sourceBankSha256)||typeof curve.closed!=='boolean'||!Array.isArray(curve.points)||curve.points.length<(curve.closed?3:2)||curve.points.length>16384||!curve.points.every(point)))return reject();
  // A single selection overlay follows one carrier rigidly. Multiple carrier
  // curves need separate groups rather than an undeclared frame conversion.
  if(g.curves?.some(curve=>curve.partId!==g.anchor?.partId))return reject();
 }
 const result={marker,parts:marker.partIds.map(id=>byId.get(id)!),byId};cache.set(atlas,result);return result;
}

export function tongueTermParts(atlas:Atlas):Map<string,Part[]>{
 const registry=tongueRegistry(atlas);
 return new Map(registry?.marker.termGroups.map(g=>[g.term,g.partIds.map(id=>registry.byId.get(id)!)])??[]);
}
export function selectedTongueFeature(atlas:Atlas,choice:HierarchyChoice|undefined):TongueTermGroup|undefined{
 const registry=tongueRegistry(atlas);if(!registry||!choice)return;
 const {marker}=registry,groups=marker.termGroups;
 const exact=groups.find(g=>g.conceptId===choice.id&&equalIds(g.partIds,choice.elements));
 if(exact)return exact;
 const native=atlas.concepts.find(c=>c.id===choice.id);
 if(native){
  // Only an explicitly declared, suppressed original source may redirect its
  // legacy concept. Native concepts often have no TA98 field; term equality is
  // neither necessary nor sufficient to establish this source-bound alias.
  const whole=groups.find(g=>g.conceptId===marker.wholeConceptId),source=marker.sourcePartId?registry.byId.get(marker.sourcePartId):undefined;
  if(marker.legacyWholeAlias===true&&source?.suppressed===true
    &&source.conceptId===native.id&&whole&&equalIds(native.elements,marker.partIds)
    &&equalIds(native.elements,whole.partIds)&&equalIds(choice.elements,native.elements)
    &&(!native.ta98Term||native.ta98Term===whole.term)&&(!native.ta98Kind||native.ta98Kind==='exact'))return whole;
  return;
 }
 return groups.find(g=>g.term===(choice.ta98Term??choice.terminology?.ta98)&&equalIds(g.partIds,choice.elements));
}
export function selectedTongueSurfaceFeature(atlas:Atlas,choice:HierarchyChoice|undefined){
 const g=selectedTongueFeature(atlas,choice);
 return g?.representation?g:undefined;
}
/** Both fibre populations occupy the same carrier; selection changes its
 * directional illustration without duplicating tissue or moving the camera. */
export function tongueFibreWeights(atlas:Atlas,choice:HierarchyChoice|undefined):[number,number]{
 const registry=tongueRegistry(atlas),groups=registry?.marker.termGroups;
 const transverse=groups?.find(g=>g.term==='A05.1.04.108'),vertical=groups?.find(g=>g.term==='A05.1.04.109');
 if(!transverse||!vertical||!equalIds(transverse.partIds,vertical.partIds))return [.5,.5];
 const selected=selectedTongueFeature(atlas,choice);
 return selected===transverse?[1,0]:selected===vertical?[0,1]:[.5,.5];
}
export function tongueLabelReference(atlas:Atlas,choice:HierarchyChoice){
 const g=selectedTongueFeature(atlas,choice);return g?.anchor?{partId:g.anchor.partId,point:g.anchor.point}:undefined;
}
export function tongueComponentChoices(atlas:Atlas,choice:HierarchyChoice|undefined):HierarchyChoice[]|undefined{
 const registry=tongueRegistry(atlas),root=selectedTongueFeature(atlas,choice);if(!registry||!root||!choice)return;
 const concepts=new Map(atlas.concepts.map(c=>[c.id,c])),groups=registry.marker.termGroups;
 const children=(key:string)=>groups.filter(g=>g.parent===key||key===root.key&&root.conceptId===registry.marker.wholeConceptId&&!g.parent&&g.key!==key);
 const pending=children(root.key),result:HierarchyChoice[]=[];
 while(pending.length){const g=pending.shift()!;pending.push(...children(g.key));const c=concepts.get(g.conceptId)!;if(g.partIds.every(id=>choice.elements.includes(id)))result.push({...c,elements:[...g.partIds]});}
 return result.length?result:undefined;
}
export function tongueReferenceChoices(atlas:Atlas){
 const registry=tongueRegistry(atlas),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 return new Map(registry?.marker.termGroups.filter(g=>!!g.representation).map(g=>{
  const p=registry.byId.get(g.anchor!.partId)!;
  return [g.term,{choice:concepts.get(g.conceptId)!,host:{id:p.conceptId,name:p.name,elements:[p.id]}}];
 })??[]);
}
