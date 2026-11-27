import type {Atlas,Part,Concept} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';

export type HepaticConnectionRole='artery'|'portal'|'hepatic_vein'|'biliary';
type ConnectionPart=Part&{hepaticConnectionRole?:HepaticConnectionRole;hepaticRegionOf?:string;hepaticSourceModel?:string;hepaticAttributeSignature?:string;hepaticRetainedParent?:boolean};
type Group={term:string;conceptId:string;parts:string[]};
type Marker={schema:string;wholePartIds:string[];regionPartIds:string[];termGroups:Group[];sourceCatalogueSha256:string};
type Registry={marker:Marker;byId:Map<string,Part>};
const roles:HepaticConnectionRole[]=['artery','portal','hepatic_vein','biliary'];
const cache=new WeakMap<Atlas,Registry|null>();
const wholeId=(role:HepaticConnectionRole,sex='male')=>`TA98HEPATIC:${sex}:${role}`;
const sourceSex=(p:ConnectionPart)=>p.hepaticSourceModel==='local-ta98'?'male':p.hepaticSourceModel==='local-female-ta98'?'female':undefined;

/** Finite communicating networks and explicitly owned surface regions are
 * different representations. A region is not a new microscopic wall layer. */
export function hepaticConnectionRegistry(atlas:Atlas):Registry|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const m=(atlas as Atlas&{hepaticConnectionIntegration?:Marker}).hepaticConnectionIntegration;
 const sex=atlas.sex;
 if(!['male','female'].includes(sex??'')||!m||m.schema!==`${sex}-hepatic-connections/v1`||!/^[a-f0-9]{64}$/.test(m.sourceCatalogueSha256)
  ||m.wholePartIds?.length!==4||!roles.every(r=>m.wholePartIds.includes(wholeId(r,sex)))
  ||!Array.isArray(m.regionPartIds)||new Set(m.regionPartIds).size!==m.regionPartIds.length||!Array.isArray(m.termGroups))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 for(const role of roles){const p=byId.get(wholeId(role,sex))as ConnectionPart|undefined;
  if(!p||p.suppressed||p.defaultHidden||!p.pairedTissue||p.normalSigns===undefined||p.hepaticConnectionRole!==role
   ||sourceSex(p)!==sex||p.tissueWallAssembly!==`hepatic-connections:${sex}:${role}`)return reject();
 }
 for(const id of m.regionPartIds){const p=byId.get(id)as ConnectionPart|undefined,source=p?.hepaticRegionOf&&byId.get(p.hepaticRegionOf);
  if(!p||p.suppressed||!source||!hepaticRegionReplaces(source,p))return reject();
 }
 for(const g of m.termGroups){const c=concepts.get(g.conceptId);
  if(!c||c.ta98Term!==g.term||!g.parts.length||g.parts.some(id=>!byId.has(id)||byId.get(id)!.suppressed))return reject();
 }
 const r={marker:m,byId};cache.set(atlas,r);return r;
}

export function hepaticConnectionTermParts(atlas:Atlas):Map<string,Part[]>{
 const r=hepaticConnectionRegistry(atlas);return new Map(r?.marker.termGroups.map(g=>[g.term,g.parts.map(id=>r.byId.get(id)!)])??[]);
}

export function hepaticConnectionRole(p:Part):HepaticConnectionRole|undefined{
 const q=p as ConnectionPart,role=q.hepaticConnectionRole,sex=sourceSex(q);
 if(!role||!roles.includes(role)||!sex||p.suppressed||p.sourceOffset)return;
 if(p.id!==wholeId(role,sex)&&!new RegExp(`^TA98HEPATIC:${sex}:region:A\\d{2}\\.\\d\\.\\d{2}\\.\\d{3}$`).test(p.id))return;
 return role;
}

export function isHepaticNativeParent(p:Part):boolean{
 const q=p as ConnectionPart;
 const sex=sourceSex(q);
 return !p.suppressed&&q.hepaticRetainedParent===true&&!!sex
  &&(sex==='male'?['TA98PART:FJ3079:A05.8.01.061','TA98HEPATIC:male:portal-distal']:['TA98HEPATIC:female:common-hepatic-parent','TA98HEPATIC:female:portal-distal']).includes(p.id);
}

/** These are identical indexed faces, not opposing anatomical neighbours.
 * Replacement remains exclusive when an explicit Explode moves the region. */
export function hepaticRegionReplaces(source:Part,target:Part):boolean{
 const p=target as ConnectionPart,s=source as ConnectionPart,role=hepaticConnectionRole(source);
 return !!role&&p.hepaticRegionOf===source.id&&hepaticConnectionRole(target)===role
  &&source.vertexCount===target.vertexCount&&!!s.hepaticAttributeSignature&&/^[a-f0-9]{64}$/.test(s.hepaticAttributeSignature)
  &&s.hepaticAttributeSignature===p.hepaticAttributeSignature&&source.normalSigns!==undefined&&target.normalSigns!==undefined
  &&JSON.stringify(source.geometryOrigin)===JSON.stringify(target.geometryOrigin)&&(source.geometryScale??1)===(target.geometryScale??1);
}

export function hepaticComponentIds(choice:HierarchyChoice,parts:Map<string,Part>,concepts:Concept[]):string[]|undefined{
 if(choice.children||choice.elements.length!==1)return;
 const source=parts.get(choice.elements[0]);if(!source||!hepaticConnectionRole(source))return;
 const c=concepts.find(c=>c.id===choice.id)as(Concept&{hepaticComponents?:string[]})|undefined,ids=c?.hepaticComponents;
 if(!ids?.length||new Set(ids).size!==ids.length||ids.some(id=>!parts.has(id)||!hepaticRegionReplaces(source,parts.get(id)!)))return;
 return ids;
}
