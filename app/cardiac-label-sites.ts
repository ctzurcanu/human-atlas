import type {Atlas,Part,Concept} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {cardiacIntegration,cardiacTermParts,type CardiacTermGroup} from './heart-wall';
import {validCardiacReferenceGeometry,type CardiacReferenceGeometry} from './cardiac-reference-types';
import {hasSourceCardiacIntegration,sourceCardiacIntegration,sourceCardiacTermParts} from './biodigital-heart';

type Point=[number,number,number];
export type CardiacLabelSite=CardiacTermGroup&{
 coordinateFrame?:'atlas-world-metres';
 referenceOnly:true;
 representation:'anatomical-point'|'regional-reference'|'collective'|CardiacReferenceGeometry['type'];
 anchor:{partId:string;point:Point};
 points:Point[]|null;
 referenceCids:number[];
 tissueGeometryEstablished:false;
 geometry?:CardiacReferenceGeometry;
};
type Pose={geometryOrigin:number[];geometryInstance:unknown;chunk:number;positions:number;indices:number;indexCount:number};
type Registry={schema:'cardiac-qualified-label-sites/v1';sex:'male'|'female';bankSha256:string;referencePlanSha256:string;termGroups:CardiacLabelSite[];carrierPoses:Record<string,Pose>;geometryBackedTermCount:number;referenceTermCount:number;labeledTermCount:number;distinctTissueCoverageUnchanged:true;typedGeometryCount?:number;typedGeometryCounts?:Record<string,number>;typedReferenceRevision?:string};
type Bundle={schema:'cardiac-qualified-label-bundle/v1';integration:Registry;concepts:Concept[]};
const cache=new WeakMap<Atlas,Registry|null>();
const kinds=new Set(['anatomical-point','opening-rim','boundary-curve','regional-reference','collective','finite-mesh','surface-region','point-set','space-region','collection']);
const point=(p:unknown):p is Point=>Array.isArray(p)&&p.length===3&&p.every(Number.isFinite);
const sha=(s:unknown)=>typeof s==='string'&&/^[a-f0-9]{64}$/.test(s);
const pose=(p:Part,sourceChunks?:number[]):Pose=>({geometryOrigin:p.geometryOrigin??[0,0,0],geometryInstance:p.geometryInstance??null,chunk:sourceChunks?sourceChunks[p.chunk]:p.chunk,positions:p.positions,indices:p.indices,indexCount:p.indexCount});

/** Typed reference geometry remains distinct from established tissue geometry.
 * Opening rims stay empty, and illustrative microscopic extents stay qualified. */
export function cardiacLabelRegistry(atlas:Atlas):Registry|null{
 if(hasSourceCardiacIntegration(atlas)){
  const source=sourceCardiacIntegration(atlas);if(!source)return null;
  return {schema:'cardiac-qualified-label-sites/v1',sex:source.sex,bankSha256:source.bankSha256,referencePlanSha256:source.sourceManifestSha256,termGroups:source.referenceSites,carrierPoses:source.carrierPoses,geometryBackedTermCount:source.termGroups.length-source.referenceSites.length,referenceTermCount:source.referenceSites.length,labeledTermCount:source.termGroups.length,distinctTissueCoverageUnchanged:true};
 }
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;},heart=cardiacIntegration(atlas);
 const m=(atlas as Atlas&{cardiacLabelIntegration?:Registry}).cardiacLabelIntegration;
 const compact=atlas as Atlas&{sourceChunkIndices?:number[];sourceChunkCount?:number},sourceChunks=compact.sourceChunkIndices;
 if(sourceChunks!==undefined&&(!Array.isArray(sourceChunks)||sourceChunks.length!==atlas.chunks.length||!Number.isSafeInteger(compact.sourceChunkCount)||compact.sourceChunkCount!<atlas.chunks.length||sourceChunks.some(i=>!Number.isSafeInteger(i)||i<0||i>=compact.sourceChunkCount!)||new Set(sourceChunks).size!==sourceChunks.length))return reject();
 if(!heart||!m||m.schema!=='cardiac-qualified-label-sites/v1'||m.sex!==heart.sex||m.bankSha256!==heart.coupledBankSha256||!sha(m.referencePlanSha256)||!Number.isSafeInteger(m.geometryBackedTermCount)||!Number.isSafeInteger(m.referenceTermCount)||m.referenceTermCount<1||m.geometryBackedTermCount<1||m.labeledTermCount!==m.geometryBackedTermCount+m.referenceTermCount||m.labeledTermCount>113||m.distinctTissueCoverageUnchanged!==true||m.termGroups.length!==m.referenceTermCount)return reject();
 const parts=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c])),body=cardiacTermParts(atlas);
 const owned=new Set([...heart.wallPartIds,...heart.valvePartIds,...heart.papillaryPartIds,...heart.chordalPartIds??[],...heart.rootPartIds??[]]);
 if(body.size!==m.geometryBackedTermCount||new Set(m.termGroups.map(g=>g.term)).size!==m.referenceTermCount||Object.entries(m.carrierPoses).some(([id,expected])=>!parts.has(id)||Object.entries(pose(parts.get(id)!,sourceChunks)).some(([key,value])=>JSON.stringify(value)!==JSON.stringify(expected[key as keyof Pose]??null))))return reject();
 for(const g of m.termGroups){const c=concepts.get(g.conceptId);
  if(!/^A12\.1\.\d{2}\.\d{3}$/.test(g.term)||body.has(g.term)||g.referenceOnly!==true||g.tissueGeometryEstablished!==false||!kinds.has(g.representation)||!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||!g.partIds.length||g.partIds.some(id=>!owned.has(id)||!parts.has(id)||parts.get(id)!.suppressed)||JSON.stringify(c.elements)!==JSON.stringify(g.partIds)||!g.partIds.includes(g.anchor.partId)||!m.carrierPoses[g.anchor.partId]||!point(g.anchor.point)||g.points!==null&&(!g.points.length||!g.points.every(point)))return reject();
  if(g.geometry&&!validCardiacReferenceGeometry(g.geometry,parts,m.bankSha256))return reject();
 }
 if(m.typedGeometryCount!==undefined){
  const counts:Record<string,number>={};for(const g of m.termGroups){if(!g.geometry)return reject();counts[g.geometry.type]=(counts[g.geometry.type]??0)+1;}
  if(m.typedGeometryCount!==m.referenceTermCount||!m.typedGeometryCounts||Object.keys(counts).length!==Object.keys(m.typedGeometryCounts).length||Object.entries(counts).some(([type,count])=>m.typedGeometryCounts![type]!==count))return reject();
 }
 cache.set(atlas,m);return m;
}

/** Label coverage, including qualified reference sites; never use this as a
 * count of separately built tissues. */
export function cardiacLabeledTermParts(atlas:Atlas):Map<string,Part[]>{
 if(hasSourceCardiacIntegration(atlas))return sourceCardiacTermParts(atlas);
 const result=cardiacTermParts(atlas),m=cardiacLabelRegistry(atlas);if(!m)return result;
 const byId=new Map(atlas.parts.map(p=>[p.id,p]));for(const g of m.termGroups)result.set(g.term,g.partIds.map(id=>byId.get(id)!));return result;
}

/** Logical component entries are distinct from their shared physical carriers.
 * Expand every installed TA98 descendant once, retaining reference limitations. */
export function cardiacComponentChoices(atlas:Atlas,choice:HierarchyChoice|undefined):HierarchyChoice[]|undefined{
 if(hasSourceCardiacIntegration(atlas)&&!sourceCardiacIntegration(atlas))return;
 const source=sourceCardiacIntegration(atlas);
 const heart=cardiacIntegration(atlas);if(!(source?.termGroups??heart?.termGroups)||!choice)return;
 // Selecting both the cardiac Heart and its donor Heart group produces a
 // selection set with the same complete membership as that donor assembly.
 // Preserve its logical heart details instead of falling back to 24 carriers.
 if(choice.id==='selection-set'){
  if(!source)return;
  const ids=new Set(choice.elements),whole=source.termGroups.find(g=>g.term==='A12.1.00.001');
  const candidates=[whole&&atlas.concepts.find(c=>c.id===whole.conceptId),atlas.concepts.find(c=>c.id===`BDM:${source.sex}:human_20_${source.sex}_cardiovascular_system-heart_ID`)];
  const assembly=candidates.find(c=>c&&c.elements.length===ids.size&&c.elements.every(id=>ids.has(id)));
  return assembly?cardiacComponentChoices(atlas,assembly):undefined;
 }
 const groups=source?.termGroups??[...heart!.termGroups!,...cardiacLabelRegistry(atlas)?.termGroups??[]];
 const term=choice.ta98Term??choice.terminology?.ta98;
 const concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 // The original donor Heart group also contains its coronary vessels. Expand
 // the validated cardiac terms on their actual carriers, then those original
 // vascular leaves; otherwise this label silently loses every typed site.
 const donorId=source&&`BDM:${source.sex}:human_20_${source.sex}_cardiovascular_system-heart_ID`;
 const donor=donorId&&choice.id===donorId?concepts.get(donorId):undefined;
 const donorChoice=donor&&choice.elements.length===donor.elements.length&&new Set(choice.elements).size===choice.elements.length&&choice.elements.every(id=>donor.elements.includes(id));
 const root=groups.find(g=>g.conceptId===choice.id)??groups.find(g=>g.kind==='exact'&&g.term===term)??(donorChoice?groups.find(g=>g.term==='A12.1.00.001'):undefined);
 if(!root||root.partIds.some(id=>!choice.elements.includes(id)))return;
 const available=cardiacLabeledTermParts(atlas);
 const result=new Map<string,HierarchyChoice>(),pending=groups.filter(g=>g.parent===root.key);
 while(pending.length){const g=pending.shift()!;pending.push(...groups.filter(child=>child.parent===g.key));
  const c=concepts.get(g.conceptId);
  if(g.kind==='exact'&&g.term!==root.term&&available.has(g.term)&&c&&!result.has(g.term)&&g.partIds.every(id=>choice.elements.includes(id)))result.set(g.term,{...c,elements:[...g.partIds]});
 }
 if(donorChoice){
  const cardiacIds=new Set(root.partIds),parts=new Map(atlas.parts.map(p=>[p.id,p]));
  for(const id of donor.elements){
   if(cardiacIds.has(id))continue;
   const p=parts.get(id),c=p&&concepts.get(p.conceptId);
   if(!p||p.suppressed||p.defaultHidden||atlas.defaultHidden?.includes(id)||!c)continue;
   const key=`source-leaf:${c.id}`,prior=result.get(key);
   if(prior){if(!prior.elements.includes(id))prior.elements.push(id);}
   else result.set(key,{...c,elements:[id]});
  }
 }
 // A reference uses carrier meshes for picking and highlighting. Those whole
 // chambers are not anatomical components contained in a surface or opening.
 return result.size||source?.referenceSites.some(s=>s.term===root.term)?[...result.values()]:undefined;
}
export function selectedCardiacLabelSite(atlas:Atlas,choice:HierarchyChoice|undefined){
 const m=cardiacLabelRegistry(atlas);if(!m||!choice)return;
 return m.termGroups.find(g=>(g.conceptId===choice.id||g.term===(choice.ta98Term??choice.terminology?.ta98))&&g.partIds.length===choice.elements.length&&g.partIds.every(id=>choice.elements.includes(id)));
}
export function cardiacSiteLabelReference(atlas:Atlas,choice:HierarchyChoice){
 const g=selectedCardiacLabelSite(atlas,choice);return g?{partId:g.anchor.partId,point:g.anchor.point}:undefined;
}
export function cardiacSiteCoverage(atlas:Atlas,choice:HierarchyChoice){
 const g=selectedCardiacLabelSite(atlas,choice);return g?atlas.concepts.find(c=>c.id===g.conceptId)?.coverageLimitation:undefined;
}
export function cardiacSiteConceptCoverage(concepts:Concept[],choice:HierarchyChoice){
 return concepts.find(c=>/^(?:TA98HEART:(male|female):label-site|BDH:(male|female):adapted-reference):A12\.1\./.test(c.id)&&(c.id===choice.id||c.ta98Term===(choice.ta98Term??choice.terminology?.ta98))&&c.elements.length===choice.elements.length&&c.elements.every(id=>choice.elements.includes(id)))?.coverageLimitation;
}

/** Install typed reference data. Parts, chunks, materials and camera state
 * remain untouched; a stale carrier pose rejects the entire bundle. */
export function withCardiacLabelBundle(atlas:Atlas,bundle:Bundle):Atlas{
 if(hasSourceCardiacIntegration(atlas))return atlas;
 if(bundle?.schema!=='cardiac-qualified-label-bundle/v1'||!bundle.integration||!Array.isArray(bundle.concepts)||bundle.concepts.length!==bundle.integration.referenceTermCount||new Set(bundle.concepts.map(c=>c.id)).size!==bundle.concepts.length)return atlas;
 const ids=new Set(bundle.concepts.map(c=>c.id));
 if(atlas.concepts.some(c=>ids.has(c.id)))return atlas;
 const candidate={...atlas,concepts:[...atlas.concepts,...bundle.concepts],cardiacLabelIntegration:bundle.integration};
 return cardiacLabelRegistry(candidate)?candidate:atlas;
}
export async function loadLocalCardiacLabels(atlas:Atlas,development:boolean,init:RequestInit):Promise<Atlas>{
 if(hasSourceCardiacIntegration(atlas))return atlas;
 const heart=cardiacIntegration(atlas);if(!development||!heart)return atlas;
 try{
  const stage=typeof window==='undefined'?null:new URLSearchParams(window.location.search).get('cardiacRefStage');
  const path=stage&&/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(stage)?`/local-models/geometry-workbench/heart/${stage}/${heart.sex}/cardiac-label-sites.json`:`/local-models/ta98-review/bowel-wall-local-atlas/${heart.sex}/cardiac-label-sites.json`;
  const response=await fetch(path,{...init,cache:'no-store'});
  return response.ok?withCardiacLabelBundle(atlas,await response.json()):atlas;
 }catch(error){if(init.signal?.aborted)throw error;return atlas;}
}
