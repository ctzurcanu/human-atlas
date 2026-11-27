import {catalogueSourceReference} from "./source-reference";
import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {AnatomyEntry} from './anatomy-hierarchy';
import {terminologyForConcept} from './anatomical-terminology';
import {nativePartFrame} from './ta98-modeled-parent';

type Point=[number,number,number];
const term='A12.2.04.008';
// The authority repeats sinus caroticus (FMA50094) beneath CCA and cervical
// ICA. Both occurrences refer to this one qualified left wall reference.
const sinusTerms=[term,'A12.2.06.003'] as const;
const hostId='TA98PART:ZA:Internal carotid artery.l:A12.2.06.002';
const featureId='TA98REFERENCE:male:carotid-sinus:left:medial-proximal-ICA';
const bankSha='e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36';
const bankUrl='/local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1/packed-left-v2/cells.bin';
const revision='left-four-cell-common-bank-v2';
const sourceUrl='/ATTRIBUTION.md#open-anatomy-source';
const hostLayout={positions:0,normals:195540,normalSigns:806208,uvs:293312,indices:480120,vertexCount:16295,indexCount:7848} as const;
const outerFirst=1270,outerEnd=2564;
export const carotidSurfaceReferenceName='Carotid sinus medial wall reference (left)';
export const carotidSurfaceReferenceDescription='Qualified reference on the available medial outer wall of the proximal left internal carotid artery. Boundaries are provisional and coverage of the first centimetre is partial. Sinus bulb morphology, lumen extent and baroreceptor histology remain unverified.';

/** The original illustrative macro-body colour is authored in unit RGB. The
 * viewer's source-material colours use byte RGB; keep that conversion local
 * to this registered organ rather than reinterpret dark colours elsewhere. */
export function carotidBodySurfaceFinish(part:Part):{color:[number,number,number];roughness:number;metalness:number}|undefined{
 if(part.id!=='TA98BODY:male:carotid-body:left:macro'||part.material!=='TA98BODY:carotid-body:illustrative-organ'||nativePartFrame(part)!=='ZA:')return;
 return {color:[165.75,107.1,71.4],roughness:.65,metalness:0};
}

export type CarotidSurfaceFeature={term:string;conceptId:string;partIds:string[];representation:'surface-region';referenceOnly:true;boundaryProvisional:true;surface:{partId:string;faces:number[];sourceBankSha256:string};anchor:{partId:string;face:number;barycentric:Point;point:Point}};
export type CarotidSurfaceMarker={schema:'carotid-surface-features/v1';side:'left';hostPartId:string;bankSha256:string;bankUrl:string;sourcePlanSha256:string;provisional:true;termGroups:CarotidSurfaceFeature[]};
type Host=Part&{carotidFiniteWallRevision?:string;carotidFiniteWallBankSha256?:string};
export type CarotidReferenceChoice={choice:HierarchyChoice;host:HierarchyChoice};
const cache=new WeakMap<Atlas,CarotidSurfaceMarker|null>();
const references=new WeakMap<Atlas,Map<string,CarotidReferenceChoice>>();
const sha=(value:unknown)=>typeof value==='string'&&/^[a-f0-9]{64}$/.test(value);
const point=(value:unknown):value is Point=>Array.isArray(value)&&value.length===3&&value.every(Number.isFinite);
const oneHost=(value:unknown)=>Array.isArray(value)&&value.length===1&&value[0]===hostId;

/** Display references on an existing carrier add no anatomical solid or
 * completed sinus geometry. Face identities bind this exact reviewed bank. */
export function carotidFeatureRegistry(atlas:Atlas):CarotidSurfaceMarker|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const m=(atlas as Atlas&{carotidSurfaceFeatureIntegration?:CarotidSurfaceMarker}).carotidSurfaceFeatureIntegration;
 const hosts=atlas.parts.filter(p=>p.id===hostId),host=hosts[0] as Host|undefined;
 if(atlas.sex!=='male'||!m||m.schema!=='carotid-surface-features/v1'||m.side!=='left'||m.hostPartId!==hostId
  ||m.bankSha256!==bankSha||m.bankUrl!==bankUrl||!sha(m.sourcePlanSha256)||m.provisional!==true
  ||!Array.isArray(m.termGroups)||m.termGroups.length!==1||hosts.length!==1||!host
  ||host.suppressed||host.defaultHidden||atlas.defaultHidden?.includes(host.id)
  ||host.system!=='arterial'||host.conceptId!==hostId||host.ta98Term!=='A12.2.06.002'||host.sourceId!=='Internal carotid artery.l'
  ||host.sourceOffset!==undefined||host.geometryInstance!==undefined||(host.geometryScale!==undefined&&host.geometryScale!==1)
  ||catalogueSourceReference(host.provenance?.url)!==sourceUrl||host.pairedTissue!==true||host.sectionSolid!==true
  ||host.sectionAssembly!=='left-carotid-lingual-coupled-finite-wall:male'
  ||host.carotidFiniteWallRevision!==revision||host.carotidFiniteWallBankSha256!==bankSha
  ||!point(host.geometryOrigin)||host.geometryOrigin.some((v,i)=>v!==[0,1.625,0][i])
  ||Object.entries(hostLayout).some(([key,value])=>host[key as keyof typeof hostLayout]!==value)
  ||!Number.isSafeInteger(host.chunk)||host.chunk<0||atlas.chunks[host.chunk]?.url!==bankUrl
  ||atlas.chunks[host.chunk]?.bytes!==822503||!Array.isArray(host.bounds)||host.bounds.length!==2
  ||!host.bounds.every(point)||host.bounds[0].some((v,i)=>v>host.bounds[1][i]))return reject();
 const groups=host.tissueDrawGroups;
 // Independent classifier replay identifies the SECOND chart as outer here.
 // The first chart is inner; the third is a shared tissue cap.
 if(!Array.isArray(groups)||groups.length!==3||groups.some(g=>!g||typeof g!=='object')
  ||groups[0].start!==0||groups[0].count!==3810||groups[0].neighbor!==undefined
  ||groups[1].start!==3810||groups[1].count!==3882||groups[1].neighbor!==undefined
  ||groups[2].start!==7692||groups[2].count!==156||groups[2].neighbor!=='TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall')return reject();
 const hostConcepts=atlas.concepts.filter(c=>c.id===hostId),featureConcepts=atlas.concepts.filter(c=>c.id===featureId);
 const hc=hostConcepts[0],fc=featureConcepts[0],g=m.termGroups[0];
 if(hostConcepts.length!==1||!hc||hc.ta98Term!=='A12.2.06.002'||hc.ta98Kind!=='exact'||!oneHost(hc.elements)
  ||featureConcepts.length!==1||!fc||fc.ta98Term!==term||fc.ta98Kind!=='exact'||!oneHost(fc.elements)
  ||!g||g.term!==term||g.conceptId!==featureId||!oneHost(g.partIds)||g.representation!=='surface-region'
  ||g.referenceOnly!==true||g.boundaryProvisional!==true)return reject();
 const s=g.surface,a=g.anchor;
 if(!s||s.partId!==hostId||s.sourceBankSha256!==bankSha||!Array.isArray(s.faces)||s.faces.length<1
  ||new Set(s.faces).size!==s.faces.length||s.faces.some(f=>!Number.isSafeInteger(f)||f<outerFirst||f>=outerEnd||f>=host.indexCount/3)
  ||!a||a.partId!==hostId||!Number.isSafeInteger(a.face)||!s.faces.includes(a.face)
  ||!point(a.point)||a.point.some((v,i)=>v<host.bounds[0][i]||v>host.bounds[1][i])
  ||!point(a.barycentric)||a.barycentric.some(w=>w<0||w>1)||Math.abs(a.barycentric.reduce((sum,w)=>sum+w,0)-1)>1e-12)return reject();
 cache.set(atlas,m);return m;
}
export const carotidSurfaceFeatureRegistry=carotidFeatureRegistry;

export function carotidSurfaceTermParts(atlas:Atlas):Map<string,Part[]>{
 const m=carotidFeatureRegistry(atlas),host=m&&atlas.parts.find(p=>p.id===m.hostPartId);
 return new Map(m&&host?sinusTerms.map(t=>[t,[host]]):[]);
}
export function selectedCarotidSurfaceFeature(atlas:Atlas,choice:HierarchyChoice|undefined){
 const m=carotidFeatureRegistry(atlas);
 if(!m||!choice||!oneHost(choice.elements)||choice.ta98Kind==='parent')return;
 const selectedTerm=choice.ta98Term??choice.terminology?.ta98;
 if(selectedTerm!==undefined&&!sinusTerms.some(t=>t===selectedTerm))return;
 return m.termGroups.find(g=>g.conceptId===choice.id||sinusTerms.some(t=>t===selectedTerm));
}
export function carotidSurfaceLabelReference(atlas:Atlas,choice:HierarchyChoice){
 const g=selectedCarotidSurfaceFeature(atlas,choice);
 return g?{partId:g.anchor.partId,point:g.anchor.point}:undefined;
}
/** A reference on an arterial carrier must not inherit “is an artery” from
 * that carrier in the details panel. Keep its qualification visible. */
export function carotidSurfaceChoiceDescription(choice:HierarchyChoice):string|undefined{
 const code=choice.ta98Term??choice.terminology?.ta98;
 return choice.id===featureId&&oneHost(choice.elements)&&choice.ta98Kind==='exact'
  &&sinusTerms.some(t=>t===code)&&choice.description===carotidSurfaceReferenceDescription
  ?carotidSurfaceReferenceDescription:undefined;
}
function referenceChoice(atlas:Atlas,g:CarotidSurfaceFeature,occurrence=g.term):HierarchyChoice{
 const c=atlas.concepts.find(c=>c.id===g.conceptId)!;
 return {...c,name:carotidSurfaceReferenceName,description:carotidSurfaceReferenceDescription,
  partOf:occurrence===term?c.partOf:hostId,
  ta98Term:occurrence,terminology:terminologyForConcept(c.id,occurrence,'exact')};
}
export function carotidReferenceChoices(atlas:Atlas):Map<string,CarotidReferenceChoice>{
 let found=references.get(atlas);if(found)return found;
 found=new Map();const m=carotidFeatureRegistry(atlas),host=m&&atlas.parts.find(p=>p.id===m.hostPartId);
 if(m&&host)for(const g of m.termGroups)for(const occurrence of sinusTerms)found.set(occurrence,{choice:referenceChoice(atlas,g,occurrence),host:{id:host.conceptId,name:host.name,elements:[host.id]}});
 references.set(atlas,found);return found;
}
export function carotidReferenceActive(reference:CarotidReferenceChoice,choice:HierarchyChoice|undefined){
 return !!choice&&choice.ta98Kind!=='parent'&&(choice.id===reference.choice.id||(choice.ta98Term??choice.terminology?.ta98)===reference.choice.ta98Term)&&oneHost(choice.elements);
}
/** Clearing selection ink retains the same existing anatomical carrier. */
export function toggleCarotidReferenceChoice(reference:CarotidReferenceChoice,choice:HierarchyChoice|undefined){
 return carotidReferenceActive(reference,choice)?reference.host:reference.choice;
}

/** Alias entries let Systems show a reference whose physical part keeps its
 * original ICA concept. The supplied hierarchy entries are never mutated. */
export function carotidSurfaceReferenceEntries(atlas:Atlas,entries:AnatomyEntry[]):AnatomyEntry[]{
 const m=carotidFeatureRegistry(atlas);if(!m)return [];
 const host=atlas.parts.find(p=>p.id===m.hostPartId)!,entry=entries.find(e=>e.id===host.conceptId&&e.parts.length===1&&e.parts[0]===host);
 if(!entry)return [];
 return m.termGroups.filter(g=>!entries.some(e=>e.id===g.conceptId)).map(g=>({...entry,id:g.conceptId,name:carotidSurfaceReferenceName,
  parts:[host],ta98Term:g.term,ta98Kind:'exact',taParentConcept:'ZA:Left common carotid artery',terminology:terminologyForConcept(g.conceptId,g.term,'exact')}));
}
