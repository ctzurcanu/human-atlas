import {catalogueSourceReference} from "./source-reference";
import type {Atlas,Concept,Part} from './anatomy';
import type {AnatomyEntry} from './anatomy-hierarchy';
import {terminologyForConcept} from './anatomical-terminology';

type Group={wholeConceptId:string;wholeTerm:string;remainderPartId:string;remainderConceptId:string;regionPartId:string;regionConceptId:string;regionTerm:string;side:'left'|'right';sourceId:string;sourceGeometrySha256:string;wholeName:string;sourceFaceCount:number;forearmFaceCount:number;remainderFaceCount:number;sharedInterfaces:{source_vertex_ids:number[];added_cap_triangles:number;same_bank_attributes:boolean}[]};
type Marker={revision:string;bankSha256:string;referenceOnly:boolean;lumenPatencyUnverified:boolean;anatomicalAcceptance:boolean;groups:Group[]};
type Registered={group:Group;whole:Concept;remainderConcept:Concept;regionConcept:Concept;remainder:Part;region:Part};
type SourcePart=Part&{nativeVenousRegionRevision?:string;nativeVenousSourcePartId?:string;nativeVenousSourceGeometrySha256?:string;nativeVenousRegionRole?:string;sourceOuterSurfaceOnly?:boolean;sourceLumenPatencyUnverified?:boolean;geometryInstance?:unknown};
const revision='male-native-forearm-source-regions-v1';
const pair=(a:string[],b:string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));

/** Exact native outer-source partitions only; no donor, wall or lumen claim. */
export function registeredNativeForearmVenousRegions(atlas:Atlas):Registered[]|undefined{
 const marker=(atlas as Atlas&{nativeVenousRegionIntegration?:Marker}).nativeVenousRegionIntegration;
 if(marker?.revision!==revision||!marker.referenceOnly||!marker.lumenPatencyUnverified||marker.anatomicalAcceptance!==false||marker.groups.length!==4||!/^[a-f0-9]{64}$/.test(marker.bankSha256))return;
 const parts=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c])),registered:Registered[]=[];
 for(const side of ['left','right'] as const)for(const kind of ['Cephalic','Basilic'] as const){
  const sourceId=`${kind} vein.${side[0]}`,wholeId=`ZA:${sourceId}`,term=kind==='Cephalic'?'A12.3.08.021':'A12.3.08.022',wholeTerm=kind==='Cephalic'?'A12.3.08.015':'A12.3.08.018';
  const rows=marker.groups.filter(g=>g.wholeConceptId===wholeId);if(rows.length!==1)return;
  const group=rows[0],regionId=`TA98PART:${wholeId}:${term}`,remainderId=`TA98SOURCE:${wholeId}:other-course`;
  if(group.side!==side||group.sourceId!==sourceId||group.wholeTerm!==wholeTerm||group.regionTerm!==term||group.regionPartId!==regionId||group.regionConceptId!==regionId||group.remainderPartId!==wholeId||group.remainderConceptId!==remainderId)return;
  const remainder=parts.get(wholeId) as SourcePart|undefined,region=parts.get(regionId) as SourcePart|undefined,whole=concepts.get(wholeId),remainderConcept=concepts.get(remainderId),regionConcept=concepts.get(regionId);
  if(!remainder||!region||!whole||!remainderConcept||!regionConcept)return;
  if(whole.ta98Kind!=='exact'||whole.ta98Term!==wholeTerm||!pair(whole.elements,[wholeId,regionId])||!pair(whole.componentElements??[],[wholeId,regionId]))return;
  if(regionConcept.ta98Kind!=='exact'||regionConcept.ta98Term!==term||regionConcept.partOf!==wholeId||!pair(regionConcept.elements,[regionId]))return;
  if(remainderConcept.ta98Term||remainderConcept.partOf!==wholeId||!pair(remainderConcept.elements,[wholeId]))return;
  if(remainder.conceptId!==wholeId||region.conceptId!==regionId||remainder.sourceId!==sourceId||region.sourceId!==`${sourceId}:forearm`||region.ta98Term!==term)return;
  const extra=(atlas.chunks[region.chunk] as typeof atlas.chunks[number]&{sha256?:string})?.sha256;
  if(extra!==marker.bankSha256||remainder.chunk!==region.chunk||remainder.vertexCount!==region.vertexCount||remainder.positions!==region.positions||remainder.normals!==region.normals||remainder.uvs!==region.uvs||remainder.indices===region.indices)return;
  if(group.sourceFaceCount!==group.forearmFaceCount+group.remainderFaceCount||region.indexCount!==group.forearmFaceCount*3||remainder.indexCount!==group.remainderFaceCount*3||group.sharedInterfaces.length!==2||group.sharedInterfaces.some(i=>i.source_vertex_ids.length!==12||i.added_cap_triangles!==0||!i.same_bank_attributes))return;
  for(const [part,role] of [[remainder,'remaining-source-course'],[region,'forearm-source-region']] as const){
   if(part.suppressed||part.system!=='venous'||part.sourceOffset||part.geometryOrigin||part.geometryInstance||part.geometryScale!==undefined||part.pairedTissue||part.sectionSolid||part.normalSigns!==undefined)return;
   if(part.nativeVenousRegionRevision!==revision||part.nativeVenousRegionRole!==role||part.nativeVenousSourcePartId!==wholeId||part.nativeVenousSourceGeometrySha256!==group.sourceGeometrySha256||!/^[a-f0-9]{64}$/.test(group.sourceGeometrySha256)||!part.sourceOuterSurfaceOnly||!part.sourceLumenPatencyUnverified)return;
   if(part.license!=='CC BY-SA 4.0'||catalogueSourceReference(part.provenance?.url)!=='/ATTRIBUTION.md#open-anatomy-source')return;
  }
  registered.push({group,whole,remainderConcept,regionConcept,remainder,region});
 }
 return registered;
}

/** Systems containment is the literal source assembly. TA98's canonical
 * parentage remains unchanged; its forearm terms are siblings of the wholes. */
export function withNativeForearmVenousEntries(atlas:Atlas,entries:AnatomyEntry[]):AnatomyEntry[]{
 const registered=registeredNativeForearmVenousRegions(atlas);if(!registered)return entries;
 const owned=new Map<string,Registered&{remainderEntry:boolean}>();
 for(const record of registered){owned.set(record.remainder.id,{...record,remainderEntry:true});owned.set(record.region.id,{...record,remainderEntry:false});}
 const selected=entries.filter(e=>e.parts.some(p=>owned.has(p.id)));
 if(selected.length!==8||selected.some(e=>e.parts.length!==1)||new Set(selected.flatMap(e=>e.parts.map(p=>p.id))).size!==8)return entries;
 return entries.map(entry=>{
  const record=owned.get(entry.parts[0]?.id);if(!record)return entry;
  const wholeName=record.group.wholeName.replace(/ \((left|right)\)$/,'');
  const concept=record.remainderEntry?record.remainderConcept:record.regionConcept;
  return {...entry,id:concept.id,name:concept.name,ta98Term:concept.ta98Term,ta98Kind:concept.ta98Kind,
   terminology:terminologyForConcept(concept.id,concept.ta98Term,concept.ta98Kind),taParentConcept:record.whole.id,
   assemblySystemPath:['Veins','Upper limbs','Superficial veins of upper limb',wholeName],
   ancestryTerms:{...entry.ancestryTerms,[wholeName]:terminologyForConcept(record.whole.id,record.group.wholeTerm,'exact')}};
 });
}
