import {catalogueSourceReference} from "./source-reference";
import {nonComponentRelationship} from './anatomical-containment';
import index from './data/ta98-parents.json';
import {taExactEntityForConcept} from './anatomical-terminology';
import type {Concept,Part} from './anatomy';
import {displayLaterality} from './laterality';
import {hepaticConnectionRole,isHepaticNativeParent} from './hepatic-connection-controls';
import {isInstalledLiverPart} from './liver-controls';
import type {HierarchyChoice} from './hierarchy-choice';
import {isInstalledCarotidBifurcationPart,isInstalledLingualArterialPart,registeredLingualFamily} from './lingual-arterial-relations';
const parents=index.parents as Record<string,string>;
const exactTerm=(concept:Concept)=>concept.ta98Kind==='parent'?undefined:concept.ta98Term??taExactEntityForConcept(concept.id);
const sides=(concept:Concept,parts:Map<string,Part>)=>{const named=displayLaterality(concept.name).side;return named?new Set([named]):new Set(concept.elements.flatMap(id=>{const part=parts.get(id),side=part&&displayLaterality(part.name).side;return side?[side]:[];}));};
const identities=new WeakMap<Concept[],Map<string,Concept>>();
const conceptIndex=new WeakMap<Concept[],Map<string,Concept[]>>();
function conceptsByTerm(concepts:Concept[]):Map<string,Concept[]>{
 const cached=conceptIndex.get(concepts);if(cached)return cached;
 const byTerm=new Map<string,Concept[]>();
 for(const concept of concepts){
  const code=exactTerm(concept);if(!code)continue;
  const group=byTerm.get(code)??[];group.push(concept);byTerm.set(code,group);
 }
 conceptIndex.set(concepts,byTerm);return byTerm;
}

/** Resolve a direct TA98 parent to its stored whole-organ assembly outside
 * the TA98 browser. Do not invent a parent from a name suffix or spatial bounds,
 * skip an absent parent, or mistake a placement mapping for exact identity.
 */
export function ta98ModeledParent(choice:HierarchyChoice,concepts:Concept[],parts:Map<string,Part>):HierarchyChoice|undefined{
 if(choice.id==='selection-set'||choice.children)return;
 let byId=identities.get(concepts);if(!byId){byId=new Map(concepts.map(concept=>[concept.id,concept]));identities.set(concepts,byId);}
 const identity=choice.ta98Term||choice.ta98Kind?choice:byId.get(choice.id)??choice;
 const code=exactTerm(identity),parent=code&&parents[code];
 if(!parent)return;
 const members=choice.elements.map(id=>parts.get(id)).filter((part):part is Part=>!!part&&!part.suppressed);
 if(!members.length)return;
 const systems=new Set(members.map(part=>part.system)),sourceSides=sides(choice,parts);
 const candidates=(conceptsByTerm(concepts).get(parent)??[]).filter(concept=>concept.id!==choice.id).flatMap(concept=>{
  if(nonComponentRelationship(identity.name,concept.name))return [];
  const elements=concept.elements.filter(id=>{const part=parts.get(id);return part&&!part.suppressed;});
  if(!elements.length||!elements.some(id=>systems.has(parts.get(id)!.system)))return [];
  // A shared TA98 parent is not evidence that separately sourced organs fit.
  // Stored assemblies containing the child remain valid across source IDs.
  if(!members.every(part=>elements.includes(part.id))){
   const frame=nativePartFrame(members[0]);
   if(!frame||!members.every(part=>nativePartFrame(part)===frame)||!elements.every(id=>nativePartFrame(parts.get(id)!)===frame))return [];
  }
  const parentSides=sides({...concept,elements},parts);
  if(parentSides.size&&[...sourceSides].some(side=>!parentSides.has(side)))return [];
  return [{...concept,elements}];
 });
 // Prefer the actual stored assembly containing the selected region. Another
 // source organ with the same term must not displace that coordinated assembly.
 candidates.sort((a,b)=>Number(members.every(part=>b.elements.includes(part.id)))-Number(members.every(part=>a.elements.includes(part.id))));
 return candidates[0];
}

/** Known native catalogue frames only. HRA donors are independently registered
 * reference organs, so a shared HRA prefix does not establish a common frame. */
const nativeFrames=new WeakMap<Part,string|undefined>();
export function nativePartFrame(part:Part):string|undefined{
 if(nativeFrames.has(part))return nativeFrames.get(part);
 let frame=/^LOCAL:(?:male|female):/.exec(part.id)?.[0]??(part.id.startsWith('ZA:')?'ZA:':undefined);
 if(!frame&&(isInstalledLingualArterialPart(part)||isInstalledCarotidBifurcationPart(part)))frame='ZA:';
 // This original macro-body was constructed directly in the installed left
 // carotid frame. Recognize only its recorded tissue and packed layout; the
 // source qualification does not establish microscopic vascular attachments.
 if(!frame&&part.id==='TA98BODY:male:carotid-body:left:macro'&&part.conceptId===part.id
  &&part.sourceId==='Constructed left carotid chemoreceptor macro-body'&&part.name==='Carotid body (left)'
  &&part.system==='arterial'&&part.tissue==='organ'&&part.pairedTissue&&part.sectionSolid
  &&part.sectionAssembly==='carotid-body:male:left'&&part.vertexCount===242&&part.indexCount===1440
  &&part.positions===0&&part.normals===2904&&part.indices===4356&&part.normalSigns===10116
  &&part.geometryOrigin?.[0]===0&&part.geometryOrigin?.[1]===1.625&&part.geometryOrigin?.[2]===0
  &&!part.sourceOffset&&part.license==='CC BY-SA 4.0'
  &&catalogueSourceReference(part.provenance?.url)==='/ATTRIBUTION.md#open-anatomy-source')frame='ZA:';
 if(isInstalledLiverPart(part)&&part.tissueWallAssembly==='liver:female:lobes-v1')frame='LOCAL:female:';
 const cardiac=/^TA98HEART:(male|female):/.exec(part.id);
 if(!frame&&cardiac&&part.tissueWallAssembly===`heart:${cardiac[1]}`&&part.pairedTissue&&part.sectionSolid&&part.normalSigns!==undefined&&!part.sourceOffset
  &&((catalogueSourceReference(part.provenance?.url)==='/ATTRIBUTION.md#open-anatomy-source'&&part.license==='CC BY-SA 4.0')
    ||(part.provenance?.url==='https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html'&&part.license==='CC BY 4.0')))frame=cardiac[1]==='male'?'ZA:':'LOCAL:female:';
 // This particular BP3D head has an explicit four-papillary registration to
 // the native male frame. A common donor label alone never establishes it;
 // the provisional registration also does not certify a tissue attachment.
 const cardiacHead=part as Part&{reviewSource?:string;registration?:{method?:string;maxLandmarkResidualMeters?:number;pairs?:string[][]}};
 if(!frame&&part.id==='TA98CARDIAC:FJ2418'&&part.sourceId==='FJ2418'&&part.conceptId==='TA98RECOVERED:A12.1.04.007'
  &&part.system==='cardiac'&&!part.sourceOffset&&part.license==='CC BY 4.0'&&part.provenance?.url==='https://lifesciencedb.jp/bp3d/'&&cardiacHead.reviewSource==='male'
  &&cardiacHead.registration?.method==='Diagonal affine fit to min/max bounds of four corresponding named papillary muscles; anatomical registration provisional'
  &&Number.isFinite(cardiacHead.registration.maxLandmarkResidualMeters)&&cardiacHead.registration.maxLandmarkResidualMeters!>=0&&cardiacHead.registration.maxLandmarkResidualMeters!<=.0003
  &&JSON.stringify(cardiacHead.registration.pairs)===JSON.stringify([['FJ2419','ZA:Anterior papillary muscle of right ventricle'],['FJ2430','ZA:Inferior papillary muscle of right ventricle'],['FJ2437','ZA:Septal papillary muscle of right ventricle'],['FJ2429','ZA:Inferior papillary muscle of left ventricle']]))frame='ZA:';
 // Native source partitions retain their baked body coordinates. Require
 // matching source identity and provenance, rather than accepting any donor
 // with a TA98 term or a similar anatomical name. Regional source IDs append
 // their local region after a colon; that suffix does not create a new frame.
 const source=part.sourceId?.split(':')[0];
 const renal=part as Part&{renalSourcePartId?:string;renalRegion?:string};
 if(!frame&&source&&/^Kidney\.[lr]$/.test(source)&&part.id.startsWith(`TA98PART:ZA:${source}:`)
  &&renal.renalSourcePartId===`ZA:${source}`&&part.tissueWallAssembly===`renal:male:${source.endsWith('.l')?'left':'right'}`
  &&['capsule','renal-tissue'].includes(renal.renalRegion??'')&&part.pairedTissue&&part.sectionSolid&&part.normalSigns!==undefined
  &&!part.sourceOffset&&catalogueSourceReference(part.provenance?.url)==='/ATTRIBUTION.md#open-anatomy-source'&&part.license==='CC BY-NC 4.0')frame='ZA:';
 if(!frame&&/^VEINWALL:male:portal:(externa|media|intima)$/.test(part.id)&&part.system==='venous'
  &&source==='Hepatic portal vein'&&part.tissueWallAssembly==='portal:male:native-trunk'&&part.pairedTissue&&part.sectionSolid&&!part.sourceOffset
  &&catalogueSourceReference(part.provenance?.url)==='/ATTRIBUTION.md#open-anatomy-source'&&part.license==='CC BY-SA 4.0')frame='ZA:';
 if(!frame&&/^TA98WALL:male:thoracic-aorta:(externa|media|intima)$/.test(part.id)&&part.system==='arterial'
  &&source==='Thoracic aorta'&&part.tissueWallAssembly==='thoracic-aorta:male:native-trunk'&&part.pairedTissue&&part.sectionSolid
  &&part.normalSigns!==undefined&&part.uvs!==undefined&&!part.sourceOffset
  &&catalogueSourceReference(part.provenance?.url)==='/ATTRIBUTION.md#open-anatomy-source'&&part.license==='CC BY-SA 4.0')frame='ZA:';
 if(!frame&&source&&part.id.startsWith(`TA98PART:ZA:${source}:`)
  &&catalogueSourceReference(part.provenance?.url)==='/ATTRIBUTION.md#open-anatomy-source'
  &&part.license==='CC BY-SA 4.0')frame='ZA:';
 if(!frame&&(hepaticConnectionRole(part)||isHepaticNativeParent(part)))frame=(part as Part&{hepaticSourceModel?:string}).hepaticSourceModel==='local-female-ta98'?'LOCAL:female:':'ZA:';
 nativeFrames.set(part,frame);return frame;
}
const componentIndex=new WeakMap<Concept[],WeakMap<Map<string,Part>,Map<string,string[]>>>();
/** Invert the exact direct-parent resolver. Shared terminology alone cannot
 * bring an independently registered donor into a native organ assembly. */
export function ta98ModeledComponentIds(choice:HierarchyChoice,concepts:Concept[],parts:Map<string,Part>):string[]{
 if(choice.id==='selection-set'||choice.children)return [];
 let byParts=componentIndex.get(concepts);if(!byParts){byParts=new WeakMap();componentIndex.set(concepts,byParts);}
 let graph=byParts.get(parts);if(!graph){
  graph=new Map();
  for(const child of concepts){
   const parent=ta98ModeledParent(child,concepts,parts);if(!parent)continue;
   const parentIds=new Set(parent.elements),frames=new Set(parent.elements.flatMap(id=>{const part=parts.get(id),frame=part&&nativePartFrame(part);return frame?[frame]:[];}));
   const ids=child.elements.filter(id=>{const part=parts.get(id);if(!part||part.suppressed)return false;if(parentIds.has(id))return true;const frame=nativePartFrame(part);return !!frame&&frames.has(frame);});
   if(!ids.length)continue;const prior=graph.get(parent.id)??[];graph.set(parent.id,[...new Set([...prior,...ids])]);
  }
  byParts.set(parts,graph);
 }
 const ids=graph.get(choice.id)??[];
 // This finite parent is an owned trunk cell, not a duplicate outer envelope.
 // Preserve it when exposing its five child cells. Names or TA98 identity alone
 // must not assemble an unrelated donor into this reviewed native network.
 const trunk=parts.get(choice.id),side=/^TA98ARTERIAL:male:lingual:(left|right):finite-wall$/.exec(choice.id)?.[1] as 'left'|'right'|undefined;
 if(trunk&&side&&isInstalledLingualArterialPart(trunk)&&choice.elements.length===1&&choice.elements[0]===trunk.id){
  const family=registeredLingualFamily([...parts.values()],side),stored=concepts.find(c=>c.id===choice.id);
  const expected=family&&[...family.values()].filter(part=>part!==trunk);
  if(stored?.ta98Kind==='exact'&&stored.ta98Term==='A12.2.05.015'&&stored.partOf===`ZA:External carotid artery.${side[0]}`
   &&stored.elements.length===1&&stored.elements[0]===trunk.id&&expected?.length===5&&ids.length===5
   &&expected.every(part=>ids.includes(part.id)&&concepts.some(c=>c.id===part.conceptId&&c.partOf===stored.id&&c.ta98Kind==='exact'&&c.ta98Term===part.ta98Term&&c.elements.includes(part.id))))return [trunk.id,...ids];
 }
 return ids;
}
