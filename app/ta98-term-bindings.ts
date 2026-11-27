import type {Atlas,Part} from './anatomy';
import ta98Parents from './data/ta98-parents.json';
import {collectionReviewMetadata} from './source-reference';

/** Explicit tissue, source-surface or anatomical source-identity controls.
 * Identity controls do not assert a finite wall, patency or anatomical fit. */
export interface Ta98SourceCollectionIdentity {
 id:string;term:string;sourceId:string;sourceAnchorConceptId:string;identityOnly:true;
 taParentTerm?:string;taParentConceptId?:string;
 members:{partId:string;sourceId:string;conceptId:string;ta98Term?:string}[];
}
/** A singular, reviewed course selects its canonical source carrier; it adds no tissue mesh. */
export interface Ta98SourceCourseIdentity {
 schema:'source-course-identity/v1';id:string;term:string;taParentTerm:string;taParentConceptId:string;sourceClinicalParentConceptId:string;
 identityOnly:true;dimension:'source-surface-envelope';finiteWallEstablished:false;lumenPatencyEstablished:false;anatomicalAcceptance:false;
 sourcePartId:string;sourceId:string;sourceBank:{url:string;bytes:number;sha256:string};
 logicalCarrierPose:Record<string,unknown>;evidence:Record<string,unknown>;carrier:Record<string,unknown>;
}
export interface Ta98TermBinding {term:string;conceptId:string;partIds:string[];kind?:'exact'|'parent';representation?:'finite-tissue'|'source-surface-reference'|'source-identity'|'source-parent-placement'|'source-collection-identity'|'source-course-identity';sourceCollectionId?:string;sourceCourseId?:string}

// Deliberately scoped to reviewed source courses, not a general name-based alias mechanism.
const REVIEWED_SOURCE_COURSE:Ta98SourceCourseIdentity={"schema":"source-course-identity/v1","id":"TA98COURSE:male:pulmonary-right-posterior-ascending","term":"A12.2.01.108","taParentTerm":"A12.2.01.107","taParentConceptId":"BDM:male:human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","sourceClinicalParentConceptId":"BDM:male:human_20_male_cardiovascular_system-branches_of_right_pulmonary_artery_ID","identityOnly":true,"dimension":"source-surface-envelope","finiteWallEstablished":false,"lumenPatencyEstablished":false,"anatomicalAcceptance":false,"sourcePartId":"BDM:male:human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","sourceId":"human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","sourceBank":{"url":"/local-models/geometry-workbench/biodigital-male-circulation-v1/atlas-cardiac-blood-v2/geometry.bin","bytes":36619168,"sha256":"b456e897d11d7e5a91958f3b6e85b93b6b948503b3fb606c77dc136583eabfca"},"logicalCarrierPose":{"geometryOrigin":[-0.05043221126410182,1.4672558650182153,-0.029527814715413454],"geometryScale":null,"sourceOffset":null,"geometryInstance":null},"evidence":{"crosswalkSha256":"7d19d5b27c469343382d198d4726b4008f3e95abc181295f64bd31055ba453ac","sourceCatalogueSha256":"caad8116d7dec25273ee64c89f654afb19bf2ed48e507409afc38d38eaf9ba39","originalAssetSha256":"6f54cbd1ed6711d37e7184b7f34b282e0ad966410237d43c1ce1d711e3cb47a3","bronchialSourceId":"human_20_male_respiratory_system-right_posterior_segmental_bronchus_ID","bronchialRawBdaSha256":"30508438a4db8470a4e1cfaa6b79bf56ad8ea41202661be9d489e69136b09b9e","feedingSourceId":"human_20_male_cardiovascular_system-right_interlobar_artery_ID","sourceFaces":1860,"matchedFeedingMouth":1},"carrier":{"id":"BDM:male:human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","name":"Right posterior segmental artery","conceptId":"BDM:male:human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","sourceId":"human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","sourceGroups":["BDM:male:human_20_male_cardiovascular_system-root","BDM:male:human_20_male_cardiovascular_system-cardiovascular_system_ID","BDM:male:human_20_male_cardiovascular_system-arterial_system_ID","BDM:male:human_20_male_cardiovascular_system-arteries_of_thorax_ID","BDM:male:human_20_male_cardiovascular_system-pulmonary_arteries_ID","BDM:male:human_20_male_cardiovascular_system-branches_of_right_pulmonary_artery_ID"],"system":"arterial","positions":5576308,"normals":5588068,"normalSigns":5593948,"indices":5594928,"vertexCount":980,"indexCount":5580,"geometryOrigin":[-0.05043221126410182,1.4672558650182153,-0.029527814715413454],"bounds":[[-0.08552743351722213,1.4285702220922885,-0.05298101369589896],[-0.01533698901098152,1.505941507944142,-0.00607461573492795]],"material":"BDM:review:vein_pbr_mtl:arterial","license":"Source redistribution rights unverified; local review only","provenance":{"label":"Registered collection, original named vascular structure","url":"/ATTRIBUTION.md#registered-collection-source","detail":"Original hierarchy and source arrays retained. One common regional review field; no accepted TA98 binding. Source Heart label is organizational provenance."},"ta98Term":"A12.2.01.107","coverageLimitation":"Named original right pulmonary source vessel; original branches-of-right-pulmonary-artery ancestry retained. This corrects a hepatic-name homonym. Fine branch completion, lumen continuity and body fit remain unaccepted."}};
const REVIEWED_SOURCE_COURSE_CONCEPT={"id":"TA98COURSE:male:pulmonary-right-posterior-ascending","name":"Ascending branch of posterior segmental artery (right)","elements":["BDM:male:human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_posterior_segmental_artery_ID","ta98Term":"A12.2.01.108","ta98Kind":"exact","sourceCourseId":"TA98COURSE:male:pulmonary-right-posterior-ascending","coverageLimitation":"Complete original BioDigital posterior ascending arterial source course to the right upper-lobe posterior segment. This selects the existing whole source envelope; body fit, finite vessel walls, perfusion and patent recipient joins remain unaccepted."};
const REVIEWED_PULMONARY_VENOUS_COURSES:{course:Ta98SourceCourseIdentity;concept:Record<string,unknown>;parentConcept:Record<string,unknown>;sourceConcept:Record<string,unknown>}[]=[{"course":{"schema":"source-course-identity/v1","id":"TA98COURSE:male:pulmonary-right-middle-vein-lateral","term":"A12.3.02.112","taParentTerm":"A12.3.02.111","taParentConceptId":"BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","sourceClinicalParentConceptId":"BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID","identityOnly":true,"dimension":"source-surface-envelope","finiteWallEstablished":false,"lumenPatencyEstablished":false,"anatomicalAcceptance":false,"sourcePartId":"BDM:male:human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID","sourceId":"human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID","sourceBank":{"url":"/local-models/geometry-workbench/biodigital-male-circulation-v1/atlas-cardiac-blood-v2/geometry.bin","bytes":36619168,"sha256":"b456e897d11d7e5a91958f3b6e85b93b6b948503b3fb606c77dc136583eabfca"},"logicalCarrierPose":{"geometryOrigin":[-0.06586144255227175,1.4181025773689075,0.030140180592065344],"geometryScale":null,"sourceOffset":null,"geometryInstance":null},"evidence":{"crosswalkSha256":"faf5fd30391613e08cb510e7923aa79ef16e7f4f5c122d1ba11289272a47e30d","sourceCatalogueSha256":"caad8116d7dec25273ee64c89f654afb19bf2ed48e507409afc38d38eaf9ba39","originalAssetSha256":"0212cc07d208c22a10bef3c93c1df2a3bcc01b8e24303a1d00caa966d036f26e","bronchialSourceId":"human_20_male_respiratory_system-right_lateral_segmental_bronchus_ID","bronchialRawBdaSha256":"b435d26fe22cc054eec7fd82f8d4cd909b4c6753f82613dfcde5f70164e3d995","drainingSourceId":"human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","sourceFaces":3288,"proximalSourceMouthGapSourceUnits":0.4104411292677784,"matchedDrainingMouth":1,"distalMouths":38},"carrier":{"id":"BDM:male:human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID","name":"Right lateral segmental vein","conceptId":"BDM:male:human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID","sourceId":"human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID","sourceGroups":["BDM:male:human_20_male_cardiovascular_system-root","BDM:male:human_20_male_cardiovascular_system-cardiovascular_system_ID","BDM:male:human_20_male_cardiovascular_system-venous_system_ID","BDM:male:human_20_male_cardiovascular_system-veins_of_thorax_ID","BDM:male:human_20_male_cardiovascular_system-pulmonary_veins_ID","BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID"],"system":"venous","positions":22780984,"normals":22801672,"normalSigns":22812016,"indices":22813740,"vertexCount":1724,"indexCount":9864,"geometryOrigin":[-0.06586144255227175,1.4181025773689075,0.030140180592065344],"bounds":[[-0.10262049960721642,1.3811841584639406,-0.004029958813277128],[-0.029102385497327066,1.4550209962738747,0.06431031999740781]],"material":"BDM:review:artery_pbr_mtl:venous","license":"Source redistribution rights unverified; local review only","provenance":{"label":"Registered collection, original named vascular structure","url":"/ATTRIBUTION.md#registered-collection-source","detail":"Original hierarchy and source arrays retained. One common regional review field; no accepted TA98 binding. Source Heart label is organizational provenance."}}},"concept":{"id":"TA98COURSE:male:pulmonary-right-middle-vein-lateral","name":"Lateral part of middle lobe vein (right)","elements":["BDM:male:human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","ta98Term":"A12.3.02.112","ta98Kind":"exact","sourceCourseId":"TA98COURSE:male:pulmonary-right-middle-vein-lateral","coverageLimitation":"Complete original BioDigital named lateral segmental venous source course of the right middle lobe. This selects the existing whole source envelope. The common middle-lobe stem is separate; proximal source mouth gap, body fit, finite vessel walls, drainage completeness and patent recipient joins remain unaccepted."},"parentConcept":{"id":"BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","name":"Right middle lobe vein","elements":["BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID","ta98Term":"A12.3.02.111","ta98Kind":"exact","coverageLimitation":"Source vessel identity established; fit and lumen continuity remain under review."},"sourceConcept":{"id":"BDM:male:human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID","name":"Right lateral segmental vein","elements":["BDM:male:human_20_male_cardiovascular_system-right_lateral_segmental_vein_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID"}},{"course":{"schema":"source-course-identity/v1","id":"TA98COURSE:male:pulmonary-right-middle-vein-medial","term":"A12.3.02.113","taParentTerm":"A12.3.02.111","taParentConceptId":"BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","sourceClinicalParentConceptId":"BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID","identityOnly":true,"dimension":"source-surface-envelope","finiteWallEstablished":false,"lumenPatencyEstablished":false,"anatomicalAcceptance":false,"sourcePartId":"BDM:male:human_20_male_cardiovascular_system-right_medial_segmental_vein_ID","sourceId":"human_20_male_cardiovascular_system-right_medial_segmental_vein_ID","sourceBank":{"url":"/local-models/geometry-workbench/biodigital-male-circulation-v1/atlas-cardiac-blood-v2/geometry.bin","bytes":36619168,"sha256":"b456e897d11d7e5a91958f3b6e85b93b6b948503b3fb606c77dc136583eabfca"},"logicalCarrierPose":{"geometryOrigin":[-0.06942034385074734,1.3734110201355534,0.04452415186853877],"geometryScale":null,"sourceOffset":null,"geometryInstance":null},"evidence":{"crosswalkSha256":"faf5fd30391613e08cb510e7923aa79ef16e7f4f5c122d1ba11289272a47e30d","sourceCatalogueSha256":"caad8116d7dec25273ee64c89f654afb19bf2ed48e507409afc38d38eaf9ba39","originalAssetSha256":"82ed5d3ac38627656ca84f043a50162c9caa5b50ed6acf4e9cd7b73a4f223732","bronchialSourceId":"human_20_male_respiratory_system-right_medial_segmental_bronchus_ID","bronchialRawBdaSha256":"324cb5a2f368a6450a2755063c21f1faf7f71c97ae73ad8383c4b94c29493370","drainingSourceId":"human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","sourceFaces":2520,"proximalSourceMouthGapSourceUnits":0.4104205601325894,"matchedDrainingMouth":1,"distalMouths":29},"carrier":{"id":"BDM:male:human_20_male_cardiovascular_system-right_medial_segmental_vein_ID","name":"Right medial segmental vein","conceptId":"BDM:male:human_20_male_cardiovascular_system-right_medial_segmental_vein_ID","sourceId":"human_20_male_cardiovascular_system-right_medial_segmental_vein_ID","sourceGroups":["BDM:male:human_20_male_cardiovascular_system-root","BDM:male:human_20_male_cardiovascular_system-cardiovascular_system_ID","BDM:male:human_20_male_cardiovascular_system-venous_system_ID","BDM:male:human_20_male_cardiovascular_system-veins_of_thorax_ID","BDM:male:human_20_male_cardiovascular_system-pulmonary_veins_ID","BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID"],"system":"venous","positions":22853196,"normals":22869060,"normalSigns":22876992,"indices":22878316,"vertexCount":1322,"indexCount":7560,"geometryOrigin":[-0.06942034385074734,1.3734110201355534,0.04452415186853877],"bounds":[[-0.1137074998768541,1.33812865654831,0.007975037999984648],[-0.025133187824640576,1.4086933837227968,0.08107326573709288]],"material":"BDM:review:artery_pbr_mtl:venous","license":"Source redistribution rights unverified; local review only","provenance":{"label":"Registered collection, original named vascular structure","url":"/ATTRIBUTION.md#registered-collection-source","detail":"Original hierarchy and source arrays retained. One common regional review field; no accepted TA98 binding. Source Heart label is organizational provenance."}}},"concept":{"id":"TA98COURSE:male:pulmonary-right-middle-vein-medial","name":"Medial part of middle lobe vein (right)","elements":["BDM:male:human_20_male_cardiovascular_system-right_medial_segmental_vein_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","ta98Term":"A12.3.02.113","ta98Kind":"exact","sourceCourseId":"TA98COURSE:male:pulmonary-right-middle-vein-medial","coverageLimitation":"Complete original BioDigital named medial segmental venous source course of the right middle lobe. This selects the existing whole source envelope. The common middle-lobe stem is separate; proximal source mouth gap, body fit, finite vessel walls, drainage completeness and patent recipient joins remain unaccepted."},"parentConcept":{"id":"BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID","name":"Right middle lobe vein","elements":["BDM:male:human_20_male_cardiovascular_system-right_middle_lobe_vein_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID","ta98Term":"A12.3.02.111","ta98Kind":"exact","coverageLimitation":"Source vessel identity established; fit and lumen continuity remain under review."},"sourceConcept":{"id":"BDM:male:human_20_male_cardiovascular_system-right_medial_segmental_vein_ID","name":"Right medial segmental vein","elements":["BDM:male:human_20_male_cardiovascular_system-right_medial_segmental_vein_ID"],"partOf":"BDM:male:human_20_male_cardiovascular_system-right_pulmonary_veins_ID"}}];
const canonical=(value:unknown):string=>JSON.stringify(value,(_key,item)=>item&&typeof item==='object'&&!Array.isArray(item)?Object.fromEntries(Object.keys(item).sort().map(key=>[key,item[key]])):item);
function sourceCoursePart(atlas:Atlas,binding:Ta98TermBinding,parts:Map<string,Part>,concepts:Map<string,Atlas['concepts'][number]>,hidden:Set<string>):Part|null{
 const reviewed=REVIEWED_PULMONARY_VENOUS_COURSES.find(r=>r.course.id===binding.sourceCourseId);
 const expected=reviewed?.course??REVIEWED_SOURCE_COURSE;
 const expectedConcept=reviewed?.concept??REVIEWED_SOURCE_COURSE_CONCEPT;
 const manifests=atlas.sourceCourseIdentities?.filter(m=>m?.id===binding.sourceCourseId)??[];
 if(binding.sourceCourseId!==expected.id||binding.term!==expected.term||binding.kind!=='exact'||binding.conceptId!==expected.id||binding.partIds.length!==1||binding.partIds[0]!==expected.sourcePartId||manifests.length!==1||canonical(manifests[0])!==canonical(expected))return null;
 const profile=collectionReviewMetadata(atlas as Atlas&{registeredCollectionReview?:Record<string,unknown>});
 if(atlas.sex!=='male'||profile?.optInOnly!==true||profile.installedAsDefault!==false||profile.anatomicalAcceptance!==false||profile.sourceBaseCatalogueSha256!=='6682b34ac0f314cddb919b6c385e1a4ea7c3deaa0e68177d019dc07b217af7bd'||profile.sourceMasterSha256!=='46262d8ab07682fcae3ff0ee634f652abad766c1e010f6e2ec311755aea9abf5')return null;
 const p=parts.get(expected.sourcePartId),c=concepts.get(expected.id),parent=concepts.get(expected.taParentConceptId);
 if(!p||p.suppressed||p.defaultHidden||hidden.has(p.id)||!c||canonical(c)!==canonical(expectedConcept)||parent?.ta98Term!==expected.taParentTerm||parent.ta98Kind!=='exact'||parent.partOf!==expected.sourceClinicalParentConceptId||parent.elements.length!==1||parent.elements[0]!==expected.taParentConceptId||!concepts.has(expected.sourceClinicalParentConceptId))return null;
 if(reviewed&&(canonical(parent)!==canonical(reviewed.parentConcept)||canonical(concepts.get(p.conceptId))!==canonical(reviewed.sourceConcept)))return null;
 const {chunk:_chunk,...record}=p;
 if(canonical(record)!==canonical(expected.carrier)||canonical(atlas.chunks[p.chunk])!==canonical(expected.sourceBank))return null;
 return p;
}

export function catalogueTa98TermParts(atlas:Atlas):Map<string,Part[]>{
 const result=new Map<string,Part[]>();
 if(!Array.isArray(atlas.ta98TermBindings))return result;
 const parts=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const hidden=new Set(atlas.defaultHidden??[]);
 // An Eye component can retain its independent, explicitly recorded vascular
 // provenance. Source-collection validation must follow that source parent.
 type OcularParents={eyeOcularCirculationMembership?:{groups:{eyeConceptId:string;sourceParents:{conceptId:string;sourceParentConceptId:string}[]}[]};
  eyeCanonicalConceptMembership?:{groups:{eyeConceptId:string;eyeballConceptId:string;activeComponentIds:string[]}[]}};
 const canonicalEyes=(atlas as Atlas&OcularParents).eyeCanonicalConceptMembership?.groups??[];
 const ocularParents=new Map(((atlas as Atlas&OcularParents).eyeOcularCirculationMembership?.groups??[])
  .flatMap(group=>group.sourceParents.map(parent=>[parent.conceptId,{eye:group.eyeConceptId,source:parent.sourceParentConceptId}] as const)));
 const sourceParentId=(concept:Atlas['concepts'][number])=>{
  const declared=ocularParents.get(concept.id),original=(concept as typeof concept&{sourceParentConceptId?:string}).sourceParentConceptId;
  const canonical=declared&&canonicalEyes.find(group=>group.eyeballConceptId===declared.eye&&group.eyeConceptId===concept.partOf&&
   concepts.get(declared.eye)?.partOf===group.eyeConceptId&&concept.elements.every(id=>group.activeComponentIds.includes(id)));
  const parent=concepts.get(concept.partOf??'');
  return declared&&(concept.partOf===declared.eye||canonical)&&original===declared.source&&
   concept.elements.every(id=>parent?.elements.includes(id))?declared.source:concept.partOf;
 };
 for(const binding of atlas.ta98TermBindings){
  if(!binding||typeof binding.term!=='string'||!/^A\d{2}\.\d\.\d{2}\.\d{3}[FM]?$/.test(binding.term)||typeof binding.conceptId!=='string'||!Array.isArray(binding.partIds)||!binding.partIds.length||binding.partIds.some(id=>typeof id!=='string')||new Set(binding.partIds).size!==binding.partIds.length)continue;
  const c=concepts.get(binding.conceptId);
  const kind=binding.kind??'exact';
  if(!['exact','parent'].includes(kind)||c?.ta98Kind!==kind||c.ta98Term!==binding.term||c.elements.length!==binding.partIds.length||!c.elements.every(id=>binding.partIds.includes(id)))continue;
  if(kind==='parent'&&!c.coverageLimitation?.trim())continue;
  const sourceReference=binding.representation==='source-surface-reference';
  const sourceIdentity=binding.representation==='source-identity';
  const sourceParent=binding.representation==='source-parent-placement';
  const sourceCollection=binding.representation==='source-collection-identity';
  const sourceCourse=binding.representation==='source-course-identity';
  if(binding.representation&&!['finite-tissue','source-surface-reference','source-identity','source-parent-placement','source-collection-identity','source-course-identity'].includes(binding.representation))continue;
  if(sourceCourse){
   const p=sourceCoursePart(atlas,binding,parts,concepts,hidden);
   if(p)result.set(binding.term,[...new Map([...(result.get(binding.term)??[]),p].map(p=>[p.id,p])).values()]);
   continue;
  }
  if((sourceReference||sourceIdentity||sourceParent||sourceCollection)&&!c.coverageLimitation?.trim())continue;
  if((sourceIdentity||sourceCollection)&&kind!=='exact')continue;
  // A named source subdivision can be placed under a broader TA98 term
  // without claiming that it represents that whole term or its exact identity.
  if(sourceParent&&(kind!=='parent'||binding.partIds.length!==1))continue;
  const members=binding.partIds.map(id=>parts.get(id));
  if(sourceCollection){
   const manifests=atlas.sourceCollectionIdentities?.filter(m=>m.id===binding.sourceCollectionId)??[];
   if(!binding.sourceCollectionId||c.sourceCollectionId!==binding.sourceCollectionId||manifests.length!==1)continue;
   const manifest=manifests[0];
   if(manifest.identityOnly!==true||manifest.term!==binding.term||!manifest.sourceId?.trim()||!manifest.sourceAnchorConceptId?.trim()||!concepts.has(manifest.sourceAnchorConceptId)||!Array.isArray(manifest.members)||manifest.members.length!==members.length)continue;
   // Source organizers establish provenance; a reviewed TA98 parent may
   // establish anatomical containment without rewriting that provenance.
   const distinctParent=manifest.taParentTerm!==undefined||manifest.taParentConceptId!==undefined;
   if(distinctParent){
    const parent=concepts.get(manifest.taParentConceptId??''),term=manifest.taParentTerm;
    const expected=(ta98Parents.parents as Record<string,string>)[binding.term];
    if(!term||expected!==term||!parent||parent.id===c.id||c.partOf!==parent.id||parent.ta98Kind!=='exact'||parent.ta98Term!==term)continue;
    const system=binding.term.startsWith('A12.2.')?'arterial':binding.term.startsWith('A12.3.')?'venous':undefined;
    const parentParts=parent.elements.map(id=>parts.get(id));
    if(!system||!parentParts.length||parentParts.some(p=>!p||p.suppressed||p.defaultHidden||hidden.has(p.id)||p.system!==system)||!parentParts.some(p=>p?.ta98Term===term)||members.some(p=>p?.system!==system))continue;
    const seen=new Set([c.id]);let ancestor:typeof parent|undefined=parent;
    while(ancestor&&!seen.has(ancestor.id)){seen.add(ancestor.id);const sourceParent=sourceParentId(ancestor);ancestor=sourceParent?concepts.get(sourceParent):undefined;}
    if(ancestor)continue;
   }else if(sourceParentId(c)!==manifest.sourceAnchorConceptId)continue;
   if(manifest.members.some(m=>!m||typeof m.partId!=='string'||typeof m.sourceId!=='string'||!m.sourceId.trim()||typeof m.conceptId!=='string'||m.ta98Term!==undefined&&typeof m.ta98Term!=='string')||new Set(manifest.members.map(m=>m.partId)).size!==members.length||new Set(manifest.members.map(m=>m.sourceId)).size!==members.length)continue;
   const declared=new Map(manifest.members.map(m=>[m.partId,m]));
   if(members.some(p=>{
    if(!p||p.suppressed||p.defaultHidden||hidden.has(p.id)||!['arterial','venous'].includes(p.system)||p.conceptId===c.id)return true;
    const member=declared.get(p.id),leaf=concepts.get(p.conceptId);
    if(!member||p.sourceId!==member.sourceId||p.conceptId!==member.conceptId||p.ta98Term!==member.ta98Term||!leaf||leaf.ta98Term!==p.ta98Term||!leaf.elements.includes(p.id)||!p.sourceGroups?.includes(manifest.sourceAnchorConceptId))return true;
    return !Number.isInteger(p.vertexCount)||p.vertexCount<=0||!Number.isInteger(p.indexCount)||p.indexCount<=0||p.indexCount%3!==0||!Number.isInteger(p.chunk)||!atlas.chunks[p.chunk]||[p.positions,p.normals,p.indices].some(offset=>!Number.isInteger(offset)||offset<0);
   }))continue;
   result.set(binding.term,[...new Map([...(result.get(binding.term)??[]),...members as Part[]].map(p=>[p.id,p])).values()]);
   continue;
  }
  if(members.some(p=>{
   if(!p||p.suppressed||p.defaultHidden||hidden.has(p.id)||p.conceptId!==c.id||p.ta98Term!==binding.term)return true;
   if(sourceParent)return !p.sourceId?.trim()||!['arterial','venous'].includes(p.system)
    ||!Number.isInteger(p.vertexCount)||p.vertexCount<=0||!Number.isInteger(p.indexCount)||p.indexCount<=0||p.indexCount%3!==0
    ||!Number.isInteger(p.chunk)||!atlas.chunks[p.chunk]||[p.positions,p.normals,p.indices].some(offset=>!Number.isInteger(offset)||offset<0);
   if(sourceIdentity)return false;
   if(!sourceReference)return !p.sectionSolid;
   const reference=p as Part&{sourceOuterSurfaceOnly?:boolean;sourceLumenPatencyUnverified?:boolean;sourceOwner?:string};
   const owner=reference.sourceOwner?parts.get(reference.sourceOwner):undefined;
   const fields=['geometryOrigin','sourceOffset','geometryInstance','normalSigns'];
   const record=p as unknown as Record<string,unknown>,source=owner as unknown as Record<string,unknown>|undefined;
   return !reference.sourceOuterSurfaceOnly||!reference.sourceLumenPatencyUnverified||!owner||owner===p||owner.suppressed||
    fields.some(key=>JSON.stringify(record[key])!==JSON.stringify(source?.[key]))||
    p.chunk!==owner.chunk||p.positions!==owner.positions||p.normals!==owner.normals||p.uvs!==owner.uvs||p.vertexCount!==owner.vertexCount||
    p.indexCount<=0||p.indexCount%3!==0||p.indices!==owner.indices+owner.indexCount*4;
  }))continue;
  result.set(binding.term,[...new Map([...(result.get(binding.term)??[]),...members as Part[]].map(p=>[p.id,p])).values()]);
 }
 return result;
}
