import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,existsSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const base='.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';
const sourcePath='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json',bodyPath='.local-models/organ-models/spinal-cord/male-body.json';
const sha=b=>createHash('sha256').update(b).digest('hex'),expected='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de';
const sourceBytes=readFileSync(sourcePath),bodyBytes=readFileSync(bodyPath);assert.equal(sha(sourceBytes),expected);
const source=JSON.parse(sourceBytes),body=JSON.parse(bodyBytes),inputPath=`${base}/nearby-derivative-context-inventory-v24.json`,input=JSON.parse(readFileSync(inputPath));
assert.equal(input.source_catalogue_sha256,expected);
const ids=[...new Set(input.unregistered_or_unresolved_candidates.map(r=>r.part_id))].sort();assert.equal(ids.length,20);
const built=await build({stdin:{contents:`
 export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations.ts';
 export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities.ts';
 export {partVisible,sectionPartVisible,readViewUrl,resolveSelection} from './app/viewer-state.ts';
 export {isSkinPart} from './app/depth-layers.ts';
 export {nativePartFrame} from './app/ta98-modeled-parent.ts';
 export {DEFAULT_VISIBLE,SYSTEMS} from './app/anatomy.ts';
`,resolveDir:process.cwd(),sourcefile:'near-face-runtime-inventory.ts'},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const reviewed=api.withReviewedCatalogueIdentities(body),preferred=api.withPreferredSourceRepresentations(reviewed);
const standard=api.withHiddenReviewAlternatives(preferred,false),alternatives=api.withHiddenReviewAlternatives(preferred,true);
const byId=atlas=>new Map(atlas.parts.map(p=>[p.id,p])),sp=byId(source),bp=byId(body),rp=byId(reviewed),pp=byId(preferred),np=byId(standard),ap=byId(alternatives);
const baseState={focus:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:true,explode:0,visible:[...api.DEFAULT_VISIBLE,'integumentary'],selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0,skinOpacity:0};
const initial=api.readViewUrl('?model=local-ta98',standard,baseState),skinShown={...initial,skinOpacity:1};
const rows=ids.map(id=>{
 const original=sp.get(id),raw=bp.get(id),normal=np.get(id),inspection=ap.get(id);assert.ok(original&&raw&&normal&&inspection);
 assert.deepEqual({...raw,chunk:body.sourceChunkIndices[raw.chunk]},original,'Ordinary record must match current installed source exactly');
 const selected=api.resolveSelection(standard,[id]),normalSelection={...initial,selected,isolate:true,hidden:(initial.hidden??[]).filter(i=>!selected.includes(i))};
 const inspectSelection={...baseState,visible:api.SYSTEMS.map(s=>s.id),skinOpacity:1,selected:[id],isolate:true};
 const pref=pp.get(id),plainDefaultHidden=body.defaultHidden?.includes(id)??false;
 let exclusion='none';if(pref.suppressed&&!raw.suppressed)exclusion='specific-preferred-source-retirement';else if(normal.suppressed&&!pref.suppressed)exclusion='default-review-alternative-suppression';else if(normal.suppressed)exclusion='already-suppressed-source';
 return {id,name:normal.name,system:normal.system,raw_default_hidden:plainDefaultHidden,raw_part_default_hidden:raw.defaultHidden??false,
  reviewed_system:rp.get(id).system,preferred_suppressed:!!pref.suppressed,preferred_replacements:pref.replacementParts??[],normal_runtime_suppressed:!!normal.suppressed,
  normal_default_visible_skin0:api.partVisible(normal,initial,standard.parts),normal_visible_skin1:api.partVisible(normal,skinShown,standard.parts),
  normal_section_visible_skin1:api.sectionPartVisible(normal,skinShown),normal_explicit_selection_ids:selected,
  original_renders_when_explicitly_selected_in_normal_mode:api.partVisible(normal,normalSelection,standard.parts),
  inspection_mode_explicit_selection_renders_original:api.partVisible(inspection,inspectSelection,alternatives.parts),
  exclusion_basis:exclusion,is_skin_part:api.isSkinPart(normal),generic_nativePartFrame:api.nativePartFrame(normal)??null,
  generic_frame_identifier_is_not_registration_certificate:true,
  drawing_kind:normal.hairField?'instanced-metric-hair-shafts-not-carrier-triangles':'indexed-triangle-mesh',hairField:normal.hairField??null,
  declared_point_or_line_reference:!!normal.landmarkType||['point','line','curve'].includes(normal.referenceType),
  finite_tissue_metadata:!!normal.pairedTissue&&!!normal.sectionSolid,
  explicit_record_pose:{geometryOrigin:normal.geometryOrigin??[0,0,0],geometryScale:normal.geometryScale??1,sourceOffset:normal.sourceOffset??[0,0,0],geometryInstance:normal.geometryInstance??null},
  material_pose:{geometryOrigin:body.materials?.[normal.material]?.geometryOrigin??null,geometryInstance:body.materials?.[normal.material]?.geometryInstance??null},
  record_sha256:sha(Buffer.from(JSON.stringify(original))),sourceId:normal.sourceId??null,provenance:normal.provenance??null,reviewStatus:normal.reviewStatus??null,
  default_exclusion_is_viewer_rule_proof_not_anatomical_invalidity:true};
});
assert.equal(rows.filter(r=>r.normal_runtime_suppressed).length,16);
const ownershipFiles=['app/preferred-source-representations.ts','app/data/preferred-source-representations.json','app/page.tsx','app/local-ta98-layer-catalogue.ts','app/viewer-state.ts','app/depth-layers.ts','app/metric-hair-field.ts','app/organ-frame.ts','app/scene.tsx'];
const output={schema:'near-face-runtime-participation/v1',source_catalogue_sha256:expected,ordinary_body_sha256:sha(bodyBytes),
 source_inventory:inputPath,source_inventory_sha256:sha(readFileSync(inputPath)),registry_revision:input.registry_revision,
 code_bindings:Object.fromEntries(ownershipFiles.map(p=>[p,sha(readFileSync(p))])),
 pipeline:['ordinary male-body JSON','withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)','readViewUrl local default','partVisible/sectionPartVisible'],
 explicit_inspection_mode:'review=alternatives; source preferences still execute before alternatives are exposed',
 no_cardiac_label_bundle_needed_for_these_noncardiac_records:true,
 count:rows.length,default_excluded:rows.filter(r=>r.normal_runtime_suppressed).length,
 retained_runtime_ids:rows.filter(r=>!r.normal_runtime_suppressed).map(r=>r.id),normal_visible_skin0_ids:rows.filter(r=>r.normal_default_visible_skin0).map(r=>r.id),
 rows,no_catalogue_or_visibility_edit:true,no_collision_or_occupancy_conclusion:true};
const target=`${base}/nearby-unresolved-runtime-v25.json`;assert.ok(!existsSync(target));writeFileSync(target,JSON.stringify(output,null,2)+'\n');
assert.equal(sha(readFileSync(sourcePath)),expected);assert.equal(sha(readFileSync(bodyPath)),sha(bodyBytes));
console.log(JSON.stringify({output:target,count:output.count,default_excluded:output.default_excluded,retained_runtime_ids:output.retained_runtime_ids,normal_visible_skin0_ids:output.normal_visible_skin0_ids,rows:rows.map(r=>({id:r.id,exclusion:r.exclusion_basis,original_available_in_inspection:r.inspection_mode_explicit_selection_renders_original,drawing_kind:r.drawing_kind}))}));
