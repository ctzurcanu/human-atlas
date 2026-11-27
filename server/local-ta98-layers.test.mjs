import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/viewer-state';export * from './app/view-part-codes';export * from './app/ta98-term-bindings';export * from './app/anatomical-components';export * from './app/reviewed-catalogue-identities';export * from './app/tissue-interface-visibility';export * from './app/pancreatic-venous-relations';export * from './app/anatomy-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {localTa98LayerCatalogue,preserveLocalTa98LayerCatalogueQuery,readViewUrl,partFingerprint,packPartIds,compatiblePartAtlas,catalogueTa98TermParts,anatomicalComponentIds,withReviewedCatalogueIdentities,tissueNeighborOccludes,tissuePosesCoincide,pancreaticVenousRelations,reviewedPancreaticVenousPair,hierarchyEntries,systemPathFor}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const read=p=>JSON.parse(fs.readFileSync(p));
test('Fixed-context vein delivery preserves surrounding anatomy and its source movement bounds',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/',before=read(base+'blood-installed-v1334/circulation-atlas.json'),after=read(base+'blood-installed-v1339/circulation-atlas.json'),r=after.bloodSurfaceRepairs;
 assert.equal(r.revision,'v1339');assert.equal(r.changedParts,3);assert.equal(r.removedIntrinsicPairs,35);assert.equal(r.remainingSelfCrossingParts,108);assert.equal(r.skeletonChanged,false);assert.equal(r.facesRemoved,0);
 const changed=new Set(r.rows.map(p=>p.partId)),old=new Map(before.parts.map(p=>[p.id,p]));assert.ok(changed.has('BDM:male:human_20_male_cardiovascular_system-left_deep_circumflex_iliac_vein_ID'));assert.ok(changed.has('BDM:male:human_20_male_cardiovascular_system-right_deep_circumflex_iliac_vein_ID'));
 for(const p of after.parts){const prior=old.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);assert.equal(p.ta98Term,prior.ta98Term);}}
 for(const key of ['genitalCourseRepair','genitalProgress','genitalCurrentVascularFit','handPalmDelivery','aorticContinuity','ta98TermBindings','sourceCollectionIdentities','concepts','materials','defaultHidden'])assert.deepEqual(after[key],before[key]);
 const scope=read(base+'blood-installed-v1339/packed-context-check.json');assert.equal(scope.passed,true);assert.equal(scope.rows.length,3);assert.equal(scope.anatomicalAcceptance,false);assert.ok(scope.rows.every(c=>c.contextTrianglesLiteral&&c.geometryChangesInsideProvedMovementDomain&&c.newPrimaryContextPairs===0&&c.intrinsicPairsAfter===0));
 const packed=read(base+'blood-installed-v1339/packed-neighbor-check.json');assert.equal(packed.passedParts,3);assert.ok(read(packed.report).checks.every(c=>c.passed&&Object.values(c.contacts).every(s=>s.newCrossingFacePairs===0)));assert.ok(read(base+'blood-installed-v1339/topology-retention.json').rows.every(c=>c.same));
});
test('Radicular vein repair uses the paired source inventory while preserving the artery and fixed anatomy',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/',before=read(base+'blood-installed-v1328/circulation-atlas.json'),after=read(base+'blood-installed-v1334/circulation-atlas.json'),r=after.bloodSurfaceRepairs;
 assert.equal(r.revision,'v1334');assert.equal(r.changedParts,1);assert.equal(r.removedIntrinsicPairs,12);assert.equal(r.remainingSelfCrossingParts,111);assert.equal(r.skeletonChanged,false);assert.equal(r.facesRemoved,0);
 const row=r.rows[0];assert.equal(row.partId,'BDM:male:human_20_male_cardiovascular_system-right_radicular_veins_ID');assert.equal(row.repair.targetPairCount,12);assert.equal(row.repair.remainingTargetPairs,0);assert.equal(row.repair.allSourceDefectsRemoved,false);
 const old=new Map(before.parts.map(p=>[p.id,p]));for(const p of after.parts){const prior=old.get(p.id);if(p.id!==row.partId)assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);}}
 for(const key of ['genitalCourseRepair','genitalProgress','genitalCurrentVascularFit','handPalmDelivery','aorticContinuity','ta98TermBindings','sourceCollectionIdentities','concepts','materials','defaultHidden'])assert.deepEqual(after[key],before[key]);
 const full=read(base+'blood-installed-v1334/packed-global-defects.json');assert.equal(full.passed,true);assert.equal(full.allPairsCompared,true);assert.equal(full.globalPairsBefore,6388);assert.equal(full.globalPairsAfter,6376);assert.equal(full.newGlobalPairs,0);assert.equal(full.remainingIntrinsicTargets,0);assert.equal(full.anatomicalAcceptance,false);
 const packed=read(base+'blood-installed-v1334/packed-neighbor-check.json');assert.equal(packed.passedParts,1);assert.ok(read(packed.report).checks.every(c=>c.passed&&Object.values(c.contacts).every(s=>s.newCrossingFacePairs===0)));assert.ok(read(base+'blood-installed-v1334/topology-retention.json').rows.every(c=>c.same));
});
test('Confluence repair uses all 26 source defects and retains the complete previous circulation',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/',before=read(base+'blood-installed-v1323/circulation-atlas.json'),after=read(base+'blood-installed-v1328/circulation-atlas.json'),r=after.bloodSurfaceRepairs;
 assert.equal(r.revision,'v1328');assert.equal(r.changedParts,1);assert.equal(r.removedIntrinsicPairs,26);assert.equal(r.remainingSelfCrossingParts,112);assert.equal(r.skeletonChanged,false);assert.equal(r.facesRemoved,0);
 const row=r.rows[0];assert.equal(row.partId,'BDM:male:human_20_male_cardiovascular_system-confluence_of_sinuses_ID');assert.equal(row.sourceDefectPairsComplete,true);assert.equal(row.sourceDefectPairs.length,26);assert.equal(new Set(row.sourceDefectPairs.map(p=>p.join(','))).size,26);
 const old=new Map(before.parts.map(p=>[p.id,p]));for(const p of after.parts){const prior=old.get(p.id);if(p.id!==row.partId)assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);assert.equal(p.ta98Term,prior.ta98Term);}}
 for(const key of ['genitalCourseRepair','genitalProgress','genitalCurrentVascularFit','handPalmDelivery','aorticContinuity','ta98TermBindings','sourceCollectionIdentities','concepts','materials','defaultHidden'])assert.deepEqual(after[key],before[key]);
 const packed=read(base+'blood-installed-v1328/packed-neighbor-check.json');assert.equal(packed.passedParts,1);assert.equal(packed.noNewCrossingFacePairs,true);assert.ok(read(packed.report).checks.every(c=>c.passed&&Object.values(c.contacts).every(s=>s.newCrossingFacePairs===0)));assert.ok(read(base+'blood-installed-v1328/topology-retention.json').rows.every(c=>c.same));
});
test('Larger vascular source repair batch preserves the complete installed circulation and source topology',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/',before=read(base+'genital-course-installed-v1317/circulation-atlas.json'),after=read(base+'blood-installed-v1323/circulation-atlas.json'),r=after.bloodSurfaceRepairs;
 assert.equal(r.revision,'v1323');assert.equal(r.changedParts,3);assert.equal(r.removedIntrinsicPairs,15);assert.equal(r.remainingSelfCrossingParts,113);assert.equal(r.skeletonChanged,false);assert.equal(r.facesRemoved,0);
 const changed=new Set(r.rows.map(p=>p.partId)),old=new Map(before.parts.map(p=>[p.id,p]));assert.equal(changed.size,3);assert.ok(changed.has('BDM:male:human_20_male_cardiovascular_system-splenic_artery_ID'));
 for(const p of after.parts){const prior=old.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);}}
 for(const key of ['genitalCourseRepair','genitalProgress','genitalCurrentVascularFit','handPalmDelivery','aorticContinuity','ta98TermBindings','sourceCollectionIdentities','concepts','materials','defaultHidden'])assert.deepEqual(after[key],before[key]);
 assert.equal(after.bloodSurfaceRepairHistory.at(-1).revision,'v1302');
 const proof=read(base+'blood-installed-v1323/packed-neighbor-check.json');assert.equal(proof.passedParts,3);assert.equal(proof.noNewCrossingFacePairs,true);assert.ok(read(proof.report).checks.every(c=>c.passed&&Object.values(c.contacts).every(x=>x.newCrossingFacePairs===0)));assert.ok(read(base+'blood-installed-v1323/topology-retention.json').rows.every(c=>c.same));
});
test('Bilateral cord course delivery retains source detail, fixed anatomy and verified branch interfaces',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const before=read(base+'blood-installed-v1302/circulation-atlas.json'),after=read(base+'genital-course-installed-v1317/circulation-atlas.json');
 const r=after.genitalCourseRepair,candidate=read(r.candidate),changed=new Set(candidate.rows.map(p=>p.partId)),old=new Map(before.parts.map(p=>[p.id,p]));
 assert.equal(r.revision,'v1317');assert.equal(r.changedParts,3);assert.equal(r.activeParts,2);assert.equal(changed.size,3);
 assert.equal(r.sourceFacesRetained,true);assert.equal(r.skeletonChanged,false);assert.equal(r.veinsChanged,false);assert.equal(r.organsChanged,false);assert.equal(r.classifierCapsDelivered,false);assert.equal(r.wholeAnatomicalAcceptance,false);
 for(const p of after.parts){const prior=old.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);assert.equal(p.suppressed,prior.suppressed);}}
 for(const key of ['concepts','materials','defaultHidden','ta98TermBindings','sourceCollectionIdentities','genitalProgress','genitalCurrentVascularFit','handPalmDelivery','aorticContinuity','bloodSurfaceRepairs','bloodSurfaceRepairHistory'])assert.deepEqual(after[key],before[key]);
 assert.equal(after.bloodSurfaceRepairs.remainingSelfCrossingParts,116);
 const prepare=read(base+'genital-course-installed-v1317/prepare-report.json'),approval=read(base+'genital-course-installed-v1317/packed-approval.json');
 assert.equal(prepare.complete,true);assert.equal(prepare.sourceClassifiersAbsentFromBank,true);assert.deepEqual(prepare.regionalArteryVeinPairs,{before:4753,after:4713});assert.ok(prepare.packedChecks.every(c=>c.selfCrossings===0&&c.packingError_m<1e-8));
 assert.equal(approval.complete,true);assert.equal(approval.passed,true);assert.equal(approval.analysisClassifierCapsDelivered,false);assert.equal(approval.wholeAnatomicalAcceptance,false);
 for(const result of approval.results){assert.equal(result.passed,true);assert.ok(result.maximumRelativeSectionContourAreaDifference<.05);for(const overlap of Object.values(result.wholeExtentOverlaps))assert.ok(overlap.passed&&overlap.after_mm3<=overlap.before_mm3+approval.overlapTolerance_mm3);}
 assert.equal(approval.interfaces.length,5);assert.ok(approval.interfaces.every(c=>c.positionsAndPackedNormalsLiteral));
});
test('Fourth vessel repair delivery retains the previous complete assembly and reciprocal aortic interface',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/',before=read(base+'blood-installed-v1292/circulation-atlas.json'),after=read(base+'blood-installed-v1302/circulation-atlas.json'),r=after.bloodSurfaceRepairs;
 assert.equal(r.revision,'v1302');assert.equal(r.changedParts,4);assert.equal(r.removedIntrinsicPairs,49);assert.equal(r.remainingSelfCrossingParts,116);assert.equal(r.skeletonChanged,false);assert.equal(r.facesRemoved,0);assert.equal(after.bloodSurfaceRepairHistory.at(-1).revision,'v1292');
 const changed=new Set(r.rows.map(p=>p.partId)),old=new Map(before.parts.map(p=>[p.id,p]));assert.equal(changed.size,4);assert.ok(changed.has('BDM:male:human_20_male_cardiovascular_system-ascending_branch_of_right_deep_circumflex_iliac_artery_ID'));
 for(const p of after.parts){const prior=old.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);}}
 for(const key of ['aorticContinuity','genitalProgress','handPalmDelivery','ta98TermBindings','sourceCollectionIdentities','concepts','materials','defaultHidden'])assert.deepEqual(after[key],before[key]);
 const packed=read(base+'blood-installed-v1302/packed-neighbor-check.json');assert.equal(packed.passedParts,4);assert.equal(packed.noNewCrossingFacePairs,true);assert.ok(read(packed.report).checks.every(c=>c.passed&&Object.values(c.contacts).every(s=>s.newCrossingFacePairs===0)));assert.ok(read(base+'blood-installed-v1302/topology-retention.json').rows.every(r=>r.same));
});
test('Current vessel surface delivery preserves source detail, fixed anatomy and all prior TA98 controls',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const before=read(base+'aortic-installed-v1283/circulation-atlas.json'),after=read(base+'blood-installed-v1292/circulation-atlas.json'),delivery=after.bloodSurfaceRepairs;
 assert.equal(delivery.revision,'v1292');assert.equal(delivery.changedParts,23);assert.equal(delivery.removedIntrinsicPairs,118);assert.equal(delivery.remainingSelfCrossingParts,120);assert.equal(delivery.facesRemoved,0);assert.equal(delivery.skeletonChanged,false);
 const changed=new Set(delivery.rows.map(p=>p.partId)),old=new Map(before.parts.map(p=>[p.id,p]));assert.equal(changed.size,23);assert.ok(changed.has('BDM:male:human_20_male_cardiovascular_system-right_proper_palmar_digital_arteries_ID'));
 for(const p of after.parts){const prior=old.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,prior);else{assert.equal(p.vertexCount,prior.vertexCount);assert.equal(p.indexCount,prior.indexCount);assert.equal(p.material,prior.material);}}
 for(const key of ['aorticContinuity','genitalProgress','handPalmDelivery','ta98TermBindings','sourceCollectionIdentities','concepts','materials','defaultHidden'])assert.deepEqual(after[key],before[key]);
 const packed=read(base+'blood-installed-v1292/packed-neighbor-check.json');assert.equal(packed.passedParts,23);assert.equal(packed.noNewCrossingFacePairs,true);const contacts=read(packed.report);assert.ok(contacts.checks.every(c=>c.passed&&Object.values(c.contacts).every(s=>s.newCrossingFacePairs===0)));
 const topo=read(base+'blood-installed-v1292/topology-retention.json');assert.equal(topo.parts,23);assert.ok(topo.rows.every(r=>r.same));
});
test('Original BioDigital heart review is local, explicit, shareable and takes precedence over blood-only review',()=>{
 const expected='/local-models/geometry-workbench/biodigital-male-circulation-v1/blood-installed-v1567/circulation-atlas.json';
 const query=new URLSearchParams('cardiac=biodigital-review&blood=biodigital-review');
 assert.equal(localTa98LayerCatalogue('local-ta98',query,true),expected);
 assert.equal(localTa98LayerCatalogue('local-ta98',query,false),undefined);
 assert.notEqual(localTa98LayerCatalogue('local-female-ta98',query,true),expected);
 assert.notEqual(localTa98LayerCatalogue('local-ta98',new URLSearchParams(),true),expected);
 for(const [model,development,preserved]of [['local-ta98',true,true],['local-ta98',false,false],['local-female-ta98',true,false]]){
  const url=new URL('https://example.invalid/?model='+model);
  preserveLocalTa98LayerCatalogueQuery(url,query,model,development);
  assert.equal(url.searchParams.get('cardiac'),preserved?'biodigital-review':null);
 }
});
test('BioDigital fit review is explicitly selected, male-only and development-only',()=>{
 const query=new URLSearchParams('blood=biodigital-review');
 const expected='/local-models/geometry-workbench/biodigital-male-circulation-v1/blood-installed-v1567/circulation-atlas.json';
 assert.equal(localTa98LayerCatalogue('local-ta98',query,true),expected);
 assert.equal(localTa98LayerCatalogue('local-ta98',query,false),undefined);
 assert.notEqual(localTa98LayerCatalogue('local-female-ta98',query,true),expected);
 assert.notEqual(localTa98LayerCatalogue('local-ta98',new URLSearchParams(),true),expected);
 assert.notEqual(localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=unknown'),true),expected);
 for(const [model,development,preserved]of [['local-ta98',true,true],['local-ta98',false,false],['local-female-ta98',true,false]]){
  const url=new URL('https://example.invalid/?model='+model);
  preserveLocalTa98LayerCatalogueQuery(url,query,model,development);
  assert.equal(url.searchParams.get('blood'),preserved?'biodigital-review':null);
 }
});
test('BioDigital circulation excludes retained legacy vessels from active anatomy and collision counts',()=>{
 const path=localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true);
 const file='.local-models/'+path.replace('/local-models/',''),atlas=read(file),hidden=new Set(atlas.defaultHidden);
 const blood=atlas.parts.filter(p=>p.system==='arterial'||p.system==='venous');
 const legacy=blood.filter(p=>!p.id.startsWith('BDM:'));
 assert.ok(legacy.length>0,'Legacy source geometry remains available for guidance');
 for(const p of legacy){assert.equal(p.suppressed,true,p.id);assert.equal(p.defaultHidden,true,p.id);assert.ok(hidden.has(p.id),p.id);}
 const active=blood.filter(p=>!p.suppressed&&!p.defaultHidden&&!hidden.has(p.id)),ids=new Set(active.map(p=>p.id));
 assert.ok(active.length>0);assert.ok(active.every(p=>p.id.startsWith('BDM:')));
 const genital=atlas.concepts.find(c=>c.id==='GENITAL:male:circulation');
 assert.ok(genital.elements.every(id=>ids.has(id)),'The genital circulation contains only active replacement sources');
 const context=file.slice(0,file.lastIndexOf('/'))+'/packed-context/report.json';
 if(fs.existsSync(context))for(const state of read(context).checks)for(const pair of state.arteryVein.pairs){
  assert.ok(ids.has(pair.artery)&&ids.has(pair.vein),'A collision references hidden guidance');
  assert.ok(pair.artery.startsWith('BDM:')&&pair.vein.startsWith('BDM:'));
 }
});
for(const sex of ['female','male'])test(`${sex} ordinary local atlas retains the installed intestinal coats and restores camera intent`,()=>{
 const model=sex==='male'?'local-ta98':'local-female-ta98',path=localTa98LayerCatalogue(model,new URLSearchParams(),true),atlas=read('.local-models/'+path.replace('/local-models/','')),prior=atlas,report=read(`reports/ta98-${sex}/small-intestine-local-atlas-pack.json`);
 assert.ok(atlas.bowelWallPrototype.ordinaryLocalAtlas);assert.equal(report.partIds.length,36);assert.ok(report.partIds.every(id=>atlas.bowelWallPrototype.partIds.includes(id)));
 const old=report.organs.find(o=>o.organ==='jejunum').oldSourceParts[0],camera=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5],q=new URLSearchParams({p:partFingerprint(prior),s:packPartIds([old],prior,false),frame:camera.join(','),cut:'coronal',slice:'.227'}),base={selected:[],visible:[],view:'three-quarter'};
 const restored=readViewUrl(q.toString(),atlas,base);assert.deepEqual(new Set(restored.selected),new Set(report.organs.find(o=>o.organ==='jejunum').newLayerParts));assert.deepEqual(restored.camera,camera);assert.equal(restored.section.position,.227);
 assert.equal(localTa98LayerCatalogue(model,new URLSearchParams('bowel=source'),true),undefined);assert.equal(localTa98LayerCatalogue(model,new URLSearchParams('delivery=0'),true),undefined);assert.equal(localTa98LayerCatalogue(model,new URLSearchParams(),false),undefined);
});
test('The user’s existing female URL fingerprint still resolves its exact hidden structures',()=>{
 const atlas=read('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'),prior=compatiblePartAtlas(atlas,'260.12sxffm.uizlzc');assert.ok(prior);assert.equal(partFingerprint(prior),'260.12sxffm.uizlzc');
 const q='p=260.12sxffm.uizlzc&h=jr.js.la-ms.mv-o3.o5.o7-oz.p1-p3.p5-px.pz-y8.153-15a.15c-15k.15m-161.163-169.1lq-1ls.1mb-1mw.1wk-1yh&l=136&skin=0&cut=coronal&slice=0.227&flip=1&frame=0.432952,0.603637,0.984524,0.565065,0.281339,0.505912,0.39865,0.480956,0.205037,0.5,1,0.5';
 const base={selected:[],visible:[],view:'three-quarter'},old=readViewUrl(q,prior,base),now=readViewUrl(q,atlas,base);assert.ok(old.hidden.length>0);assert.deepEqual(now.camera,old.camera);assert.deepEqual(now.visible,old.visible);for(const id of old.hidden){const p=atlas.parts.find(p=>p.id===id);for(const next of p.replacementParts??[id])assert.ok(now.hidden.includes(next));}
});


test('explicit packed pulmonary-root candidate is local-only and survives copying a review URL',()=>{
 const query=new URLSearchParams('blood=biodigital-root-review');
 const path=localTa98LayerCatalogue('local-ta98',query,true);
 assert.equal(path,'/local-models/geometry-workbench/biodigital-male-circulation-v1/blood-installed-v1567/circulation-atlas.json');
 assert.equal(localTa98LayerCatalogue('local-ta98',query,false),undefined);
 assert.notEqual(localTa98LayerCatalogue('local-female-ta98',query,true),path);
 const url=new URL('http://127.0.0.1:3016/?model=local-ta98');
 preserveLocalTa98LayerCatalogueQuery(url,query,'local-ta98',true);
 assert.equal(url.searchParams.get('blood'),'biodigital-root-review');
});

test('current triangulation delivery preserves every packed control and each changed patch boundary',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const path='.local-models/'+localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true).replace('/local-models/','');
 const current=read(path),proof=current.bloodTriangulationRepairs;
 const after=read(base+proof.revision+'/circulation-atlas.json'),before=read(proof.sourceAtlas??base+'blood-installed-v1451/circulation-atlas.json');
 const checks=read(base+proof.revision+'/packed-checks.json');
 assert.deepEqual(current.bloodTriangulationRepairs,after.bloodTriangulationRepairs);
 assert.ok(checks.complete&&checks.passed&&checks.packed);assert.equal(proof.changedParts,checks.rows.length);assert.ok(proof.changedParts>0);
 assert.equal(proof.removedIntrinsicPairs,checks.rows.reduce((sum,row)=>sum+row.intrinsicBefore-row.intrinsicAfter,0));assert.ok(proof.removedIntrinsicPairs>0);
 assert.ok(checks.rows.every(r=>r.newIntrinsicCrossingPairs===0&&r.newContextPairs===0&&r.intrinsicAfter<r.intrinsicBefore));
 const oldParts=new Map(before.parts.map(p=>[p.id,p])),selected=new Map(checks.rows.map(r=>[r.partId,r]));
 const banks=new Map(),bank=(atlas,part)=>{const url=atlas.chunks[part.chunk].url;assert.ok(url.startsWith('/local-models/'));if(!banks.has(url))banks.set(url,fs.readFileSync('.local-models/'+url.slice('/local-models/'.length)));return banks.get(url);};
 const restored=new Set((current.bloodHandSourceBaselineHistory??[]).flatMap(h=>h.partIds));
 for(const row of checks.rows){const now=current.parts.find(p=>p.id===row.partId),then=after.parts.find(p=>p.id===row.partId);
  if(!restored.has(row.partId)){assert.deepEqual(now,then,'retained triangulation '+row.partId);continue;}
  // A later reviewed source-shape restoration may move these controls, while
  // the original repaired diagonal and its raw control ancestry remain literal.
  const geometry=new Set(['positions','normals','normalSigns','uvs','colors','indices','chunk','bounds']);
  for(const key of Object.keys(then))if(!geometry.has(key))assert.deepEqual(now[key],then[key]);
  assert.deepEqual(bank(current,now).subarray(now.indices,now.indices+now.indexCount*4),bank(after,then).subarray(then.indices,then.indices+then.indexCount*4));
 }
 const boundary=(buffer,offset,faces)=>{
  const edges=new Map();
  for(const face of faces){const ids=[0,1,2].map(k=>buffer.readUInt32LE(offset+(face*3+k)*4));
   for(let i=0;i<3;i++){const a=ids[i],b=ids[(i+1)%3],key=a+','+b,reverse=b+','+a;
    if(edges.has(reverse)){const count=edges.get(reverse)-1;if(count)edges.set(reverse,count);else edges.delete(reverse);}
    else edges.set(key,(edges.get(key)||0)+1);
   }
  }return [...edges].sort(([a],[b])=>a.localeCompare(b));
 };
 const identityIds=new Set(after.bloodSourceIdentityAdditions?.at(-1)?.partIds??[]),withoutTags=p=>Object.fromEntries(Object.entries(p).filter(([k])=>!['ta98Term','ta98Kind','coverageLimitation'].includes(k)));
 for(const part of after.parts){const old=oldParts.get(part.id);assert.ok(old);
  if(!selected.has(part.id)){if(identityIds.has(part.id))assert.deepEqual(withoutTags(part),withoutTags(old));else assert.deepEqual(part,old);continue;}
  assert.equal(part.vertexCount,old.vertexCount);assert.equal(part.indexCount,old.indexCount);
  const original=bank(before,old),packed=bank(after,part);
  for(const [key,width]of [['positions',12],['normals',6],['normalSigns',1],['uvs',8],['colors',3]]){
   assert.equal(key in part,key in old);if(key in part)assert.deepEqual(packed.subarray(part[key],part[key]+part.vertexCount*width),original.subarray(old[key],old[key]+old.vertexCount*width),part.id+' '+key);
  }
  const faces=selected.get(part.id).changedFaces,changed=new Set(faces);
  for(let i=0;i<part.indexCount/3;i++)if(!changed.has(i))assert.deepEqual(packed.subarray(part.indices+i*12,part.indices+(i+1)*12),original.subarray(old.indices+i*12,old.indices+(i+1)*12));
  assert.deepEqual(boundary(packed,part.indices,faces),boundary(original,old.indices,faces),part.id+' patch boundary');
 }
 assert.equal(before.parts.length,after.parts.length);assert.deepEqual(after.materials,before.materials);assert.deepEqual(after.defaultHidden,before.defaultHidden);
 assert.deepEqual(after.bloodTriangulationRepairHistory.at(-1),before.bloodTriangulationRepairs);
});

test('partial hand course delivery preserves source branches, the complete inlet, fixed components and all other anatomy',()=>{
 const path='.local-models/'+localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true).replace('/local-models/','');
 const current=read(path),revision=current.bloodComponentCourseRepairHistory.at(-1).revision;
 const historical=path.replace(/blood-installed-v\d+\/circulation-atlas.json$/,revision+'/circulation-atlas.json');
 const after=read(historical),info=read(historical.replace('/circulation-atlas.json','/prepare-report.json'));
 const before=read(info.sourceCatalogue),contract=read(info.sourceReport.replace(/report\.json$/,'contract.json'));
 assert.ok(info.complete&&info.packed&&info.passedForLocalPreview);assert.equal(info.anatomicalAcceptance,false);
 assert.equal(info.newIntrinsicPairs,0);assert.equal(info.removedIntrinsicPairs,219);assert.ok(info.remainingLimits.length);
 const changed=new Set(info.changedPartIds),oldParts=new Map(before.parts.map(p=>[p.id,p]));assert.equal(changed.size,1);
 const banks=new Map(),bank=(atlas,p)=>{const url=atlas.chunks[p.chunk].url;if(!banks.has(url))banks.set(url,fs.readFileSync('.local-models/'+url.slice('/local-models/'.length)));return banks.get(url);};
 for(const p of after.parts){const old=oldParts.get(p.id);assert.ok(old);
  if(!changed.has(p.id)){assert.deepEqual(p,old);continue;}
  for(const key of ['vertexCount','indexCount','ta98Term','conceptId','material','sourceGroups','geometryOrigin','geometryScale','sourceOffset'])assert.deepEqual(p[key],old[key]);
  const a=bank(before,old),b=bank(after,p),moving=new Set(contract.movingRawIds);
  // This source uses every raw index. Bind compact geometry control IDs to
  // packed IDs before using the recorded source-component ownership.
  const used=new Set();for(let i=0;i<p.indexCount;i++)used.add(b.readUInt32LE(p.indices+i*4));assert.equal(used.size,p.vertexCount);
  assert.deepEqual(b.subarray(p.indices,p.indices+p.indexCount*4),a.subarray(old.indices,old.indices+old.indexCount*4));
  for(const [key,width]of [['positions',12],['normals',6],['normalSigns',1],['uvs',8],['colors',3]]){
   assert.equal(key in p,key in old);if(!(key in p))continue;
   const protectedIds=['positions','normals','normalSigns'].includes(key)?Array.from({length:p.vertexCount},(_,i)=>i).filter(i=>!moving.has(i)):Array.from({length:p.vertexCount},(_,i)=>i);
   for(const id of [...protectedIds,...contract.fixedMatchedInletRawIds])assert.deepEqual(b.subarray(p[key]+id*width,p[key]+(id+1)*width),a.subarray(old[key]+id*width,old[key]+(id+1)*width),p.id+' '+key+' protected '+id);
  }
 }
 for(const key of ['concepts','materials','defaultHidden','ta98TermBindings','sourceCollectionIdentities'])assert.deepEqual(after[key],before[key]);
 const last=after.bloodComponentCourseRepairHistory.at(-1);assert.equal(last.removedIntrinsicPairs,219);assert.equal(last.repairedComponentSelfPairs,0);assert.equal(last.skeletonChanged,false);assert.equal(last.anatomicalAcceptance,false);
});

test('coupled arterial inlet delivery preserves both owners, source detail and the fixed skeleton',()=>{
 const path='.local-models/'+localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true).replace('/local-models/','');
 const current=read(path),last=current.bloodCoupledCourseRepairHistory.at(-1);
 const snapshot=path.replace(/blood-installed-v\d+\/circulation-atlas.json$/,last.revision+'/circulation-atlas.json');
 const after=read(snapshot);
 const info=read(snapshot.replace('/circulation-atlas.json','/install-report.json')),before=read(info.sourceCatalogue);
 assert.ok(info.complete&&info.installed&&info.fixedSkeletonPartRecordsLiteral&&info.unchangedPartsLiteral);
 assert.equal(info.removedIntrinsicPairs,7);assert.equal(info.newIntrinsicPairs,0);assert.equal(info.anatomicalAcceptance,false);
 assert.equal(info.sourceMouth.maximumGap_m,0);assert.ok(info.sourceMouth.completeSourceBoundaries);
 assert.equal(info.sourceMouth.wholeLumenPatencyEstablished,false);assert.equal(last.partIds.length,2);
 const changed=new Set(last.partIds),original=new Map(before.parts.map(p=>[p.id,p]));
 for(const p of after.parts){const old=original.get(p.id);assert.ok(old);if(!changed.has(p.id)){assert.deepEqual(p,old);continue;}
  for(const k of ['vertexCount','indexCount','material','conceptId','ta98Term','sourceGroups','geometryOrigin','geometryScale','sourceOffset'])assert.deepEqual(p[k],old[k]);
  const bank=(atlas,part)=>fs.readFileSync('.local-models/'+atlas.chunks[part.chunk].url.slice('/local-models/'.length));
  const a=bank(before,old),b=bank(after,p);
  assert.deepEqual(b.subarray(p.indices,p.indices+p.indexCount*4),a.subarray(old.indices,old.indices+old.indexCount*4));
  for(const [k,width]of [['uvs',8],['colors',3]]){assert.equal(k in p,k in old);if(k in p)assert.deepEqual(b.subarray(p[k],p[k]+p.vertexCount*width),a.subarray(old[k],old[k]+old.vertexCount*width));}
 }
 for(const k of ['concepts','materials','defaultHidden','ta98TermBindings','sourceCollectionIdentities','bloodComponentCourseRepairHistory'])assert.deepEqual(after[k],before[k]);
 const packed=read(snapshot.replace('/circulation-atlas.json','/packed-review.json'));assert.ok(packed.complete&&packed.actualPackedContactsMatchReviewedProposal);
 for(const context of packed.contextReports)for(const system of ['skeletal','arterial','venous','integumentary'])assert.ok(context.contactsBySystem[system][1]<=context.contactsBySystem[system][0]);
});

test('bilateral source restoration retains all detailed hand networks, proximal words and fixed anatomy',()=>{
 const path='.local-models/'+localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true).replace('/local-models/','');
 const current=read(path),last=current.bloodHandSourceBaselineHistory.at(-1);
 const snapshot=path.replace(/blood-installed-v\d+\/circulation-atlas.json$/,last.revision+'/circulation-atlas.json');
 const after=read(snapshot),info=read(snapshot.replace('/circulation-atlas.json','/install-report.json')),before=read(info.sourceCatalogue),proposal=read(info.sourceReport);
 assert.ok(info.complete&&info.installed&&info.restorationOnly&&info.fixedSkeletonLiteral&&info.allOtherAnatomyLiteral);
 assert.equal(info.beforeSelfPairs,10518);assert.equal(info.afterSelfPairs,0);assert.equal(info.bodyPackingAccepted,false);assert.equal(info.anatomicalAcceptance,false);
 const rows=new Map(proposal.rows.map(r=>[r.partId,r])),old=new Map(before.parts.map(p=>[p.id,p]));assert.equal(rows.size,52);
 for(const side of ['left','right'])assert.equal(proposal.rows.filter(r=>r.side===side).length,26);
 const banks=new Map(),bank=(atlas,p)=>{const url=atlas.chunks[p.chunk].url;if(!banks.has(url))banks.set(url,fs.readFileSync('.local-models/'+url.slice('/local-models/'.length)));return banks.get(url);};
 for(const p of after.parts){const original=old.get(p.id);assert.ok(original);const r=rows.get(p.id);if(!r){assert.deepEqual(p,original);continue;}
  for(const key of ['vertexCount','indexCount','conceptId','material','ta98Term','ta98Kind','sourceGroups','geometryOrigin','geometryScale','sourceOffset'])assert.deepEqual(p[key],original[key]);
  const a=bank(before,original),b=bank(after,p);assert.deepEqual(b.subarray(p.indices,p.indices+p.indexCount*4),a.subarray(original.indices,original.indices+original.indexCount*4));
  const used=[...new Set(Array.from({length:p.indexCount},(_,i)=>b.readUInt32LE(p.indices+i*4)))].sort((a,b)=>a-b);
  for(const id of r.fixedRawIds){const raw=used[id];assert.notEqual(raw,undefined);for(const [key,width]of [['positions',12],['normals',6],['normalSigns',1]])assert.deepEqual(b.subarray(p[key]+raw*width,p[key]+(raw+1)*width),a.subarray(original[key]+raw*width,original[key]+(raw+1)*width));}
  for(const [key,width]of [['uvs',8],['colors',3]]){assert.equal(key in p,key in original);if(key in p)assert.deepEqual(b.subarray(p[key],p[key]+p.vertexCount*width),a.subarray(original[key],original[key]+original.vertexCount*width));}
 }
 for(const key of ['concepts','materials','defaultHidden','ta98TermBindings','sourceCollectionIdentities','bloodCoupledCourseRepairHistory','bloodComponentCourseRepairHistory'])assert.deepEqual(after[key],before[key]);
 const packed=read(snapshot.replace('/circulation-atlas.json','/packed-candidate-report.json'));assert.equal(packed.rows.length,52);assert.ok(packed.rows.every(r=>r.afterSelfPairs===0));
 const ports=read(snapshot.replace('/circulation-atlas.json','/source-mouth-replay.json'));assert.equal(ports.summary.verified,23);assert.equal(ports.summary.pending,0);
});

test('compact hand clearance retains imported detail, hierarchy and fixed anatomy with packed contact evidence',()=>{
 const path='.local-models/'+localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true).replace('/local-models/','');
 const current=read(path),revision=current.bloodCompactClearanceHistory.at(-1).revision;
 const snapshot=path.replace(/blood-installed-v\d+\/circulation-atlas.json$/,revision+'/circulation-atlas.json');
 const after=read(snapshot),info=read(snapshot.replace('/circulation-atlas.json','/install-report.json')),before=read(info.sourceCatalogue);
 assert.ok(info.complete&&info.packed&&info.installed&&info.allOtherAnatomyLiteral&&info.sourceFacesUvsColorsMaterialsLiteral);
 assert.equal(info.anatomicalAcceptance,false);assert.equal(info.sourceMouthSummary.verified,23);assert.equal(info.sourceMouthSummary.pending,0);
 const proposal=read(info.sourceReport);assert.ok(info.packedBody.skinOutside<proposal.initial.skinOutside);assert.ok(info.packedBody.boneInside<=proposal.initial.boneInside);assert.ok(info.packedBody.energy_m2<proposal.initial.energy_m2);
 const changed=new Set(info.changedPartIds),old=new Map(before.parts.map(p=>[p.id,p]));assert.equal(changed.size,4);
 const banks=new Map(),bank=(atlas,p)=>{const url=atlas.chunks[p.chunk].url;if(!banks.has(url))banks.set(url,fs.readFileSync('.local-models/'+url.slice('/local-models/'.length)));return banks.get(url);};
 for(const p of after.parts){const q=old.get(p.id);assert.ok(q);if(!changed.has(p.id)){assert.deepEqual(p,q);continue;}
  const geometry=new Set(['positions','normals','normalSigns','uvs','colors','indices','chunk','bounds']);
  for(const key of Object.keys(q))if(!geometry.has(key))assert.deepEqual(p[key],q[key]);
  const a=bank(before,q),b=bank(after,p);assert.deepEqual(b.subarray(p.indices,p.indices+p.indexCount*4),a.subarray(q.indices,q.indices+q.indexCount*4));
  for(const [key,width]of [['uvs',8],['colors',3]]){assert.equal(key in p,key in q);if(key in p)assert.deepEqual(b.subarray(p[key],p[key]+p.vertexCount*width),a.subarray(q[key],q[key]+q.vertexCount*width));}
 }
 for(const key of Object.keys(before))if(!['parts','chunks','bloodCompactClearanceHistory'].includes(key))assert.deepEqual(after[key],before[key]);
 const review=read(snapshot.replace('/circulation-atlas.json','/packed-review.json'));assert.equal(review.contexts.length,4);
 for(const ref of review.contexts){const evidence=read(ref.path);assert.ok(evidence.complete&&evidence.candidateContextConsistent);assert.ok(evidence.contextLimit>=evidence.activeContextCount);assert.ok(evidence.rows.every(x=>x.after<=x.before));}
});

test('reviewed pancreatic source resolves the splenic occurrence and remains selectable under its recipient',()=>{
 const path='.local-models/'+localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true).replace('/local-models/','');
 const atlas=withReviewedCatalogueIdentities(read(path)),id='BDM:male:human_20_male_cardiovascular_system-pancreatic_veins_ID';
 const part=atlas.parts.find(p=>p.id===id),concept=atlas.concepts.find(c=>c.id===part.conceptId);
 assert.equal(part.ta98Term,'A12.3.12.029');assert.equal(part.ta98Kind,'exact');assert.equal(concept.ta98Term,part.ta98Term);
 const recipient=atlas.concepts.find(c=>c.id===concept.partOf);assert.equal(recipient.ta98Term,'A12.3.12.028');
 const bindings=catalogueTa98TermParts(atlas);assert.ok(bindings.get('A12.3.12.029')?.some(p=>p.id===id));
 assert.ok(!bindings.get('A12.3.12.022')?.some(p=>p.id===id));
 const parts=new Map(atlas.parts.filter(p=>!p.suppressed).map(p=>[p.id,p]));
 assert.ok(!anatomicalComponentIds(recipient,[],parts,atlas.concepts).includes(id));
 const entries=hierarchyEntries(atlas),entry=entries.find(e=>e.id===concept.id),parentEntry=entries.find(e=>e.id===recipient.id);
 assert.equal(entry.taParentConcept,parentEntry.id);assert.deepEqual(systemPathFor(entry),systemPathFor(parentEntry));assert.equal(systemPathFor(entry).at(-1),parentEntry.name);
 const routes=pancreaticVenousRelations(atlas.parts);
 assert.ok(routes.get(id).some(r=>r.kind==='after'&&r.target.id===recipient.elements[0]));
 assert.ok(routes.get(recipient.elements[0]).some(r=>r.kind==='before'&&r.target.id===id));
 assert.equal(reviewedPancreaticVenousPair(atlas.parts.map(p=>p.id===id?{...p,ta98Term:undefined}:p)),undefined);
 assert.equal(reviewedPancreaticVenousPair(atlas.parts.filter(p=>p.id!==recipient.elements[0])),undefined);
 const camera=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5];
 const restored=readViewUrl('select='+encodeURIComponent(id)+'&frame='+camera.join(','),atlas,{selected:[],visible:[],view:'front'});
 assert.deepEqual(restored.selected,[id]);assert.deepEqual(restored.camera,camera);
});

test('genital branch partitions expose six verified controls, preserve old selections and hide duplicate sources',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const atlas=withReviewedCatalogueIdentities(read(base+'genital-installed-v1219/circulation-atlas.json'));
 const previous=read(base+'genital-installed-v1196/circulation-atlas.json');
 const proof=atlas.genitalCurrentVascularFit.ta98Controls,bindings=catalogueTa98TermParts(atlas);
 assert.equal(proof.resolvedTerms.length,6);assert.equal(proof.newPartIds.length,17);
 for(const term of proof.resolvedTerms)assert.ok(bindings.get(term)?.length,term+' must resolve to actual active meshes');
 assert.deepEqual(proof.unresolvedTerms,['A12.2.12.087M','A12.2.16.016']);
 const parts=new Map(atlas.parts.filter(p=>!p.suppressed).map(p=>[p.id,p]));
 const camera=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5];
 for(const [id,replacements]of Object.entries(proof.replacements)){
  const source=atlas.parts.find(p=>p.id===id);assert.ok(source.suppressed&&atlas.defaultHidden.includes(id));
  const concept=atlas.concepts.find(c=>c.id===id);assert.deepEqual(concept.elements,replacements);
  assert.deepEqual(new Set(anatomicalComponentIds(concept,[],parts,atlas.concepts)),new Set(replacements));
  const restored=readViewUrl('select='+encodeURIComponent(id)+'&frame='+camera.join(','),atlas,{selected:[],visible:[],view:'front'});
  assert.deepEqual(new Set(restored.selected),new Set(replacements));assert.deepEqual(restored.camera,camera);
  assert.equal(replacements.reduce((n,r)=>n+parts.get(r).indexCount,0),source.indexCount);
 }
 const region=atlas.concepts.find(c=>c.id==='GENITAL:male:circulation');assert.equal(region.elements.length,94);assert.equal(new Set(region.elements).size,94);assert.ok(region.elements.every(id=>parts.has(id)));
 assert.ok(compatiblePartAtlas(atlas,partFingerprint(previous)));
 const originals=new Map(previous.parts.map(p=>[p.id,p]));for(const p of atlas.parts.filter(p=>p.system==='skeletal'))assert.deepEqual(p,originals.get(p.id));
 assert.deepEqual(atlas.materials,previous.materials);
});

test('constructed bilateral genital rami have actual shared parent portals, four selectable meshes and both new TA98 controls',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const raw=read(base+'genital-installed-v1226/circulation-atlas.json'),atlas=withReviewedCatalogueIdentities(raw),prior=read(base+'genital-installed-v1219/circulation-atlas.json');
 const proof=atlas.genitalBranchDelivery,bindings=catalogueTa98TermParts(atlas),parts=new Map(atlas.parts.filter(p=>!p.suppressed).map(p=>[p.id,p]));
 assert.equal(proof.newPartIds.length,4);assert.equal(proof.branches.length,4);
 for(const term of ['A12.2.12.087M','A12.2.16.016'])assert.equal(bindings.get(term)?.filter(p=>proof.newPartIds.includes(p.id)).length,2,term);
 const bank=fs.readFileSync(base+'genital-installed-v1226/geometry.bin');
 for(const branch of proof.branches){
  const p=parts.get(branch.partId),d=parts.get(branch.branchPartId);assert.ok(d.vertexCount>12&&d.indexCount>60);assert.equal(d.ta98Term,branch.term);
  assert.match(d.coverageLimitation,/Constructed/);assert.equal(branch.lumenDisksCreated,0);
  for(let j=0;j<6;j++){
   const pi=branch.parentPortalRawIds[j],di=branch.daughterPortalRawIds[j];
   assert.deepEqual(bank.subarray(p.positions+pi*12,p.positions+(pi+1)*12),bank.subarray(d.positions+di*12,d.positions+(di+1)*12));
   assert.deepEqual(bank.subarray(p.normals+pi*6,p.normals+(pi+1)*6),bank.subarray(d.normals+di*6,d.normals+(di+1)*6));
  }
  const parent=atlas.concepts.find(c=>c.id===p.conceptId);assert.ok(anatomicalComponentIds(parent,[],parts,atlas.concepts).includes(d.id));
 }
 const group=atlas.concepts.find(c=>c.id==='GENITAL:male:circulation');assert.equal(group.elements.length,98);assert.equal(new Set(group.elements).size,98);assert.ok(proof.newPartIds.every(id=>group.elements.includes(id)));
 const old=new Map(prior.parts.map(p=>[p.id,p]));for(const p of atlas.parts.filter(p=>p.system==='skeletal'))assert.deepEqual(p,old.get(p.id));
 assert.deepEqual(atlas.handPalmDelivery,prior.handPalmDelivery);assert.ok(compatiblePartAtlas(atlas,partFingerprint(prior)));
});

test('current genital parent controls keep parent arteries and all eight added terms selectable without duplicate geometry',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const a=withReviewedCatalogueIdentities(read(base+'genital-installed-v1229/circulation-atlas.json')),prior=withReviewedCatalogueIdentities(read(base+'genital-installed-v1226/circulation-atlas.json'));
 const bindings=catalogueTa98TermParts(a),old=catalogueTa98TermParts(prior),parts=new Map(a.parts.filter(p=>!p.suppressed).map(p=>[p.id,p]));
 for(const term of ['A12.2.12.086M','A12.2.16.014'])assert.equal(bindings.get(term)?.length,2,term+' retains both parent arterial courses');
 for(const [term,meshes]of old)assert.ok(meshes.every(p=>bindings.get(term)?.some(q=>q.id===p.id)),term+' existing binding retained');
 for(const term of [...a.genitalCurrentVascularFit.ta98Controls.resolvedTerms,'A12.2.12.087M','A12.2.16.016'])assert.ok(bindings.get(term)?.length,term);
 for(const branch of a.genitalBranchDelivery.branches){const child=a.concepts.find(c=>c.id===branch.branchPartId),parent=a.concepts.find(c=>c.id===child.partOf);assert.ok(anatomicalComponentIds(parent,[],parts,a.concepts).includes(branch.branchPartId));}
 assert.deepEqual(a.chunks,prior.chunks);assert.equal(a.parts.length,prior.parts.length);assert.equal(new Set(a.parts.map(p=>p.id)).size,a.parts.length);
 for(let i=0;i<a.parts.length;i++){const {conceptId:ignored,...p}=a.parts[i],{conceptId:oldIgnored,...q}=prior.parts[i];assert.deepEqual(p,q);}
});

test('blood term batch resolves its 42 reviewed meshes through ten TA98 controls while retaining previous controls',t=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const a=withReviewedCatalogueIdentities(read(base+'blood-installed-v1237/circulation-atlas.json')),prior=withReviewedCatalogueIdentities(read(base+'genital-installed-v1229/circulation-atlas.json'));
 const current=catalogueTa98TermParts(a),previous=catalogueTa98TermParts(prior),batch=a.bloodTermBatch;
 assert.deepEqual(batch.terms.filter(term=>!previous.has(term)),['A12.2.08.058','A12.3.10.010','A12.3.11.014']);
 t.diagnostic(JSON.stringify({newlyResolvedTerms:batch.terms.filter(term=>!previous.has(term)),membersByTerm:batch.terms.map(term=>({term,before:previous.get(term)?.length??0,after:current.get(term)?.length??0}))}));
 assert.equal(batch.terms.length,10);assert.equal(batch.mappedPartIds.length,42);assert.equal(batch.selfCrossingsBefore,23);assert.equal(batch.selfCrossingsAfter,0);
 for(const binding of batch.bindings)for(const id of binding.partIds)assert.ok(current.get(binding.term)?.some(p=>p.id===id),binding.term+' '+id);
 for(const [term,meshes]of previous)for(const p of meshes)assert.ok(current.get(term)?.some(q=>q.id===p.id),term+' previous member retained: '+p.id);
 const foot=current.get('A12.3.11.014').filter(p=>batch.mappedPartIds.includes(p.id));assert.equal(foot.length,2);assert.ok(foot.every(p=>p.sourceGroups.some(g=>g.includes('dorsal_veins_of_')&&g.includes('_foot_ID'))));
 assert.ok(batch.heldAnteriorSaphenousSourceIds.length===2);assert.ok(batch.heldAnteriorSaphenousSourceIds.every(id=>!batch.mappedPartIds.includes(id)));
 const parts=new Map(a.parts.map(p=>[p.id,p]));for(const binding of batch.bindings){const c=a.concepts.find(c=>c.id===binding.conceptId);if(binding.partIds.length>1)assert.deepEqual(new Set(anatomicalComponentIds(c,[],parts,a.concepts)),new Set(binding.partIds));}
 const original=new Map(prior.parts.map(p=>[p.id,p]));assert.equal(a.parts.length,prior.parts.length);assert.equal(a.materials.length,prior.materials.length);
 for(const p of a.parts.filter(p=>p.system==='skeletal'))assert.deepEqual(p,original.get(p.id));
 const genital=prior.concepts.find(c=>c.id==='GENITAL:male:circulation').elements;for(const id of genital)assert.deepEqual(parts.get(id),original.get(id));
 assert.deepEqual(a.handPalmDelivery,prior.handPalmDelivery);assert.deepEqual(a.defaultHidden,prior.defaultHidden);assert.ok(compatiblePartAtlas(a,partFingerprint(prior)));
});

test('genital progression retains all previous controls, adds right epididymal mesh and correctly distinguishes cerebral branch parents',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const a=withReviewedCatalogueIdentities(read(base+'genital-installed-v1257/circulation-atlas.json')),prior=withReviewedCatalogueIdentities(read(base+'blood-installed-v1237/circulation-atlas.json')),now=catalogueTa98TermParts(a),old=catalogueTa98TermParts(prior),delivery=a.genitalProgress;
 for(const [term,meshes]of old)for(const p of meshes)assert.ok(now.get(term)?.some(q=>q.id===p.id),term+' retained '+p.id);
 for(const term of ['A12.2.07.042','A12.2.07.099'])assert.equal(now.get(term)?.filter(p=>p.id.startsWith('BDM:')).length,2,term);
 const branch=a.parts.find(p=>p.id===delivery.branch.branchPartId);assert.ok(now.get('A12.2.12.088M').some(p=>p.id===branch.id));assert.equal(branch.name,'Epididymal branch (right)');assert.match(branch.coverageLimitation,/Constructed/);
 assert.ok(a.parts.some(p=>p.ta98Term==='A12.2.12.088M'&&p.id.includes('left_')));
 const parts=new Map(a.parts.map(p=>[p.id,p]));for(const side of ['left','right']){const child=a.concepts.find(c=>c.id.includes('set_of_'+side+'_paracentral_branches_ID'));assert.equal(child.ta98Term,'A12.2.07.042');assert.equal(child.partOf,'BDM:male:human_20_male_cardiovascular_system-'+side+'_callosomarginal_artery_ID');const parent=a.concepts.find(c=>c.id===child.partOf);assert.equal(parent.ta98Term,'A12.2.07.037');assert.ok(parts.has(child.elements[0]));assert.equal(now.get('A12.2.07.042').filter(p=>p.id===child.elements[0]).length,1);}
 const group=a.concepts.find(c=>c.id==='GENITAL:male:circulation');assert.equal(group.elements.length,99);assert.equal(new Set(group.elements).size,99);assert.ok(group.elements.includes(branch.id));
 const previous=new Map(prior.parts.map(p=>[p.id,p]));for(const p of a.parts.filter(p=>p.system==='skeletal'))assert.deepEqual(p,previous.get(p.id));assert.deepEqual(a.materials,prior.materials);assert.deepEqual(a.defaultHidden,prior.defaultHidden);assert.deepEqual(a.handPalmDelivery,prior.handPalmDelivery);assert.ok(compatiblePartAtlas(a,partFingerprint(prior)));
 const bank=fs.readFileSync(base+'genital-installed-v1257/geometry.bin'),parent=parts.get(delivery.branch.partId);for(let j=0;j<6;j++){const pi=delivery.branch.parentPortalRawIds[j],di=delivery.branch.daughterPortalRawIds[j];assert.deepEqual(bank.subarray(parent.positions+12*pi,parent.positions+12*(pi+1)),bank.subarray(branch.positions+12*di,branch.positions+12*(di+1)));}
 assert.ok(delivery.metrics.skinAfter<delivery.metrics.skinBefore);assert.ok(delivery.metrics.arteryVeinAfter<delivery.metrics.arteryVeinBefore);
});


test('reviewed blood surface batch preserves all TA98 controls, source detail, fixed anatomy and old URLs',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const before=withReviewedCatalogueIdentities(read(base+'genital-installed-v1257/circulation-atlas.json')),after=withReviewedCatalogueIdentities(read(base+'blood-installed-v1269/circulation-atlas.json')),delivery=after.bloodSurfaceRepairs;
 assert.equal(delivery.changedParts,16);assert.equal(delivery.removedIntrinsicPairs,43);assert.equal(delivery.facesRemoved,0);assert.equal(delivery.skeletonChanged,false);assert.equal(delivery.wholeCirculationAccepted,false);
 assert.deepEqual(after.parts.map(p=>p.id),before.parts.map(p=>p.id));assert.deepEqual(after.concepts,before.concepts);assert.deepEqual(after.defaultHidden,before.defaultHidden);assert.deepEqual(after.materials,before.materials);assert.deepEqual(after.genitalProgress,before.genitalProgress);assert.deepEqual(after.handPalmDelivery,before.handPalmDelivery);assert.ok(compatiblePartAtlas(after,partFingerprint(before)));
 const old=catalogueTa98TermParts(before),now=catalogueTa98TermParts(after);assert.deepEqual([...now.keys()].sort(),[...old.keys()].sort());for(const [term,members]of old)assert.deepEqual(now.get(term).map(p=>p.id).sort(),members.map(p=>p.id).sort());
 const changed=new Set(delivery.rows.map(r=>r.partId)),previous=new Map(before.parts.map(p=>[p.id,p]));for(const p of after.parts){const old=previous.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,old);else{assert.equal(p.vertexCount,old.vertexCount);assert.equal(p.indexCount,old.indexCount);assert.equal(p.conceptId,old.conceptId);assert.equal(p.ta98Term,old.ta98Term);assert.equal(p.system,old.system);}}
 const proof=read(base+'blood-installed-v1269/prepare-report.json');assert.equal(proof.packedChecks.length,16);assert.ok(proof.packedChecks.every(p=>p.selfCrossings===0&&p.packingError_m<1e-8));
});
test('second reviewed blood batch retains prior repairs, TA98 controls and fixed anatomy',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const before=withReviewedCatalogueIdentities(read(base+'blood-installed-v1269/circulation-atlas.json')),after=withReviewedCatalogueIdentities(read(base+'blood-installed-v1276/circulation-atlas.json')),delivery=after.bloodSurfaceRepairs;
 assert.equal(delivery.changedParts,22);assert.equal(delivery.removedIntrinsicPairs,133);assert.equal(delivery.facesRemoved,0);assert.equal(delivery.skeletonChanged,false);assert.equal(delivery.wholeCirculationAccepted,false);
 assert.deepEqual(after.parts.map(p=>p.id),before.parts.map(p=>p.id));assert.deepEqual(after.concepts,before.concepts);assert.deepEqual(after.defaultHidden,before.defaultHidden);assert.deepEqual(after.materials,before.materials);assert.deepEqual(after.genitalProgress,before.genitalProgress);assert.deepEqual(after.handPalmDelivery,before.handPalmDelivery);assert.ok(compatiblePartAtlas(after,partFingerprint(before)));
 const old=catalogueTa98TermParts(before),now=catalogueTa98TermParts(after);assert.deepEqual([...now.keys()].sort(),[...old.keys()].sort());for(const [term,members]of old)assert.deepEqual(now.get(term).map(p=>p.id).sort(),members.map(p=>p.id).sort());
 const changed=new Set(delivery.rows.map(r=>r.partId)),previous=new Map(before.parts.map(p=>[p.id,p]));for(const p of after.parts){const old=previous.get(p.id);if(!changed.has(p.id))assert.deepEqual(p,old);else{assert.equal(p.vertexCount,old.vertexCount);assert.equal(p.indexCount,old.indexCount);assert.equal(p.conceptId,old.conceptId);assert.equal(p.ta98Term,old.ta98Term);assert.equal(p.system,old.system);}}
 assert.equal(after.bloodSurfaceRepairHistory.at(-1).revision,'v1269');assert.equal(delivery.remainingSelfCrossingParts,143);const proof=read(base+'blood-installed-v1276/prepare-report.json');assert.equal(proof.packedChecks.length,22);assert.ok(proof.packedChecks.every(p=>p.selfCrossings===0&&p.packingError_m<1e-8));
});


test('renal cortical-radiate controls resolve the four imported networks without changing their geometry or prior repairs',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
 const before=withReviewedCatalogueIdentities(read(base+'blood-installed-v1339/circulation-atlas.json')),a=withReviewedCatalogueIdentities(read(base+'blood-installed-v1346/circulation-atlas.json'));
 const original=new Map(before.parts.map(p=>[p.id,p])),concepts=new Map(a.concepts.map(c=>[c.id,c])),terms=catalogueTa98TermParts(a),ids=new Set(a.renalVascularIdentityDelivery.mappedPartIds);
 assert.equal(ids.size,4);assert.deepEqual(a.chunks,before.chunks);assert.deepEqual(a.defaultHidden,before.defaultHidden);assert.deepEqual(a.materials,before.materials);assert.deepEqual(a.bloodSurfaceRepairs,before.bloodSurfaceRepairs);assert.deepEqual(a.genitalProgress,before.genitalProgress);
 for(const p of a.parts){const old=original.get(p.id);if(!ids.has(p.id)){assert.deepEqual(p,old);continue;}
  const term=p.system==='arterial'?'A08.1.03.004':'A08.1.04.004';assert.equal(p.ta98Term,term);assert.equal(concepts.get(p.conceptId).ta98Term,term);assert.ok(terms.get(term).some(q=>q.id===p.id));
  for(const key of ['chunk','positions','normals','normalSigns','indices','vertexCount','indexCount','uvs','colors','geometryOrigin','sourceGroups'])assert.deepEqual(p[key],old[key]);
  assert.equal(concepts.get(p.conceptId).partOf,before.concepts.find(c=>c.id===p.conceptId).partOf);
 }
});

test('aortic source cells keep their lumen open, hide only reciprocal annular interfaces and retain every prior selector',()=>{
 const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';const before=withReviewedCatalogueIdentities(read(base+'blood-installed-v1276/circulation-atlas.json')),a=withReviewedCatalogueIdentities(read(base+'aortic-installed-v1283/circulation-atlas.json')),delivery=a.aorticContinuity;
 assert.equal(delivery.junctionCrossingsBefore,141);assert.equal(delivery.junctionCrossingsAfter,0);assert.ok(delivery.minimumLocalGuideClearance_m>.004);assert.equal(delivery.remainingBodyPlacementDefects,true);assert.equal(delivery.skeletonChanged,false);assert.equal(delivery.otherVesselsChanged,false);
 const old=catalogueTa98TermParts(before),now=catalogueTa98TermParts(a);assert.deepEqual([...now.keys()].sort(),[...old.keys()].sort());for(const [term,members]of old)assert.deepEqual(now.get(term).map(p=>p.id).sort(),members.map(p=>p.id).sort());assert.deepEqual(a.concepts,before.concepts);assert.deepEqual(a.parts.map(p=>p.id),before.parts.map(p=>p.id));assert.ok(compatiblePartAtlas(a,partFingerprint(before)));
 const changed=new Set(delivery.parts.map(p=>p.partId)),previous=new Map(before.parts.map(p=>[p.id,p]));for(const p of a.parts)if(!changed.has(p.id))assert.deepEqual(p,previous.get(p.id));assert.deepEqual(a.bloodSurfaceRepairs,before.bloodSurfaceRepairs);assert.deepEqual(a.bloodSurfaceRepairHistory,before.bloodSurfaceRepairHistory);assert.deepEqual(a.genitalProgress,before.genitalProgress);
 const bank=fs.readFileSync(base+'aortic-installed-v1283/geometry.bin'),parts=delivery.parts.map(row=>a.parts.find(p=>p.id===row.partId));const surfaceKeys=[];
 for(const p of parts){assert.equal(p.sectionSolid,true);assert.equal(p.pairedTissue,true);assert.equal(p.tissueDrawGroups.length,2);const cap=p.tissueDrawGroups[1],other=parts.find(q=>q!==p);assert.equal(cap.neighbor,other.id);assert.equal(cap.count,delivery.sharedInterfaceFaces*3);assert.equal(p.tissueDrawGroups[0].count+cap.count,p.indexCount);
  const keys=[];for(let j=cap.start;j<cap.start+cap.count;j+=3){const tri=[];for(let k=0;k<3;k++){const id=bank.readUInt32LE(p.indices+4*(j+k));tri.push([0,1,2].map(axis=>bank.readFloatLE(p.positions+12*id+4*axis)+p.geometryOrigin[axis]));}keys.push(tri.map(v=>v.join(',')).sort().join(';'));
   const [x,z]=[.010,-.007],u=tri.map(v=>[v[0]-x,v[2]-z]);const cross=(a,b)=>a[0]*b[1]-a[1]*b[0],s=[cross(u[0],u[1]),cross(u[1],u[2]),cross(u[2],u[0])];assert.ok(!(s.every(v=>v>=-1e-18)||s.every(v=>v<=1e-18)),'No interface triangle fills the measured lumen center');}
  surfaceKeys.push(keys.sort());
 }
 assert.deepEqual(surfaceKeys[0],surfaceKeys[1]);const matrix=p=>({elements:[1,0,0,0,0,1,0,0,0,0,1,0,...p.geometryOrigin,1]});assert.ok(tissuePosesCoincide(matrix(parts[0]),parts[0].geometryOrigin,matrix(parts[1]),parts[1].geometryOrigin));assert.ok(tissueNeighborOccludes(1,1));assert.equal(tissueNeighborOccludes(1,0),false);const separated=matrix(parts[1]);separated.elements[12]+=.02;assert.equal(tissuePosesCoincide(matrix(parts[0]),parts[0].geometryOrigin,separated,parts[1].geometryOrigin),false);
});
