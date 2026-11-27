import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {build} from './atlas-test-build.mjs';

const built=await build({stdin:{contents:"export * from './app/viewer-state';export * from './app/anatomical-components';export * from './app/anatomy-label-plan';export * from './app/eye-membership';export {requestedSelectionChoice} from './app/requested-selection-choice';export {resolveGuestHierarchy} from './app/guest-hierarchy';export * from './app/ta98-term-bindings';export * from './app/local-ta98-layer-catalogue';export * from './app/laterality';export {buildSystemTree} from './app/system-tree';export {hierarchyEntries} from './app/anatomy-hierarchy';export {ocularMouthRelations,withVerifiedAngularIovMouth} from './app/ocular-mouth-relations';",resolveDir:process.cwd()},bundle:true,write:false,platform:'node',format:'esm'});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const route=api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true);
const file=process.env.ATLAS_FAMILY_CANDIDATE??'.local-models'+route.slice('/local-models'.length);
const a=JSON.parse(fs.readFileSync(file));
const parts=new Map(a.parts.map(p=>[p.id,p])),concepts=new Map(a.concepts.map(c=>[c.id,c]));
const membership=a.eyeVesselMembership;
const physical=a.eyePhysicalComponentMembership;
const physicalIds=eye=>physical?.groups.find(g=>g.eyeConceptId===eye.id)?.activeComponentIds??eye.elements;

test('Verified angular connection is added only for the installed bank and remains reciprocal',()=>{
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 const ids=[prefix+'right_angular_vein_ID',prefix+'right_inferior_ophthalmic_vein_ID'];
 const bank='d65fa6eeaca3bad2d54a2b8988a559a7b8e9343c0d07e99902fdc0530f492482';
 const source={parts:ids.map(id=>({id,system:'venous',chunk:0})),chunks:[{sha256:bank}]};
 const result=api.withVerifiedAngularIovMouth(source);
 assert.equal(api.ocularMouthRelations(result.parts[0],result.parts)[0].target.id,ids[1]);
 assert.equal(api.ocularMouthRelations(result.parts[1],result.parts)[0].target.id,ids[0]);
 assert.equal(source.parts[0].ocularMouthConnections,undefined);
 assert.equal(api.withVerifiedAngularIovMouth({...source,chunks:[{sha256:'other'}]}).parts[0].ocularMouthConnections,undefined);
});

test('The installed review catalogue exposes the verified angular connection',()=>{
 const atlas=api.withVerifiedAngularIovMouth(a);
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 const angular=atlas.parts.find(part=>part.id===prefix+'right_angular_vein_ID');
 const iov=atlas.parts.find(part=>part.id===prefix+'right_inferior_ophthalmic_vein_ID');
 assert.ok(angular&&iov);
 assert.ok(api.ocularMouthRelations(angular,atlas.parts).some(relation=>relation.kind==='continuous'&&relation.target.id===iov.id));
 assert.ok(api.ocularMouthRelations(iov,atlas.parts).some(relation=>relation.kind==='continuous'&&relation.target.id===angular.id));
});

test('Vorticose openings are reciprocal and depend on the verified geometry banks',()=>{
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 const atlas=api.withVerifiedAngularIovMouth(a);
 for(const side of ['left','right'])for(const role of ['superior','inferior']){
  const vortex=atlas.parts.find(part=>part.id===prefix+side+'_vorticose_veins_ID');
  const vein=atlas.parts.find(part=>part.id===prefix+side+'_'+role+'_ophthalmic_vein_ID');
  assert.ok(api.ocularMouthRelations(vortex,atlas.parts).some(relation=>relation.target.id===vein.id));
  assert.ok(api.ocularMouthRelations(vein,atlas.parts).some(relation=>relation.target.id===vortex.id));
  const altered={...a,chunks:a.chunks.map((chunk,index)=>index===vortex.chunk?{...chunk,sha256:'changed'}:chunk)};
  const fresh=api.withVerifiedAngularIovMouth(altered);
  assert.ok(!api.ocularMouthRelations(fresh.parts.find(part=>part.id===vortex.id),fresh.parts).some(relation=>relation.target.id===vein.id));
 }
});

test('Verified ophthalmic connections are reciprocal on both sides',()=>{
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 const atlas=api.withVerifiedAngularIovMouth(a);
 for(const side of ['left','right'])for(const [vessel,peer] of [
  ['superior_ophthalmic_vein_ID','cavernous_sinus_ID'],
  ['inferior_ophthalmic_vein_ID','cavernous_sinus_ID'],
  ['superior_ophthalmic_vein_ID','central_retinal_vein_ID:region:A12.3.06.112'],
  ['inferior_ophthalmic_vein_ID','angular_vein_ID'],
  ['inferior_ophthalmic_vein_ID','pterygoid_plexus_ID'],
 ]){
  const first=atlas.parts.find(part=>part.id===prefix+side+'_'+vessel);
  const second=atlas.parts.find(part=>part.id===prefix+side+'_'+peer);
  assert.ok(api.ocularMouthRelations(first,atlas.parts).some(r=>r.target.id===second.id));
  assert.ok(api.ocularMouthRelations(second,atlas.parts).some(r=>r.target.id===first.id));
 }
});

test('Ordinary male catalogue exposes the same verified eye-vessel connections',()=>{
 const route=api.localTa98LayerCatalogue('local-ta98',new URLSearchParams(),true);
 const file='.local-models'+route.slice('/local-models'.length);
 const ordinary=api.withVerifiedAngularIovMouth(JSON.parse(fs.readFileSync(file)));
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 for(const side of ['left','right'])for(const role of ['superior','inferior']){
  const vortex=ordinary.parts.find(part=>part.id===prefix+side+'_vorticose_veins_ID');
  const vein=ordinary.parts.find(part=>part.id===prefix+side+'_'+role+'_ophthalmic_vein_ID');
  assert.ok(vortex&&vein);
  assert.ok(api.ocularMouthRelations(vortex,ordinary.parts).some(r=>r.target.id===vein.id));
  assert.ok(api.ocularMouthRelations(vein,ordinary.parts).some(r=>r.target.id===vortex.id));
 }
});

test('Infraorbital and adjacent vein openings are reciprocal on both sides',()=>{
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 const atlas=api.withVerifiedAngularIovMouth(a);
 for(const side of ['left','right'])for(const role of ['angular_vein','anterior_superior_alveolar_vein','middle_superior_alveolar_vein']){
  const vein=atlas.parts.find(part=>part.id===prefix+side+'_infraorbital_vein_ID');
  const peer=atlas.parts.find(part=>part.id===prefix+side+'_'+role+'_ID');
  assert.ok(api.ocularMouthRelations(vein,atlas.parts).some(r=>r.target.id===peer.id));
  assert.ok(api.ocularMouthRelations(peer,atlas.parts).some(r=>r.target.id===vein.id));
 }
});
const legacyEyeCodes='3f.3g.8c.8d.ii.ij.1ai.1aj.1fu.1fv.1gl.1gm.1gp.1gq.1mf.1mg.22c.22d.272.273.2jz.2k0.2qv.2qw.2uq.2ur.37t.37u.39o.39p.3mc.3md.3ml.3mm';
const legacyEyeIds=legacyEyeCodes.split('.').map(code=>a.parts[parseInt(code,36)].id);
const eyeBloodIds=physical.groups.flatMap(group=>group.ocularVesselComponentIds);

test('legacy generic Systems Eye retains all native adnexa and includes its registered circulation',()=>{
 const nervous=api.buildSystemTree(api.hierarchyEntries(a),'systems',a).find(branch=>branch.id==='nervous');
 const find=nodes=>{for(const node of nodes){if(node.name==='Eye'&&node.parts.some(part=>eyeBloodIds.includes(part.id)))return node;if(node.kind==='group'){const result=find(node.nodes);if(result)return result;}}};
 const eye=find(nervous.nodes);assert.ok(eye);assert.equal(eye.parts.length,58);
 assert.ok(legacyEyeIds.every(id=>eye.parts.some(part=>part.id===id)));
 assert.ok(eyeBloodIds.every(id=>eye.parts.some(part=>part.id===id)));
 assert.ok(!eye.parts.some(part=>part.id.includes('pterygoid_plexus')));
 assert.equal(new Set(eye.parts.map(part=>part.id)).size,58);
});

test('the actual saved 34-native Eye URL restores its 24 vessels without moving the camera',()=>{
 const frame=[.892865,.460691,.806777,.592623,.934391,.799421,.645823,.450056,.448219,.5,1,.5];
 const base={visible:['sensory','arterial','venous'],selected:[],isolate:false,view:'front',rotate:false,explode:0,reset:4,focus:7};
 const restored=api.readViewUrl('?model=local-ta98&p=6rv.1i5aj2l.9t7hdf&s='+legacyEyeCodes+'&context=0&skin=0&frame='+frame.join(','),a,base);
 assert.equal(restored.selected.length,58);assert.ok(legacyEyeIds.every(id=>restored.selected.includes(id)));assert.ok(eyeBloodIds.every(id=>restored.selected.includes(id)));
 assert.deepEqual(restored.camera,frame);assert.equal(restored.reset,base.reset);assert.equal(restored.focus,0,'The saved URL requests no Center action');
 assert.deepEqual(api.restoreEyeAssemblySelection(a,legacyEyeIds.slice(0,3)),legacyEyeIds.slice(0,3));
 const choice=api.restoredEyeAssemblyChoice(a,restored.selected);assert.equal(choice.name,'Eye');
 assert.match(choice.coverageLimitation,/extraocular/);assert.equal(choice.ta98Term,undefined);
 const labels=api.selectedAnatomyLabels(a,restored.selected,choice,false);assert.equal(labels.length,1);assert.equal(labels[0].members.length,58);assert.equal(labels[0].components.length,58);
 assert.equal(api.anatomyLabelChoice(labels[0]).elements.length,58);
 assert.equal(api.anatomicalComponentIds(choice,[],parts,a.concepts).length,58);
});

test('generic Eye and TA parent use explicit ocular ownership without crossing sides or adding recipient plexus',()=>{
 for(const term of ['Eye','A15.2.00.001','TA98:A15.2.00.001']){
  const selected=api.resolveSelection(a,[term]);assert.equal(selected.length,58);assert.ok(eyeBloodIds.every(id=>selected.includes(id)));
  assert.equal(api.requestedSelectionChoice(a,[term],selected).componentElements.length,58);
 }
 const group=physical.groups.find(group=>group.side==='left'),native=group.activeComponentIds.filter(id=>!eyeBloodIds.includes(id)).map(id=>parts.get(id));
 const left=api.eyePresentationParts(a,native,{name:'Eye'});assert.equal(left.length,18);assert.ok(left.every(part=>!part.id.includes('right_')&&!/\.r$/.test(part.id)));
 const hierarchy={schema:'human-atlas-hierarchy/v2',id:'ta98',name:'Terminologia Anatomica',roots:['A15.2.00.001'],nodes:[{id:'A15.2.00.001',name:'Eye',partIds:legacyEyeIds}]};
 const [node]=api.resolveGuestHierarchy(a,hierarchy);assert.equal(node.parts.length,58);assert.ok(eyeBloodIds.every(id=>node.parts.some(part=>part.id===id)));
});

test('ocular arterial and venous continuity requires reciprocal same-side verified mouth metadata',()=>{
 const prefix='BDM:male:human_20_male_cardiovascular_system-';
 const source=structuredClone(parts.get(prefix+'left_short_posterior_ciliary_arteries_ID'));
 const parent=structuredClone(parts.get(prefix+'left_ophthalmic_artery_ID'));
 const opposite=structuredClone(parts.get(prefix+'right_ophthalmic_artery_ID'));
 const marker=targetId=>({targetId,contractSha256:'a'.repeat(64),boundaryVertices:6,opening:true});
 source.ocularMouthConnections=[marker(parent.id)];parent.ocularMouthConnections=[];
 assert.deepEqual(api.ocularMouthRelations(source,[source,parent]),[]);
 parent.ocularMouthConnections=[marker(source.id)];
 assert.deepEqual(api.ocularMouthRelations(source,[source,parent]).map(row=>[row.kind,row.target.id]),[['continuous',parent.id]]);
 assert.equal(api.ocularMouthRelations(parent,[source,parent])[0].target.id,source.id);
 parent.ocularMouthConnections[0].contractSha256='b'.repeat(64);
 assert.deepEqual(api.ocularMouthRelations(source,[source,parent]),[]);
 source.ocularMouthConnections=[marker(opposite.id)];opposite.ocularMouthConnections=[marker(source.id)];
 assert.deepEqual(api.ocularMouthRelations(source,[source,opposite]),[]);
 const vein=structuredClone(parts.get(prefix+'left_inferior_ophthalmic_vein_ID'));
 const recipient={id:prefix+'left_pterygoid_plexus_ID',system:'venous'};
 vein.ocularMouthConnections=[marker(recipient.id)];recipient.ocularMouthConnections=[marker(vein.id)];
 assert.deepEqual(api.ocularMouthRelations(vein,[vein,recipient]).map(row=>[row.kind,row.target.id]),[['continuous',recipient.id]]);
 assert.equal(api.ocularMouthRelations(recipient,[vein,recipient])[0].target.id,vein.id);
 const wrongSide={...recipient,id:prefix+'right_pterygoid_plexus_ID'};
 vein.ocularMouthConnections=[marker(wrongSide.id)];
 assert.deepEqual(api.ocularMouthRelations(vein,[vein,wrongSide]),[]);
 vein.ocularMouthConnections=[marker(recipient.id)];recipient.system='arterial';
 assert.deepEqual(api.ocularMouthRelations(vein,[vein,recipient]),[]);
 recipient.system='venous';recipient.ocularMouthConnections[0].opening=false;
 assert.deepEqual(api.ocularMouthRelations(vein,[vein,recipient]),[]);
});

test('an eye selects and expands its own vortex mesh and supplies its anatomical parent',()=>{
 assert.equal(membership.groups.length,2);
 for(const g of membership.groups){
  const eye=concepts.get(g.eyeConceptId),vein=concepts.get(g.vesselConceptId);
  const opposite=membership.groups.find(n=>n.side!==g.side).vesselConceptId;
  const selected=api.resolveSelection(a,[eye.id]);assert.deepEqual(new Set(selected),new Set(physicalIds(eye)));assert.ok(selected.includes(g.vesselConceptId));assert.ok(selected.includes(g.referenceScleraId));assert.ok(!selected.includes(opposite));
  const components=api.anatomicalComponentIds(eye,[],parts,a.concepts);assert.ok(components.includes(g.vesselConceptId));assert.ok(!components.includes(opposite));
  const directEye=physical?.groups.find(group=>group.eyeConceptId===eye.id)?.canonicalEyeConceptId??eye.id;
  assert.equal(vein.partOf,directEye);assert.equal(vein.sourceParentConceptId,g.sourceParentConceptId);
  assert.equal(api.storedAssemblyParent(vein,a.concepts,parts).id,directEye);
  assert.deepEqual(api.resolveSelection(a,[vein.id]),[g.vesselConceptId]);
  const labels=api.selectedAnatomyLabels(a,selected,eye,false);assert.equal(labels.length,1);assert.ok(labels[0].components.includes(g.vesselConceptId));
 }
 const both=api.resolveSelection(a,[membership.wholeEyeConceptId]);assert.deepEqual(new Set(both),new Set(membership.groups.flatMap(g=>physicalIds(concepts.get(g.eyeConceptId)))));
 for(const g of membership.groups)assert.ok(both.includes(g.vesselConceptId));
});

test('eye ownership preserves every physical mesh and every circulatory TA98 binding',()=>{
 const history=a.bloodEyeMembershipHistory.at(-1),base='.local-models/geometry-workbench/biodigital-male-circulation-v1';
 const release=path.join(base,history.revision),installedMembership=JSON.parse(fs.readFileSync(path.join(release,'circulation-atlas.json')));
 const prep=JSON.parse(fs.readFileSync(path.join(release,'prepare-report.json')));
 const prior=JSON.parse(fs.readFileSync(prep.sourceCatalogue));
 for(const field of ['parts','chunks','materials','defaultHidden','ta98TermBindings','sourceCollectionIdentities','sourceVesselFamilies'])assert.deepEqual(installedMembership[field],prior[field],field);
 const oldTerms=api.catalogueTa98TermParts(prior),newTerms=api.catalogueTa98TermParts(a);
 assert.deepEqual([...newTerms.keys()].sort(),[...oldTerms.keys()].sort());
 // Current geometry may be repaired after the metadata-only membership
 // release.  TA98 membership is the set of physical IDs, not buffer offsets,
 // normals or bounds; the historical release's unchanged meshes are checked
 // independently above.
 for(const [term,members] of oldTerms)assert.deepEqual(newTerms.get(term)?.map(p=>p.id),members.map(p=>p.id),term);
 assert.equal(membership.geometryChanged,false);assert.equal(membership.placementAcceptance,false);
});

test('active physical eye includes its retina and records unfit unique source anatomy as hidden guidance',{skip:!physical},()=>{
 const history=a.bloodEyePhysicalComponentsHistory.at(-1),base='.local-models/geometry-workbench/biodigital-male-circulation-v1';
 const prep=JSON.parse(fs.readFileSync(path.join(base,history.revision,'prepare-report.json')));
 const prior=JSON.parse(fs.readFileSync(prep.sourceCatalogue)),oldConcepts=new Map(prior.concepts.map(c=>[c.id,c]));
 const hidden=new Set(a.defaultHidden);
 for(const group of physical.groups){
  const eye=concepts.get(group.eyeConceptId),before=oldConcepts.get(eye.id),retina=parts.get(group.retinaPartId);
  assert.ok(retina&&retina.id===`ZA:Retina.${group.side==='left'?'l':'r'}`&&!retina.suppressed&&!hidden.has(retina.id));
  const expected=new Set(before.elements.filter(id=>!group.guidancePartIds.includes(id)).concat(retina.id,group.ocularVesselComponentIds??[]));
  assert.deepEqual(new Set(eye.elements),expected);assert.deepEqual(eye.componentElements,eye.elements);
  assert.deepEqual(new Set(group.activeComponentIds),expected);
  assert.ok(eye.elements.every(id=>parts.has(id)&&!parts.get(id).suppressed&&!hidden.has(id)));
  assert.equal(concepts.get(retina.conceptId).partOf,eye.id);
  assert.equal(api.storedAssemblyParent(concepts.get(retina.conceptId),a.concepts,parts).id,eye.id);
  assert.deepEqual(eye.guidanceElements,group.guidancePartIds);
  for(const proof of group.guidanceEvidence){
   const part=parts.get(proof.partId),sourceConcept=concepts.get(part.conceptId),oldSource=oldConcepts.get(part.conceptId);
   assert.ok(hidden.has(part.id)&&before.elements.includes(part.id)&&!eye.elements.includes(part.id));
   assert.equal(part.reviewStatus,'alignment-and-overlap-pending');assert.equal(part.provenance.label,'BodyParts3D 4.0');
   assert.equal(proof.duplicateClaimed,false);assert.equal(proof.anatomicalConceptRetained,true);
   assert.equal(sourceConcept.partOf,eye.id);assert.equal(sourceConcept.guidanceOnlyGeometry,true);
   for(const field of ['elements','ta98Term','ta98Kind'])assert.deepEqual(sourceConcept[field],oldSource[field],field);
  }
  const selected=api.resolveSelection(a,[eye.id]),expanded=api.selectedAnatomyLabels(a,selected,eye,true);
  assert.deepEqual(new Set(expanded.flatMap(label=>label.members)),expected);
  for(const label of expanded){
   assert.equal(api.displayLaterality(label.name).side,group.side);
   assert.equal(api.lateralityClass(api.displayLaterality(label.name).side),`side-${group.side}`);
   assert.deepEqual(api.anatomyLabelChoice(label).elements,label.members);
  }
 }
 const whole=concepts.get(physical.wholeEyeConceptId);
 assert.deepEqual(new Set(whole.elements),new Set(physical.groups.flatMap(g=>g.activeComponentIds)));
 assert.deepEqual(whole.componentElements,whole.elements);
 assert.equal(physical.wholeEyeAnatomyComplete,false);assert.equal(physical.placementAcceptance,false);
});

test('ocular circulation belongs to each Eye and is usable in both Systems branches',{skip:!a.eyeOcularCirculationMembership},()=>{
 const tree=api.buildSystemTree(api.hierarchyEntries(a),'systems',a);
 const nervous=tree.find(branch=>branch.id==='nervous'),blood=tree.find(branch=>branch.id==='circulatory');
 const find=(nodes,id)=>{for(const node of nodes){if(node.id===id)return node;if(node.kind==='group'){const result=find(node.nodes,id);if(result)return result;}}};
 assert.equal(new Set(nervous.parts.map(part=>part.id)).size,nervous.parts.length,'Eye assembly counts unique physical meshes');
 for(const group of a.eyeOcularCirculationMembership.groups){
  const eye=concepts.get(group.eyeConceptId),selected=api.resolveSelection(a,[eye.id]);
  const canonical=physical?.groups.find(g=>g.eyeConceptId===eye.id)?.canonicalEyeConceptId??eye.id;
  const side=find(nervous.nodes,'eye:whole:'+canonical),vessels=find(nervous.nodes,'eye:circulation:'+canonical);
  assert.ok(side&&vessels,'The Systems Eye node exposes its circulation');
  assert.deepEqual(new Set(side.parts.map(part=>part.id)),new Set(eye.elements));
  assert.deepEqual(new Set(vessels.parts.map(part=>part.id)),new Set(group.componentIds));
  const opposite=a.eyeOcularCirculationMembership.groups.find(other=>other.side!==group.side);
  assert.ok(!selected.some(id=>opposite.componentIds.includes(id)));
  for(const id of group.componentIds){assert.ok(selected.includes(id));assert.ok(blood.parts.some(part=>part.id===id),'Circulatory membership retained');}
  for(const id of group.conceptIds){const vessel=concepts.get(id);assert.equal(vessel.partOf,canonical);assert.ok(vessel.sourceParentConceptId);assert.equal(api.storedAssemblyParent(vessel,a.concepts,parts).id,canonical);}
 }
});

test('literal Eye and eyeball select the same modeled vascular subset with standard labels and parents',{skip:!a.eyeCanonicalConceptMembership},()=>{
 const registry=a.eyeCanonicalConceptMembership;
 const pair=concepts.get(registry.wholeEyeConceptId);
 assert.equal(pair.name,'Eyes');assert.equal(pair.wholeEyeAnatomyComplete,false);
 for(const group of registry.groups){
  const eye=concepts.get(group.eyeConceptId),ball=concepts.get(group.eyeballConceptId);
  assert.equal(eye.name,group.side+' eye');assert.equal(eye.partOf,pair.id);
  assert.equal(ball.partOf,eye.id);assert.equal(api.storedAssemblyParent(ball,a.concepts,parts).id,eye.id);
  assert.equal(eye.wholeEyeAnatomyComplete,false);assert.equal(eye.physicalSubsetOnly,true);
  assert.deepEqual(new Set(eye.elements),new Set(ball.elements));assert.equal(eye.elements.length,18);
  assert.equal(group.vascularComponentIds.length,12);
  for(const choice of [eye,ball]){
   const selected=api.resolveSelection(a,[choice.id]);assert.deepEqual(new Set(selected),new Set(group.activeComponentIds));
   const components=api.anatomicalComponentIds(choice,[],parts,a.concepts);
   assert.ok(group.vascularComponentIds.every(id=>components.includes(id)));
   const labels=api.selectedAnatomyLabels(a,selected,choice,false);assert.equal(labels.length,1);
   assert.equal(labels[0].id,choice.id);assert.ok(group.vascularComponentIds.every(id=>labels[0].components.includes(id)));
   assert.deepEqual(api.anatomyLabelChoice(labels[0]).elements,labels[0].members);
   const expanded=api.selectedAnatomyLabels(a,selected,choice,true);
   assert.deepEqual(new Set(expanded.flatMap(label=>label.members)),new Set(selected));
  }
  assert.deepEqual(eye.guidanceElements,group.guidancePartIds);
  assert.ok(group.guidancePartIds.every(id=>a.defaultHidden.includes(id)&&!eye.elements.includes(id)));
  for(const id of group.vascularComponentIds){
   const vessel=concepts.get(parts.get(id).conceptId);
   let parent=api.storedAssemblyParent(vessel,a.concepts,parts);
   while(parent&&parent.id!==ball.id&&parent.id!==eye.id)parent=api.storedAssemblyParent(parent,a.concepts,parts);
   assert.ok(parent?.id===ball.id||parent?.id===eye.id);
  }
 }
 assert.deepEqual(new Set(pair.elements),new Set(registry.groups.flatMap(g=>g.activeComponentIds)));
});

test('Systems exposes canonical Eye names without counting alias meshes twice',{skip:!a.eyeCanonicalConceptMembership},()=>{
 const tree=api.buildSystemTree(api.hierarchyEntries(a),'systems',a),nervous=tree.find(branch=>branch.id==='nervous');
 const find=(nodes,id)=>{for(const n of nodes){if(n.id===id)return n;if(n.kind==='group'){const v=find(n.nodes,id);if(v)return v;}}};
 assert.equal(new Set(nervous.parts.map(p=>p.id)).size,nervous.parts.length);
 for(const group of a.eyeCanonicalConceptMembership.groups){
  const eye=find(nervous.nodes,'eye:whole:'+group.eyeConceptId),blood=find(nervous.nodes,'eye:circulation:'+group.eyeConceptId);
  assert.equal(eye.name,group.side+' eye');assert.equal(eye.parts.length,18);
  assert.deepEqual(new Set(blood.parts.map(p=>p.id)),new Set(group.vascularComponentIds));
 }
});

test('physical eye metadata revision changes no geometry, hidden policy, blood collections or original ownership history',{skip:!physical},()=>{
 const history=a.bloodEyePhysicalComponentsHistory.at(-1),base='.local-models/geometry-workbench/biodigital-male-circulation-v1';
 const release=path.join(base,history.revision),candidate=JSON.parse(fs.readFileSync(path.join(release,'circulation-atlas.json')));
 const prep=JSON.parse(fs.readFileSync(path.join(release,'prepare-report.json'))),prior=JSON.parse(fs.readFileSync(prep.sourceCatalogue));
 const mutable=new Set(['concepts','eyePhysicalComponentMembership','bloodEyePhysicalComponentsHistory']);
 assert.deepEqual(Object.fromEntries(Object.entries(candidate).filter(([k])=>!mutable.has(k))),Object.fromEntries(Object.entries(prior).filter(([k])=>!mutable.has(k))));
 assert.deepEqual(candidate.bloodEyeMembershipHistory,prior.bloodEyeMembershipHistory);
 const oldConcepts=new Map(prior.concepts.map(c=>[c.id,c]));
 for(const group of candidate.eyePhysicalComponentMembership.groups){
  assert.deepEqual(candidate.concepts.find(c=>c.id===group.vesselConceptId),oldConcepts.get(group.vesselConceptId));
 }
 assert.equal(candidate.eyePhysicalComponentMembership.geometryChanged,false);
 assert.equal(candidate.eyePhysicalComponentMembership.visibilityPolicyChanged,false);
});
