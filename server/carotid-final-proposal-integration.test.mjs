import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
const base='.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1';
const future=base+'/isolated-future-relations-v3';
const futureParentsPath=future+'/app/data/ta98-parents.json';
const futureParents=fs.existsSync(futureParentsPath)?fs.readFileSync(futureParentsPath):undefined;
const sha=bytes=>createHash('sha256').update(bytes).digest('hex');
if(futureParents)assert.equal(sha(futureParents),'84799edca98ac140b22989c2d583e52f5edda0ae51954a4b1295b5f48a04bf20');
const built=await build({plugins:futureParents?[{name:'exact-future-ta98-parents',setup(builder){builder.onLoad({filter:/app\/data\/ta98-parents\.json$/},()=>({contents:futureParents.toString(),loader:'json'}));}}]:[],stdin:{contents:`
 export * from './app/lingual-arterial-relations';
 export * from './app/ta98-modeled-parent';
 export * from './app/anatomical-components';
 export * from './app/anatomical-relations';
 export * from './app/carotid-surface-features';
 export * from './app/anatomical-reference-selection';
 export * from './app/anatomy-label-plan';
 export * from './app/hierarchy-navigation';
 export {resolveGuestHierarchy} from './app/guest-hierarchy';
 export {hierarchyEntries} from './app/anatomy-hierarchy';
 export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities';
 export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const proposalPath=base+'/bilateral-body-catalogue-v1/proposed.json';
const expected='346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6';
const available=fs.existsSync(proposalPath),skip={skip:!available&&'Private local proposal is not redistributed.'};
const raw=available?fs.readFileSync(proposalPath):null;
const proposal=raw?JSON.parse(raw):null;
const prepared=proposal?api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(proposal))):null;
const bodyId='TA98BODY:male:carotid-body:left:macro';
const hostId='TA98PART:ZA:Internal carotid artery.l:A12.2.06.002';
const referenceId='TA98REFERENCE:male:carotid-sinus:left:medial-proximal-ICA';
const leftCCA='ZA:Left common carotid artery';
const walk=n=>[n,...(n.children??[]).flatMap(walk)];
const code=n=>n.ta98Term??n.terminology?.ta98;
const localAsset=url=>'.local-models/'+url.split('/local-models/')[1];
let systems;
const systemTree=()=>systems??=api.hierarchyNavigation(prepared,'systems');

test('the exact proposal registers both approved common banks and preserves the installed left records',skip,()=>{
 assert.equal(sha(raw),expected);
 const installed=JSON.parse(fs.readFileSync(base+'/catalogue-adapter-v1/proposed-v3.json'));
 assert.equal(sha(fs.readFileSync(base+'/catalogue-adapter-v1/proposed-v3.json')),'a9ae0fbaaf62fdf600a966e43641f54a4b20621944784e31ec43f43260428638');
 for(const id of installed.carotidFiniteWallIntegration.parts){assert.deepEqual(proposal.parts.find(p=>p.id===id),installed.parts.find(p=>p.id===id),'Left installed physical record changed: '+id);}
 const banks={left:'e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36',right:'7b7bcbf93627134a3e240bb10716c202a1455e4779a4f429ae309a66145eb737'};
 for(const side of ['left','right']){
  const family=api.registeredLingualFamily(prepared.parts,side);assert.equal(family?.size,6,side);
  for(const p of family.values()){assert.equal(api.isInstalledLingualArterialPart(p),true,p.id);assert.equal(api.nativePartFrame(p),'ZA:',p.id);assert.equal(p.carotidFiniteWallBankSha256,banks[side]);}
  const bif=prepared.parts.find(p=>p.id===`TA98ARTERIAL:male:carotid-bifurcation:${side}:finite-wall`);assert.ok(bif);
  assert.equal(api.isInstalledCarotidBifurcationPart(bif),true);assert.equal(api.nativePartFrame(bif),'ZA:');assert.deepEqual(bif.geometryOrigin,[0,1.625,0]);
  assert.equal(bif.carotidFiniteWallRevision,`${side}-four-cell-common-bank-v${side==='left'?2:1}`);
  assert.equal(sha(fs.readFileSync(localAsset(prepared.chunks[bif.chunk].url))),banks[side]);
 }
});
test('the single left body is authored as a chemoreceptor organ in the common-carotid Systems ancestry',skip,()=>{
 const body=prepared.parts.find(p=>p.id===bodyId),concept=prepared.concepts.find(c=>c.id===bodyId);
 assert.ok(body);assert.equal(proposal.parts.filter(p=>p.ta98Term==='A12.2.04.007').length,1);
 assert.equal(body.tissue,'organ');assert.equal(body.anatomicalEntityType,'chemoreceptor-organ');assert.match(body.description,/chemoreceptor organ/);
 assert.match(body.displaySystemQualification,/not an artery/);assert.equal(body.ta98Parent,'A12.2.04.006');assert.equal(concept.partOf,leftCCA);
 assert.equal(api.nativePartFrame(body),'ZA:');
 const choices=walk(systemTree()).filter(n=>n.id===bodyId);assert.equal(choices.length,1);assert.deepEqual(choices[0].elements,[bodyId]);
 assert.ok(api.hierarchyAncestors(systemTree(),choices[0]).some(n=>code(n)==='A12.2.04.006'));
 const entries=api.hierarchyEntries(prepared).filter(e=>e.id===bodyId);assert.equal(entries.length,1);assert.equal(entries[0].ta98Term,'A12.2.04.007');
});
test('the exact body gets byte-RGB finish and altered or unregistered sources are rejected',skip,()=>{
 const body=proposal.parts.find(p=>p.id===bodyId);
 assert.deepEqual(api.carotidBodySurfaceFinish(body),{color:[165.75,107.1,71.4],roughness:.65,metalness:0});
 for(const mutate of [p=>p.geometryOrigin=[0,0,0],p=>p.sourceId='Unregistered analytic body',p=>p.positions=4,p=>p.normalSigns++,p=>p.indexCount+=3,p=>p.provenance.url='https://example.com/donor',p=>p.sectionAssembly='carotid-body:male:right']){
  const changed=structuredClone(body);mutate(changed);assert.equal(api.nativePartFrame(changed),undefined,mutate.toString());assert.equal(api.carotidBodySurfaceFinish(changed),undefined,mutate.toString());
 }
 for(const mutate of [p=>p.id='ZA:Unregistered carotid body',p=>p.material='unregistered-unit-rgb-material']){const changed=structuredClone(body);mutate(changed);assert.equal(api.carotidBodySurfaceFinish(changed),undefined,mutate.toString());}
});
test('Left CCA Show components includes the body and bifurcation while excluding the adjoining ICA carrier',skip,()=>{
 const parts=new Map(prepared.parts.map(p=>[p.id,p])),choice=prepared.concepts.find(c=>c.id===leftCCA);
 const relations=choice.elements.flatMap(id=>api.anatomicalRelations(parts.get(id),prepared.parts));
 const ids=api.anatomicalComponentIds(choice,relations,parts,prepared.concepts);
 assert.ok(ids.includes(bodyId),'Left carotid body is absent from Show components');
 assert.ok(ids.includes('TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall'),'Left bifurcation is absent from Show components');
 // The qualified sinus patch is separately selectable; its whole ICA carrier
 // is an adjoining artery and must not become contained CCA tissue.
 assert.equal(ids.includes(hostId),false,'The entire adjoining ICA carrier became a physical CCA component');
 assert.equal(new Set(ids).size,ids.length);
 const ref=api.anatomicalReferenceChoices(prepared).get('A12.2.04.008');assert.ok(ref);
 assert.deepEqual(ref.choice.elements,[hostId]);
 assert.equal(api.ta98ModeledParent(ref.choice,prepared.concepts,parts)?.id,leftCCA);
 const componentRelations=api.anatomicalComponentRelations(choice,relations,parts,prepared.concepts);
 for(const id of [bodyId,'TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall'])assert.ok(componentRelations.some(r=>r.kind==='contains'&&r.target.id===id),id);
 assert.equal(componentRelations.some(r=>r.kind==='contains'&&r.target.id===hostId),false);
});
test('the sinus remains surface-only with qualified reference choices and the actual label anchor',skip,()=>{
 const registry=api.carotidFeatureRegistry(prepared);assert.ok(registry);const feature=registry.termGroups[0];assert.equal(feature.referenceOnly,true);assert.equal(feature.representation,'surface-region');assert.equal(feature.boundaryProvisional,true);
 assert.equal(proposal.parts.filter(p=>p.id===referenceId).length,0,'A reference acquired a duplicate anatomical mesh');
 const oldIds=new Set(JSON.parse(fs.readFileSync(base+'/catalogue-adapter-v1/proposed-v3.json')).parts.map(p=>p.id));
 assert.deepEqual(proposal.parts.filter(p=>!oldIds.has(p.id)).map(p=>p.id).sort(),[bodyId,'TA98ARTERIAL:male:carotid-bifurcation:right:finite-wall'].sort());
 const ref=api.anatomicalReferenceChoices(prepared).get('A12.2.04.008');assert.ok(ref);assert.equal(ref.choice.id,referenceId);assert.deepEqual(ref.choice.elements,[hostId]);assert.match(ref.choice.name,/medial wall reference/);assert.match(ref.choice.description,/Qualified reference/);assert.match(ref.choice.description,/remain unverified/);
 const labels=api.selectedAnatomyLabels(prepared,[hostId],ref.choice,false);assert.equal(labels.length,1);assert.equal(labels[0].id,referenceId);assert.deepEqual(labels[0].anchor,{partId:hostId,point:feature.anchor.point});assert.equal(api.anatomyLabelChoice(labels[0]).id,referenceId);
 const shown=api.toggleAnatomicalReferenceChoice(ref,undefined),cleared=api.toggleAnatomicalReferenceChoice(ref,shown);assert.equal(shown.id,referenceId);assert.equal(cleared.id,hostId);assert.deepEqual(shown.elements,cleared.elements);
 const choices=walk(systemTree()).filter(n=>n.id===referenceId);assert.equal(choices.length,1);assert.ok(api.hierarchyAncestors(systemTree(),choices[0]).some(n=>code(n)==='A12.2.04.006'));assert.deepEqual(choices[0].elements,[hostId]);
 const host=prepared.parts.find(p=>p.id===hostId),bytes=fs.readFileSync(localAsset(prepared.chunks[host.chunk].url));const g=feature.anchor;
 const triangle=[0,1,2].map(k=>{const i=bytes.readUInt32LE(host.indices+g.face*12+k*4);return [0,1,2].map(axis=>bytes.readFloatLE(host.positions+i*12+axis*4)+host.geometryOrigin[axis]);});
 const replay=[0,1,2].map(axis=>triangle.reduce((sum,p,i)=>sum+p[axis]*g.barycentric[i],0));assert.ok(Math.hypot(...replay.map((v,i)=>v-g.point[i]))<1e-12);
});
test('the exact future TA98 hierarchy retains the .006 parent and a non-solid .008 reference',skip,()=>{
 const coverage=JSON.parse(fs.readFileSync(future+'/reports/ta98/coverage.json'));assert.equal(coverage.atlasSources['local-ta98'].sha256,expected,'Future TA98 asset was built from another male proposal');
 const guest=JSON.parse(fs.readFileSync(future+'/public/assets/ta98.json'));const root=api.hierarchyNavigation(prepared,'guest:ta98',guest),nodes=walk(root);
 const body=nodes.find(n=>n.id==='hierarchy:guest:ta98:A12.2.04.007'),sinus=nodes.find(n=>n.id==='hierarchy:guest:ta98:A12.2.04.008');assert.ok(body);assert.ok(sinus);assert.deepEqual(body.elements,[bodyId]);assert.deepEqual(sinus.elements,[hostId]);
 for(const child of [body,sinus]){const ancestors=api.hierarchyAncestors(root,child);assert.equal(api.ta98PartOfChoice(child,ancestors)?.terminology?.ta98,'A12.2.04.006');}
 // The qualified app representation remains an existing wall-surface
 // reference. The source annotation does not certify bulb morphology.
 const sourceSinus=guest.nodes.find(n=>n.id==='A12.2.04.008');assert.ok(sourceSinus);
 assert.equal(sourceSinus.name,'carotid sinus');assert.match(sourceSinus.description,/Latin: sinus caroticus/);assert.match(sourceSinus.description,/Source entity type: immaterial entity/);
 const selected=api.selectedCarotidSurfaceFeature(prepared,sinus);assert.equal(selected?.referenceOnly,true);assert.equal(selected?.representation,'surface-region');
 const reference=api.anatomicalReferenceChoices(prepared).get('A12.2.04.008');assert.match(reference.choice.description,/Qualified reference/);assert.match(reference.choice.description,/remain unverified/);
 const resolvedSinus=api.resolveGuestHierarchy(prepared,guest).flatMap(walk).find(n=>n.id==='A12.2.04.008');assert.ok(resolvedSinus);
 assert.match(resolvedSinus.description,/Qualified reference/);assert.match(resolvedSinus.description,/remain unverified/);assert.doesNotMatch(resolvedSinus.description,/supplemental duct/);
 assert.deepEqual(api.carotidSurfaceLabelReference(prepared,sinus),api.carotidSurfaceLabelReference(prepared,prepared.concepts.find(c=>c.id===referenceId)));
 const cervicalSinus=nodes.find(n=>n.id==='hierarchy:guest:ta98:A12.2.06.003');assert.ok(cervicalSinus);
 assert.deepEqual(cervicalSinus.elements,[hostId]);
 assert.equal(api.ta98PartOfChoice(cervicalSinus,api.hierarchyAncestors(root,cervicalSinus))?.terminology?.ta98,'A12.2.06.002');
 assert.equal(api.selectedCarotidSurfaceFeature(prepared,cervicalSinus),selected);
 assert.deepEqual(api.carotidSurfaceLabelReference(prepared,cervicalSinus),api.carotidSurfaceLabelReference(prepared,sinus));
 const cervicalReference=api.anatomicalReferenceChoices(prepared).get('A12.2.06.003');assert.ok(cervicalReference);
 assert.equal(cervicalReference.choice.ta98Term,'A12.2.06.003');assert.equal(cervicalReference.choice.partOf,hostId);
 assert.equal(api.ta98ModeledParent(cervicalReference.choice,prepared.concepts,new Map(prepared.parts.map(p=>[p.id,p])))?.id,hostId);
 const cervicalResolved=api.resolveGuestHierarchy(prepared,guest).flatMap(walk).find(n=>n.id==='A12.2.06.003');assert.ok(cervicalResolved);
 assert.match(cervicalResolved.description,/Qualified reference/);assert.match(cervicalResolved.description,/remain unverified/);
 const cervicalLabels=api.selectedAnatomyLabels(prepared,[hostId],cervicalReference.choice,false);
 assert.equal(cervicalLabels.length,1);assert.deepEqual(cervicalLabels[0].anchor,{partId:hostId,point:selected.anchor.point});
 assert.equal(walk(systemTree()).filter(n=>n.id===referenceId).length,1,'Authority occurrences duplicated a Systems entry');
 assert.equal(sha(fs.readFileSync(proposalPath)),expected,'Tests mutated the proposed catalogue');
});
