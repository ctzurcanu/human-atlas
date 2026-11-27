import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:`
 export * from './app/liver-surface-features';
 export * from './app/liver-landmark';
 export * from './app/liver-reference-selection';
 export * from './app/anatomy-label-plan';
 export * from './app/hierarchy-navigation';
 export {hierarchyEntries} from './app/anatomy-hierarchy';
 export {buildSystemTree} from './app/system-tree';
 export * from './app/preferred-source-representations';
 export * from './app/reviewed-catalogue-identities';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const path='.local-models/geometry-workbench/liver-completion/liver-surface-controls-v143/with-landmark.json';
const source=fs.existsSync(path)?JSON.parse(fs.readFileSync(path)):null;
const skip={skip:!source&&'Local source data is not redistributed.'};
const walk=n=>[n,...(n.children??[]).flatMap(walk)];
const code=n=>n.ta98Term??n.terminology?.ta98;

test('18 surface annotations retain the carrier and reject stale or non-surface identities',skip,()=>{
 const registry=api.liverFeatureRegistry(source);assert.equal(registry.termGroups.length,18);
 const ordinary=JSON.parse(fs.readFileSync('.local-models/geometry-workbench/liver-completion/male-liver-section-install-v142/atlas.json'));
 assert.deepEqual(source.parts,ordinary.parts,'surface annotations never add duplicate anatomy');
 assert.deepEqual(source.chunks,ordinary.chunks);
 for(const mutate of [a=>a.liverFeatureIntegration.bankSha256='0'.repeat(64),a=>a.liverFeatureIntegration.termGroups[0].surface.faces.push(9999999),a=>a.liverFeatureIntegration.termGroups[0].representation='finite-tissue']){
  const wrong=structuredClone(source);mutate(wrong);assert.equal(api.liverFeatureRegistry(wrong),null);
 }
});

test('surface selections use ordinary labels with the actual facet anchor',skip,()=>{
 for(const group of source.liverFeatureIntegration.termGroups){
  const c=source.concepts.find(c=>c.id===group.conceptId);
  const item=api.selectedAnatomyLabels(source,c.elements,c,false);
  assert.equal(item.length,1);assert.equal(item[0].id,c.id);assert.deepEqual(item[0].anchor.point,group.anchor.point);
  assert.equal(api.anatomyLabelChoice(item[0]).id,c.id);
 }
});

test('a small landmark uses the standard callout without claiming a completed opening',skip,()=>{
 const r=api.liverLandmarkRegistry(source);assert.ok(r);
 const group=r.marker.termGroups[0],choice=source.concepts.find(c=>c.id===group.conceptId);
 const labels=api.selectedAnatomyLabels(source,choice.elements,choice,false);
 assert.equal(labels.length,1);assert.equal(labels[0].name,'Porta hepatis');assert.deepEqual(labels[0].anchor.point,group.anchor.point);
 assert.equal(group.completeOpeningRepresented,false);
 const wrong=structuredClone(source);wrong.liverLandmarkIntegration.termGroups[0].completeOpeningRepresented=true;assert.equal(api.liverLandmarkRegistry(wrong),undefined);
});

test('reference checkmarks toggle annotations while retaining the same anatomical host',skip,()=>{
 const references=api.liverReferenceChoices(source);
 for(const r of references.values()){
  const shown=api.toggleLiverReferenceChoice(r,undefined);assert.equal(shown.id,r.choice.id);
  const cleared=api.toggleLiverReferenceChoice(r,shown);assert.equal(cleared.id,r.host.id);
  assert.deepEqual(cleared.elements,shown.elements,'clearing an annotation never hides its liver host');
  assert.equal(api.liverReferenceActive(r,cleared),false);
 }
});

test('Systems exposes all 69 liver terms, distinguishing mapped controls from empty terms',skip,()=>{
 const atlas=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(source)));
 const nodes=walk(api.hierarchyNavigation(atlas,'systems'));
 const liver=nodes.find(n=>n.name==='Liver'&&code(n)==='A05.8.01.001');assert.ok(liver);
 const terms=new Set(walk(liver).map(code).filter(c=>c?.startsWith('A05.8.01.')));assert.equal(terms.size,69);
 const flatNode=n=>[n,...(n.nodes??[]).flatMap(flatNode)];
 const visibleTree=api.buildSystemTree(api.hierarchyEntries(atlas),'systems',atlas).flatMap(b=>b.nodes.flatMap(flatNode));
 const visibleTerms=new Set(visibleTree.map(n=>n.entry?.ta98Term??n.terminology?.ta98).filter(c=>c?.startsWith('A05.8.01.')));
 assert.equal(visibleTerms.size,69,'the actual rendered Systems tree uses the same complete inventory');
 for(const group of source.liverFeatureIntegration.termGroups)assert.deepEqual(walk(liver).find(n=>code(n)===group.term).elements,group.partIds);
 assert.deepEqual(walk(liver).find(n=>code(n)==='A05.8.01.056').elements,[],'missing microscopic anatomy never borrows the whole liver');
});
