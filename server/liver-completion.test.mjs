import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:`
 export * from './app/liver-completion-controls';
 export * from './app/liver-surface-features';
 export * from './app/liver-reference-selection';
 export * from './app/anatomy-label-plan';
 export * from './app/hierarchy-navigation';
 export {hierarchyEntries} from './app/anatomy-hierarchy';
 export {buildSystemTree} from './app/system-tree';
 export * from './app/preferred-source-representations';
 export * from './app/reviewed-catalogue-identities';
 export * from './app/viewer-state';
 export * from './app/liver-section-volumes';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const path='.local-models/geometry-workbench/liver-completion/completion-v200/atlas.json';
const raw=fs.existsSync(path)?JSON.parse(fs.readFileSync(path)):null,skip={skip:!raw&&'Private local sources are not redistributed.'};
const walk=n=>[n,...(n.children??[]).flatMap(walk)];
test('63 actual liver representations survive ordinary source cleanup',skip,()=>{
 const a=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(raw)));
 const m=api.liverCompletionRegistry(a);assert.ok(m);assert.equal(m.partIds.length,16);
 assert.ok(m.partIds.every(id=>a.parts.some(p=>p.id===id&&!p.suppressed)));
 const all=new Map([...api.liverCompletionTermParts(a),...api.liverSurfaceTermParts(a)]);
 for(const g of m.termGroups)assert.deepEqual(all.get(g.term).map(p=>p.id),g.partIds);
 const livers=walk(api.hierarchyNavigation(a,'systems')).filter(n=>n.name==='Liver'&&n.terminology?.ta98==='A05.8.01.001');assert.equal(livers.length,1,'all new liver representations share one Systems organ');const root=livers[0];
 const nodes=walk(root).filter(n=>n.terminology?.ta98?.startsWith('A05.8.01.'));
 assert.equal(new Set(nodes.map(n=>n.terminology.ta98)).size,69);
 assert.equal(nodes.filter(n=>n.elements.length).length,63);
 for(const id of [...m.functionalPartIds,...m.microscopicPartIds])assert.ok(!root.elements.includes(id),'whole-organ selection cannot include an alternate tissue presentation or magnified fragment');
});
test('boundary controls use ordinary labels and retain the same physical host when cleared',skip,()=>{
 const references=api.liverReferenceChoices(raw);
 for(const g of raw.liverCompletionIntegration.termGroups.filter(g=>g.curves)){
  const c=raw.concepts.find(c=>c.id===g.conceptId);const labels=api.selectedAnatomyLabels(raw,c.elements,c,false);
  assert.equal(labels.length,1);assert.deepEqual(labels[0].anchor.point,g.anchor.point);
  const r=references.get(g.term);assert.ok(r);assert.deepEqual(api.toggleLiverReferenceChoice(r,c).elements,c.elements);
 }
});
test('functional presentation excludes duplicate gross liver tissue in both scene and section',skip,()=>{
 const m=raw.liverCompletionIntegration,p=raw.parts.find(p=>p.id==='TA98PART:ZA:Liver:right');
 const state={selected:[m.functionalPartIds[0]],hidden:raw.defaultHidden.filter(id=>id!==m.functionalPartIds[0]),visible:['digestive'],depthHidden:[],isolate:false,region:'all',skinOpacity:0};
 assert.equal(api.partVisible(p,state,raw.parts),false);assert.equal(api.sectionPartVisible(p,state),false);
 const restored={...state,selected:[p.id],hidden:raw.defaultHidden,isolate:true};
 assert.equal(api.partVisible(p,restored,raw.parts),true);assert.equal(api.sectionPartVisible(p,restored),true);
});
test('microscopic terms have real separate structures and explicit magnified navigation',skip,()=>{
 for(const g of raw.liverCompletionIntegration.termGroups.filter(g=>g.detailModel)){
  const c=raw.concepts.find(c=>c.id===g.conceptId);assert.match(api.liverMicroscopicDetailUrl(c),/model=liver-unit/);
  assert.ok(raw.parts.find(p=>p.id===g.partIds[0]).indexCount>=120);
 }
 const m=raw.liverCompletionIntegration;assert.equal(m.microscopicPartIds.length,6);
 assert.equal(api.liverMicroscopicDetailUrl(raw.concepts.find(c=>c.id==='ZA:Liver')),null);
});
test('failed peritoneal coats and six unfinished terms are not substituted by generic liver geometry',skip,()=>{
 const terms=api.liverCompletionTermParts(raw);
 for(const n of [28,31,32,52,53,68])assert.equal(terms.has(`A05.8.01.${String(n).padStart(3,'0')}`),false);
 assert.equal(api.liverCompletionRegistry({...raw,liverCompletionIntegration:{...raw.liverCompletionIntegration,partIds:[]}}),null);
});
