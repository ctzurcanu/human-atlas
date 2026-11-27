import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:`export * from './app/liver-controls';export * from './app/liver-completion-controls';export * from './app/liver-surface-features';export * from './app/hepatic-connection-controls';export * from './app/liver-section-volumes';export * from './app/hierarchy-navigation';export * from './app/preferred-source-representations';export * from './app/reviewed-catalogue-identities';`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const walk=n=>[n,...(n.children??[]).flatMap(walk)];
for(const sex of['male','female']){
 const path=`.local-models/geometry-workbench/liver-completion/final-clearance-v270/${sex}/atlas.json`;
 test(`${sex}: all 69 liver terms have live Systems and TA98 controls`,{skip:!fs.existsSync(path)&&'Private local geometry'},()=>{
  const raw=JSON.parse(fs.readFileSync(path)),a=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(raw)));
  assert.ok(api.liverControlRegistry(a),'gross tissue registry');assert.ok(api.liverCompletionRegistry(a),'completion registry');assert.ok(api.liverSectionRegistry(a),'section registry');assert.ok(api.hepaticConnectionRegistry(a),'conduit registry');assert.ok(api.liverFeatureRegistry(a),'surface registry');
  const systems=walk(api.hierarchyNavigation(a,'systems')),liver=systems.find(n=>n.name==='Liver'&&n.terminology?.ta98==='A05.8.01.001');assert.ok(liver);
  const terms=walk(liver).filter(n=>n.terminology?.ta98?.startsWith('A05.8.01.'));assert.equal(new Set(terms.map(n=>n.terminology.ta98)).size,69);assert.equal(terms.filter(n=>n.elements.length).length,69);
  const authority=JSON.parse(fs.readFileSync('public/assets/ta98.json')),ta=walk(api.hierarchyNavigation(a,'guest:ta98',authority));
  for(const g of [...a.liverIntegration.termGroups,...a.hepaticConnectionIntegration.termGroups].filter(g=>g.term.startsWith('A05.8.01.'))){const node=ta.find(n=>n.id==='hierarchy:guest:ta98:'+g.term);assert.ok(node?.elements.length,g.term);}
  const whole=a.concepts.find(c=>c.id===a.liverIntegration.wholeConceptId);assert.equal(api.liverTissuePeelSteps(a,whole).length,4);
 });
}
