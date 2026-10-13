import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/gallbladder-wall';export * from './app/anatomical-components';export * from './app/anatomical-relations';export * from './app/guest-hierarchy';export * from './app/local-ta98-layer-catalogue';export * from './app/anatomy-label-plan';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const read=p=>JSON.parse(fs.readFileSync(p));
const staged=read('.local-models/ta98-review/female-gallbladder-regions/atlas.json');
const installed=read('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json');
for(const [phase,a] of [['staged',staged],...(installed.femaleGallbladderRegions?[['installed',installed]]:[])]){
 const parts=new Map(a.parts.map(p=>[p.id,p]));
 test(`${phase}: four exact region controls select their four disjoint wall coats and peel stages`,()=>{
  const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};api.resolveGuestHierarchy(a,api.parseGuestHierarchy(read('public/assets/ta98.json'))).forEach(walk);
  for(const suffix of ['002','003','004','005']){
   const term='A05.8.02.'+suffix,choice=a.concepts.find(c=>c.id==='GALLREGION:female:'+term),members=choice.elements.map(id=>parts.get(id));
   assert.equal(members.length,4);assert.ok(members.every(p=>p.pairedTissue&&p.sectionSolid&&!p.suppressed&&p.gallbladderRegion===term));
   assert.deepEqual(new Set(nodes.get(term).directParts.map(p=>p.id)),new Set(choice.elements));
   const steps=api.gallbladderWallPeelSteps(a,choice);assert.equal(steps.length,4);assert.ok(steps.every(s=>s.ids.length===1));
   const items=api.selectedAnatomyLabels(a,choice.elements,choice,true),targets=api.layerLabelTargets(items,a.parts,a.concepts);assert.equal(targets.size,4);
   assert.equal(new Set([...targets.values()].map(t=>t.coordinate)).size,4);
  }
 });
 test(`${phase}: every neck coat retains the cystic interface and whole-organ supply returns all sixteen cells`,()=>{
  const choice=a.concepts.find(c=>c.id==='LOCAL:female:gallbladder'),members=choice.elements.map(id=>parts.get(id));assert.equal(members.length,16);
  assert.ok(api.gallbladderWholeSource(members,a.parts));assert.equal(api.gallbladderWallPeelSteps(a,choice).length,4);
  const duct=parts.get(a.cysticFoldIntegration.ductPartId),graph=api.gallbladderWallRelations(a.parts);
  for(const p of members){const connected=p.gallbladderRegion==='A05.8.02.005';assert.equal(api.gallbladderWallInterface(p,duct),connected);if(connected)assert.ok(graph.get(p.id).some(r=>r.kind==='continuous'&&r.target.id===duct.id));}
  const links=api.anatomicalConceptFunctionalRelations(choice,members,a.parts);assert.ok(links.some(r=>r.kind==='arterial'));
  for(const row of links.filter(r=>['arterial','venous','innervation'].includes(r.kind))){const back=api.anatomicalRelations(row.target,a.parts).find(r=>r.targetConcept?.id===choice.id);assert.ok(back);assert.deepEqual(new Set(back.targetConcept.elements),new Set(choice.elements));}
 });
}
test('independent packed proof retains the macro lumen, all original controls and reciprocal normals',()=>{
 const audit=read('reports/ta98-female/gallbladder-regions-audit.json'),proof=read('reports/ta98-female/gallbladder-regions-packed-validation.json'),plan=read('.local-models/ta98-review/female-gallbladder-regions/construction.json');
 assert.equal(audit.catalogueSha256,plan.catalogueSha256);assert.equal(proof.catalogueSha256,plan.catalogueSha256);assert.equal(audit.bankSha256,plan.bankSha256);
 assert.equal(audit.cells.length,16);assert.ok(audit.cells.every(c=>c.components===1&&c.closed&&c.selfIntersections===0));assert.ok(audit.originalCoordinatesAndNormalsExact&&audit.continuousMouthCoordinatesExact&&audit.unrelatedPartsAndDefaultHiddenExact);
 assert.equal(audit.sharedCoatSignFields,4);assert.equal(audit.reciprocalContacts.length,30);assert.ok(proof.continuousMacroLumen);assert.ok(proof.macroBoundaries.every(s=>s.components===1&&s.openTerminalLoops===1&&s.euler===1&&s.selfIntersections===0));
 assert.equal(audit.continuousTissuePrismJacobianMaxima.length,4);assert.ok(audit.continuousTissuePrismJacobianMaxima.every(value=>value<0));
 assert.equal(plan.newUVBytes,0);assert.equal(plan.newImageBytes,0);assert.ok(plan.incrementalGzipBytes<150000);
});
test('regional staging cannot replace public or male anatomy',()=>{
 const q=new URLSearchParams('biliary=regions');assert.equal(api.localTa98LayerCatalogue('female',q,true),undefined);assert.equal(api.localTa98LayerCatalogue('local-female-ta98',q,false),undefined);assert.equal(api.localTa98LayerCatalogue('local-ta98',q,true),'/local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json');
});
