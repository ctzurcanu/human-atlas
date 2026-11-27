import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
import {Quaternion,Vector3} from 'three';

const built=await build({stdin:{contents:"export * from './app/cardiac-label-sites';export * from './app/heart-wall';export * from './app/heart-system';export * from './app/guest-hierarchy';export * from './app/anatomy-label-plan';export * from './app/reviewed-catalogue-identities';export * from './app/preferred-source-representations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const authority=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
const directory=process.env.ATLAS_CARDIAC_LABEL_BUNDLE_STAGE??'.local-models/ta98-review/bowel-wall-local-atlas';
const data={};
for(const sex of ['male','female']){
 const raw=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),bundle=JSON.parse(fs.readFileSync(`${directory}/${sex}/cardiac-label-sites.json`));
 const atlas=api.withCardiacLabelBundle(raw,bundle);data[sex]={raw,bundle,atlas};
 test(`${sex}: the ordinary page loading pipeline retains all 113 controls`,()=>{
  const prepared=api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(raw));
  const loaded=api.withHiddenReviewAlternatives(api.withCardiacLabelBundle(prepared,bundle));
  assert.equal(api.cardiacLabeledTermParts(loaded).size,113);
  assert.equal(api.withHeartSystem(loaded,[])[0].itemCount,113);
 });
 test(`${sex}: 113 selectable terms retain a separate count of actual tissue controls`,()=>{
  const m=api.cardiacLabelRegistry(atlas);assert.ok(m);
  assert.equal(api.cardiacLabeledTermParts(atlas).size,113);
  assert.equal(api.cardiacTermParts(atlas).size,42);
  assert.equal(m.termGroups.length,71);assert.ok(m.termGroups.every(g=>g.referenceOnly&&g.tissueGeometryEstablished===false));
  assert.strictEqual(atlas.parts,raw.parts);assert.strictEqual(atlas.chunks,raw.chunks);assert.strictEqual(atlas.materials,raw.materials);
  assert.equal(atlas.concepts.length,raw.concepts.length+71);
  for(const g of m.termGroups){const concept=atlas.concepts.find(c=>c.id===g.conceptId);assert.match(concept.coverageLimitation,/Provisional/);assert.equal(api.cardiacSiteCoverage(atlas,concept),concept.coverageLimitation);assert.equal(api.cardiacSiteConceptCoverage(atlas.concepts,concept),concept.coverageLimitation);assert.ok(g.anchor.point.every(Number.isFinite));}
 });
 test(`${sex}: all labels resolve in both the TA98 and Heart Systems trees`,()=>{
  const queue=[...api.resolveGuestHierarchy(atlas,authority)],byTerm=new Map();
  while(queue.length){const n=queue.pop();if(byTerm.has(n.id))continue;byTerm.set(n.id,n);queue.push(...n.children);}
  assert.equal(byTerm.get('A12.1.00.001').itemCount,101);
  for(const [term,parts]of api.cardiacLabeledTermParts(atlas))assert.ok(parts.every(p=>byTerm.get(term)?.directParts.some(q=>q.id===p.id)),term);
  const hearts=api.withHeartSystem(atlas,[]);assert.equal(hearts.length,1);const nodes=[hearts[0]],labels=new Set();
  assert.equal(hearts[0].itemCount,113);
  while(nodes.length){const n=nodes.pop();const term=n.kind==='entry'?n.entry.ta98Term:n.terminology?.ta98;if(term)labels.add(term);if(n.kind==='group')nodes.push(...n.nodes);}
  for(const term of api.cardiacLabeledTermParts(atlas).keys())assert.ok(labels.has(term),term);
  const foramen=atlas.concepts.find(c=>c.ta98Term==='A12.1.01.007'&&c.cardiacReferenceSite);
  assert.match(foramen.coverageLimitation,/adult.*closed|closed adult/i);
 });
 test(`${sex}: each reference uses the normal label and selects only its own item`,()=>{
  for(const g of api.cardiacLabelRegistry(atlas).termGroups){
   const choice=atlas.concepts.find(c=>c.id===g.conceptId),items=api.selectedAnatomyLabels(atlas,choice.elements,choice,false);
   assert.equal(items.length,1,g.term);assert.equal(items[0].id,choice.id);assert.deepEqual(items[0].anchor,g.anchor);
   const selected=api.anatomyLabelChoice(items[0]);assert.equal(selected.id,choice.id);assert.deepEqual(selected.elements,choice.elements);
  }
 });
 test(`${sex}: opening the Heart label expands 112 distinct anatomical component labels`,()=>{
  const whole=atlas.concepts.find(c=>c.id===atlas.finiteCardiacIntegration.wholeConceptId);
  const components=api.cardiacComponentChoices(atlas,whole);
  assert.equal(components.length,112);
  assert.equal(new Set(components.map(c=>c.ta98Term)).size,112);
  const labels=api.selectedAnatomyLabels(atlas,whole.elements,whole,true);
  assert.equal(labels.length,112);
  assert.equal(new Set(labels.map(item=>item.id)).size,112);
  for(const item of labels){const selected=api.anatomyLabelChoice(item);assert.equal(selected.id,item.id);assert.deepEqual(selected.elements,item.members);}
 });
 test(`${sex}: stale carrier poses and nonfinite label anchors fail closed`,()=>{
  const bad=structuredClone(raw);bad.parts.find(p=>p.id===bundle.integration.termGroups[0].anchor.partId).geometryOrigin=[2,3,4];
  assert.strictEqual(api.withCardiacLabelBundle(bad,bundle),bad);
  const invalid=structuredClone(bundle);invalid.integration.termGroups[0].anchor.point[0]=NaN;
  assert.strictEqual(api.withCardiacLabelBundle(raw,invalid),raw);
  const inflated=structuredClone(bundle);inflated.integration.geometryBackedTermCount++;
  inflated.integration.labeledTermCount++;
  assert.strictEqual(api.withCardiacLabelBundle(raw,inflated),raw);
 });
}

test('female points use the same positive similarity as the shared cardiac bank',()=>{
 const local=(atlas,point)=>{const p=atlas.parts.find(p=>p.id===`TA98HEART:${atlas.finiteCardiacIntegration.sex}:epicardium`);return new Vector3(...point).sub(new Vector3(...p.geometryOrigin)).divideScalar(p.geometryInstance.scale).applyQuaternion(new Quaternion(...p.geometryInstance.rotation).invert());};
 const male=api.cardiacLabelRegistry(data.male.atlas),female=api.cardiacLabelRegistry(data.female.atlas);
 for(const g of male.termGroups){const other=female.termGroups.find(f=>f.term===g.term);assert.ok(local(data.male.atlas,g.anchor.point).distanceTo(local(data.female.atlas,other.anchor.point))<1e-12,g.term);}
});

test('the unfitted female pulmonary source cannot reappear through the primary trunk control',()=>{
 const a=data.female.raw,old='LOCAL:female:pulmonary_trunk:1508',current='LOCAL:female:heart-template:za-pulmonary-trunk';
 const source=a.parts.find(p=>p.id===old);assert.ok(source.suppressed&&source.defaultHidden);
 assert.ok(a.defaultHidden.includes(old));assert.ok(a.finiteCardiacIntegration.rootPartIds.includes(current));
 assert.deepEqual(a.concepts.find(c=>c.id==='LOCAL:female:pulmonary_trunk').elements,[current]);
 assert.notEqual(source.conceptId,'LOCAL:female:pulmonary_trunk');
});
