import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const built=await build({stdin:{contents:`
 export * from './app/hierarchy-navigation';
 export * from './app/anatomy-hierarchy';
 export * from './app/pulmonary-system';
 export * from './app/preferred-source-representations';
 export * from './app/reviewed-catalogue-identities';
 export * from './app/viewer-state';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const inventory=api.pulmonarySystemInventory;
const terms=new Set(inventory.terms.map(row=>row.id));
const walk=node=>[node,...(node.children??[]).flatMap(walk)];
const term=node=>node.ta98Term??node.terminology?.ta98;

test('the Systems lung inventory contains the 70 authoritative pulmonary entries, excluding mediastinum',()=>{
 // The public authority graph is the same deterministic input used by the generator.
 const nodes=new Map(JSON.parse(fs.readFileSync('public/assets/ta98.json')).nodes.map(node=>[node.id,node]));
 assert.equal(terms.size,70);assert.equal(inventory.terms.length,70);
 for(const row of inventory.terms)assert.equal(row.name,nodes.get(row.id)?.name);
 assert.ok(!inventory.terms.some(row=>/mediastinum$/.test(row.name)));
 assert.equal(inventory.terms.find(row=>row.id==='A06.5.02.015').parent,'A06.5.01.003');
 assert.equal(inventory.terms.find(row=>row.id==='A06.5.02.020').parent,'A06.5.01.003');
});

for(const sex of ['female','male']){
 const file=`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`;
 const atlas=fs.existsSync(file)?api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(JSON.parse(fs.readFileSync(file))))):null;
 const options={skip:!atlas&&'Local review geometry is not redistributed.'};

 test(`${sex}: all 70 items appear once under Systems > Respiratory system > Lungs`,options,()=>{
  const root=api.hierarchyNavigation(atlas,'systems');
  const respiratory=root.children.find(node=>node.name==='Respiratory system');assert.ok(respiratory);
  const lungs=respiratory.children.filter(node=>node.name==='Lungs');assert.equal(lungs.length,1);
  const pulmonary=walk(lungs[0]).filter(node=>terms.has(term(node)));
  assert.equal(pulmonary.length,70);assert.deepEqual(new Set(pulmonary.map(term)),terms);
  const left=pulmonary.find(node=>term(node)==='A06.5.01.003');
  const leftTerms=new Set(walk(left).map(term));
  assert.ok(leftTerms.has('A06.5.02.015'));assert.ok(leftTerms.has('A06.5.02.020'));
  const unresolved=['A06.5.02.026','A06.5.02.027','A07.1.02.005','A07.1.01.001'];
  for(const code of unresolved)assert.deepEqual(pulmonary.find(node=>term(node)===code).elements,[],code+' must not borrow a whole lung');
 });

 test(`${sex}: existing 19 territories and 27 bronchi remain selectable without scattered duplicate leaves`,options,()=>{
  const root=api.hierarchyNavigation(atlas,'systems');
  const lungs=walk(root).find(node=>node.name==='Lungs'&&term(node)==='A06.5.01.001');
  const expected=new Set([...atlas.pulmonarySegmentIntegration.partIds,...atlas.pulmonaryBranchIntegration.partIds]);
  for(const id of expected)assert.ok(lungs.elements.includes(id),id);
  assert.equal(new Set(lungs.elements).size,lungs.elements.length);
  const bronchial=lungs.children.find(node=>node.name==='Bronchial tree');
  assert.deepEqual(new Set(bronchial.elements),new Set(atlas.pulmonaryBranchIntegration.partIds));
  assert.equal(bronchial.children.length,27);
  const outside=walk(root).filter(node=>!walk(lungs).includes(node)&&!node.children);
  assert.ok(!outside.some(node=>node.elements.some(id=>expected.has(id))),'no pulmonary mesh is left scattered elsewhere in Systems');
  const available=new Map(atlas.parts.filter(part=>!part.suppressed).map(part=>[part.id,part]));
  for(const node of walk(lungs))for(const id of node.elements)assert.ok(available.has(id));
  for(const code of ['A06.5.01.017','A06.5.01.018','A06.5.01.019','A06.5.01.020'])assert.ok(walk(lungs).find(node=>term(node)===code).elements.length,code+' collective selection');
 });

 test(`${sex}: Systems lung grouping preserves catalogue geometry and the separate TA98 hierarchy`,options,()=>{
  const before=JSON.stringify(atlas),authority=JSON.parse(fs.readFileSync('public/assets/ta98.json'));
  const original=api.hierarchyNavigation(atlas,'guest:ta98',authority);
  api.hierarchyNavigation(atlas,'systems');
  assert.equal(JSON.stringify(atlas),before);
  assert.deepEqual(api.hierarchyNavigation(atlas,'guest:ta98',authority),original);
 });
}
