import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const compiled=await build({stdin:{contents:`
 export * from './app/diaphragm-controls';
 export * from './app/hierarchy-navigation';
 export * from './app/preferred-source-representations';
 export * from './app/reviewed-catalogue-identities';
 export * from './app/viewer-state';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const hierarchy=JSON.parse(fs.readFileSync('public/assets/ta98.json'));
function flatten(root){const found=[],seen=new Set();const walk=n=>{if(seen.has(n))return;seen.add(n);found.push(n);n.children?.forEach(walk);};walk(root);return found;}

for(const sex of ['male','female']){
const file=`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`;
const atlas=fs.existsSync(file)?api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(JSON.parse(fs.readFileSync(file))))):null;
const registry=atlas&&api.diaphragmControlRegistry(atlas);
const whole=registry?.marker.termGroups.find(g=>g.term==='A04.4.02.001');
const ids=new Set(whole?.partIds);
const options={skip:atlas?false:'Local review catalogue is not distributed with the application.'};

test(`${sex}: Systems and TA98 select the same complete diaphragm as the direct atlas link`,options,()=>{
 assert.ok(registry,'the installed complete diaphragm is required by this regression');
 const direct=api.resolveSelection(atlas,[whole.conceptId]);
 assert.deepEqual(new Set(direct),ids);
 for(const [mode,guest]of [['systems',undefined],['guest:ta98',hierarchy]]){
  const choices=flatten(api.hierarchyNavigation(atlas,mode,guest));
  const controls=choices.filter(c=>c.children&&(c.ta98Term??c.terminology?.ta98)===whole.term);
  assert.equal(controls.length,1,mode+' has one whole diaphragm control');
  assert.deepEqual(new Set(controls[0].elements),ids,mode+' includes all muscle, tendon, ligament and fascia cells');
  const children=flatten(controls[0]).filter(c=>!c.children).flatMap(c=>c.elements);
  for(const id of ids)assert.ok(children.includes(id),mode+' retains a control for '+id);
 }
});

test(`${sex}: a partial whole-organ tree control is corrected without expanding a separately selected tissue`,options,()=>{
 assert.ok(registry);
 const group={id:'hierarchy:systems:diaphragm',name:'Diaphragm',elements:[whole.partIds[0]],children:[],terminology:{ta98:whole.term}};
 assert.deepEqual(new Set(api.diaphragmAssemblyChoice(atlas,group).elements),ids);
 const crural=atlas.parts.find(p=>p.id==='TA98PART:ZA:Diaphragm:Tendon:crural-attachment');
 const tissue={id:crural.conceptId,name:crural.name,elements:[crural.id],ta98Term:whole.term};
 assert.equal(api.diaphragmAssemblyChoice(atlas,tissue),tissue);
 const other={...group,terminology:{ta98:'A04.5.04.001'}};
 assert.equal(api.diaphragmAssemblyChoice(atlas,other),other);
});

test(`${sex}: the Systems browser assigns each diaphragm cell once while preserving tissue identity`,options,()=>{
 assert.ok(registry);
 const root=api.hierarchyNavigation(atlas,'systems');
 const leaves=flatten(root).filter(c=>!c.children);
 for(const id of ids)assert.equal(leaves.flatMap(c=>c.elements).filter(member=>member===id).length,1,id);
 assert.equal(atlas.parts.find(p=>p.id==='ZA:Diaphragmatic fascia').system,'fascia');
});
}
