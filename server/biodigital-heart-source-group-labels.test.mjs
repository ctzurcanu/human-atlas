import assert from 'node:assert/strict';
import {readFileSync,existsSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';

const profile=process.env.ATLAS_HEART_LABEL_PROFILE??'.local-models/geometry-workbench/biodigital-male-circulation-v1/heart-additional-opening-labels-v1/candidate-v1/atlas.json';
const local=existsSync(profile);
let api,atlas;
if(local){
 const out=await build({stdin:{contents:"export * from './app/cardiac-label-sites';export * from './app/anatomy-label-plan';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
 api=await import('data:text/javascript;base64,'+Buffer.from(out.outputFiles[0].contents).toString('base64'));
 atlas=JSON.parse(readFileSync(profile,'utf8'));
}
const sourceId='BDM:male:human_20_male_cardiovascular_system-heart_ID';
test('donor Heart expansion exposes every installed cardiac term and its original vessel leaves',{skip:!local},()=>{
 const m=atlas.sourceCardiacIntegration,whole=m.termGroups.find(g=>g.term==='A12.1.00.001'),choice=atlas.concepts.find(c=>c.id===sourceId);
 const entries=api.cardiacComponentChoices(atlas,choice),terms=new Set(entries.map(c=>c.ta98Term).filter(Boolean));
 for(const g of m.termGroups.filter(g=>g.term!==whole.term))assert.ok(terms.has(g.term),g.term);
 for(const s of m.referenceSites)assert.equal(entries.filter(c=>c.ta98Term===s.term).length,1,s.term);
 const physical=new Set(whole.partIds),extras=choice.elements.filter(id=>!physical.has(id));
 for(const id of extras)assert.equal(entries.filter(c=>!c.ta98Term?.startsWith('A12.1.')&&c.elements.includes(id)).length,1,id);
 assert.equal(entries.length,m.termGroups.length-1+extras.length);
 assert.equal(new Set(entries.map(c=>c.id)).size,entries.length);
 const labels=api.selectedAnatomyLabels(atlas,choice.elements,choice,true);
 assert.equal(labels.length,entries.length);
 for(const s of m.referenceSites)assert.deepEqual(labels.find(l=>l.id===s.conceptId)?.anchor,s.anchor);
});
test('source group expansion never invents a donor parent for a partial or unrelated selection',{skip:!local},()=>{
 const choice=atlas.concepts.find(c=>c.id===sourceId);
 assert.equal(api.cardiacComponentChoices(atlas,{...choice,elements:choice.elements.slice(1)}),undefined);
 assert.equal(api.cardiacComponentChoices(atlas,{...choice,id:'some-unrelated-source'}),undefined);
 assert.equal(api.cardiacComponentChoices(atlas,{...choice,elements:[...choice.elements,choice.elements[0]]}),undefined);
});
test('canonical heart references keep their individual carriers and reject stale cardiac poses',{skip:!local},()=>{
 const m=atlas.sourceCardiacIntegration,whole=m.termGroups.find(g=>g.term==='A12.1.00.001'),choice=atlas.concepts.find(c=>c.id===whole.conceptId);
 assert.equal(api.cardiacComponentChoices(atlas,choice).length,m.termGroups.length-1);
 for(const s of m.referenceSites)assert.deepEqual(api.cardiacComponentChoices(atlas,atlas.concepts.find(c=>c.id===s.conceptId)),[]);
 const broken=structuredClone(atlas),p=broken.parts.find(p=>p.id===whole.partIds[0]);p.geometryOrigin[0]+=.001;
 assert.equal(api.cardiacComponentChoices(broken,broken.concepts.find(c=>c.id===sourceId)),undefined);
});
