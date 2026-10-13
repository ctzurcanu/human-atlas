import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
const read=path=>readFile(path,'utf8').then(JSON.parse);
test('Superior duodenal source is a regional identity, while whole duodenum retains its source',async t=>{
 const metadata=await read('app/data/ta98-metadata.json');
 assert.equal(metadata.byConcept['HRA:VH_F_duodenum_superior'].ta98,'A05.6.02.002');
 assert.deepEqual(metadata.reviewedCatalogueIdentities['HRA:VH_F_duodenum_superior'],{term:'A05.6.02.002',kind:'exact'});
 let atlas,coverage;
 try{atlas=await read('.local-models/ta98-female-runtime.json');coverage=await read('reports/ta98/coverage.json');}
 catch(e){if(e.code==='ENOENT'){t.skip('Local source review unavailable');return;}throw e;}
 const source=atlas.concepts.find(c=>c.id==='HRA:VH_F_duodenum_superior');
 assert.equal(source.ta98Term,'A05.6.02.002');
 assert.equal(atlas.concepts.find(c=>c.id==='TA98FEMALE:A05.6.02.002').name,'Superior part of duodenum');
 assert.ok(source.elements.includes('HRA:VH_F_duodenal_ampulla'),'ampulla remains a superior-part component');
 assert.equal(atlas.concepts.find(c=>c.id==='LOCAL:female:intestine_duodenum').ta98Term,'A05.6.02.001');
 const superior=coverage.terms.find(t=>t.id==='A05.6.02.002'),whole=coverage.terms.find(t=>t.id==='A05.6.02.001');
 for(const model of ['female','local-female-ta98']){
  assert.ok(superior.models[model].exact.some(c=>c.id===source.id));
  assert.ok(!whole.models[model].exact.some(c=>c.id===source.id),'a regional donor cannot establish the whole term');
 }
});
