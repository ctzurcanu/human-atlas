import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import vm from 'node:vm';
const source=fs.readFileSync('scripts/import-local-zygote.mjs','utf8'),start=source.indexOf('function systemFor('),end=source.indexOf('\nfunction decode(',start),systemFor=vm.runInNewContext('('+source.slice(start,end)+')',{systems:{2:'integumentary',3:'muscular',7:'arterial'}});
test('All seven explicit native male large-intestinal leaves override inherited genital routes',()=>{
 const a=JSON.parse(fs.readFileSync('.local-models/male.json'));const leaves=a.parts.filter(p=>/^large_intestine_(cecum|ascending_colon|transverse_colon|descending_colon|sigmoid_colon|rectum|anus)$/.test(p.sourceId));assert.equal(leaves.length,7);for(const p of leaves)assert.equal(systemFor(p.sourceId,{layer:6,route:['Reproductive','Male genital organs']}),'digestive',p.sourceId);
});
test('Explicit organ correction preserves vessel/muscle classes and real reproductive organs',()=>{
 assert.equal(systemFor('large_intestine_rectum_artery',{layer:7,route:['Reproductive']}),'arterial');assert.equal(systemFor('puborectalis',{layer:3,route:['Reproductive']}),'muscular');for(const n of ['uterus','prostate','testis'])assert.equal(systemFor(n,{layer:6,route:['Reproductive']}),'reproductive');
});
