import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/gallbladder-relations'; export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const models=Object.fromEntries(['male','female'].map(sex=>[sex,JSON.parse(readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`))]));
if(models.female.gallbladderWallIntegration){const install=JSON.parse(readFileSync('reports/ta98-female/gallbladder-current-coats-install.json'));models.female=JSON.parse(gunzipSync(readFileSync(install.snapshot)));}
test('the original native gallbladder sources retain reciprocal supply, drainage and vagal contracts',()=>{
 for(const [sex,atlas] of Object.entries(models)){
  const pieces=atlas.parts.filter(api.nativeGallbladderPart);assert.equal(pieces.length,sex==='male'?4:1);
  const graph=api.gallbladderRelations(atlas.parts);assert.equal(graph,api.gallbladderRelations(atlas.parts));
  for(const p of pieces){
   const links=api.anatomicalRelations(p,atlas.parts);
   for(const [kind,back] of [['arterial','supplies'],['venous','drains'],['innervation','innervates']]){
    const expected=graph.get(p.id).filter(r=>r.kind===kind);assert.ok(expected.length);
    for(const r of expected){assert.ok(links.some(q=>q.kind===r.kind&&q.target.id===r.target.id));assert.ok(api.anatomicalRelations(r.target,atlas.parts).some(q=>q.kind===back&&q.target.id===p.id&&q.note===r.note));assert.match(r.note,/unverified/);}
   }
   const supply=graph.get(p.id).find(r=>r.kind==='arterial');
   assert.deepEqual(supply.via,sex==='male'?['Right hepatic artery','Cystic artery']:['Cystic artery']);
   assert.equal(supply.target.id,sex==='male'?'ZA:Proper hepatic artery':'LOCAL:female:r_hepatic_artery:1719');
   if(sex==='male')assert.ok(links.some(r=>r.kind==='lymphaticDrainage'&&r.target.id==='ZA:Cystic node'));
   assert.ok(graph.get(p.id).every(r=>!r.target.id.startsWith('HRA:')&&r.target.id.startsWith(sex==='male'?'ZA:':'LOCAL:female:')));
  }
 }
});
test('source veins do not imply a universal cystic vein, and a duct or wall does not inherit whole-organ supply',()=>{
 const atlas=models.female,p=atlas.parts.find(api.nativeGallbladderPart);
 const g=api.gallbladderRelations(atlas.parts).get(p.id);
 assert.ok(g.some(r=>r.kind==='venous'&&r.via.includes('Hepatic sinusoids')));
 assert.ok(g.some(r=>r.kind==='venous'&&r.via.includes('Cystic venous tributaries')));
 assert.ok(!g.some(r=>r.target.id==='HRA:VH_F_cystic_vein'));
 for(const bad of [{...p,suppressed:true},{...p,id:'HRA:VH_F_gallbladder'},{...p,sourceId:'other'},{...p,name:'Tunica mucosa'},{...p,system:'connective'},{...p,id:'LOCAL:female:bile_duct:1031',sourceId:'bile_duct',name:'Bile duct'}])assert.equal(api.nativeGallbladderPart(bad),false);
 const right=atlas.parts.find(p=>p.id==='LOCAL:female:r_hepatic_artery:1719');
 assert.equal(api.gallbladderRelations([p,{...right,suppressed:true}]).size,0);
 assert.equal(api.gallbladderRelations([p,{...right,id:'HRA:VH_F_right_hepatic_artery'}]).size,0);
 assert.equal(api.gallbladderRelations([p,{...right,sourceId:'other'}]).size,0);
});
