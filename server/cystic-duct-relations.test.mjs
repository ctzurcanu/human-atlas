import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/digestive-duct-relations'; export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const atlas=JSON.parse(readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
const cystic=atlas.parts.find(p=>p.id==='FJ3080'),trunks=atlas.parts.filter(p=>p.sourceId==='FJ3079'&&!p.suppressed);
test('current male cystic duct exposes one reciprocal hepatic/bile junction without tissue containment',()=>{
 const relations=api.anatomicalRelations(cystic,atlas.parts);
 assert.equal(trunks.length,2);
 for(const trunk of trunks){
  const kind=atlas.biliaryJunctionIntegration?.continuousConfluence?'continuous':'adjacent';
  assert.ok(relations.some(r=>r.kind===kind&&r.target.id===trunk.id&&/bidirectional/.test(r.note)));
  assert.ok(api.anatomicalRelations(trunk,atlas.parts).some(r=>r.kind===kind&&r.target.id===cystic.id));
  assert.ok(!relations.some(r=>r.target.id===trunk.id&&['contains','partOf'].includes(r.kind)));
 }
 assert.match(api.digestiveDuctDescription(cystic).text,atlas.cysticFoldIntegration?/Three raised crescentic/:/spiral fold.*under review/);
});
test('suppressed, displaced, foreign or misidentified duct sources cannot inherit these routes',()=>{
 for(const change of [{suppressed:true},{sourceOffset:[0,0,0]},{sourceId:'foreign'},{conceptId:'wrong'},{system:'reproductive'},{id:'HRA:cystic'}])
  assert.equal(api.digestiveDuctRelations([{...cystic,...change},...trunks]).has(cystic.id),false);
 for(const change of [{suppressed:true},{sourceOffset:[0,0,0]},{sourceId:'foreign'},{conceptId:'wrong'},{system:'reproductive'}])
  assert.equal(api.digestiveDuctRelations([cystic,...trunks.map(p=>({...p,...change}))]).size,0);
});
test('female compound source never gains a guessed separate cystic duct from hidden donor geometry',()=>{
 const female=JSON.parse(readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
 const graph=api.digestiveDuctRelations(female.parts);
 for(const links of graph.values())assert.ok(links.every(r=>r.target.id!=='FJ3080'&&!/^HRA:.*cystic/.test(r.target.id)));
});
