import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/colonic-taenia';export {anatomicalRelations} from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {nativeColonicTaenia,colonicTaeniaFinish,colonicTaeniaRelations,anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const a=JSON.parse(await fs.readFile('.local-models/ta98-male.json'));const ids=['ZA:Free taenia','ZA:Mesocolic taenia','ZA:Omental taenia'];
test('Only three verified native identities receive an economical smooth-muscle finish',()=>{
 const matches=a.parts.filter(nativeColonicTaenia);assert.deepEqual(matches.map(p=>p.id).sort(),ids.sort());
 for(const p of matches){assert.deepEqual(colonicTaeniaFinish(p).color,[174,91,111]);assert.equal(colonicTaeniaFinish(p).map,undefined);for(const invalid of [{...p,suppressed:true},{...p,system:'connective'},{...p,sourceId:'different'},{...p,id:'HRA:band'},{...p,name:'Other taenia'}])assert.equal(colonicTaeniaFinish(invalid),undefined);}
 for(const p of a.parts.filter(p=>p.material==='ZA:digestive:Ligament'&&!ids.includes(p.id)))assert.equal(colonicTaeniaFinish(p),undefined);
});
test('Same-source associations are reciprocal and retain missing regional portions',()=>{
 const graph=colonicTaeniaRelations(a.parts);for(const id of ids){const band=a.parts.find(p=>p.id===id),links=graph.get(id);assert.equal(links.filter(r=>r.via?.[0]==='Regional portion of the longitudinal muscle layer').length,4);for(const link of links){assert.match(link.note,/unverified/);assert.ok(graph.get(link.target.id).some(r=>r.target.id===id&&r.note===link.note&&r.kind==='adjacent'));assert.ok(anatomicalRelations(band,a.parts).some(r=>r.target.id===link.target.id&&r.note===link.note));}}
 assert.ok(graph.get('ZA:Mesocolic taenia').some(r=>r.target.id==='ZA:Mesocolon'));
 assert.ok(graph.get('ZA:Omental taenia').some(r=>r.target.id==='ZA:Greater omentum'));
 assert.ok(!graph.get('ZA:Free taenia').some(r=>/omentum|mesocolon/i.test(r.target.name)));
});
test('Absent, renamed, suppressed or independently registered sources cannot restore band associations',()=>{
 const band=a.parts.find(p=>p.id==='ZA:Free taenia'),colon=a.parts.find(p=>p.id==='ZA:Ascending colon');
 for(const other of [{...colon,suppressed:true},{...colon,id:'LOCAL:female:ascending_colon:76'},{...colon,sourceId:'other'},{...colon,name:'Other colon'}])assert.equal(colonicTaeniaRelations([band,other]).size,0);
 assert.equal(colonicTaeniaRelations([{...band,id:'LOCAL:female:taenia'}]).size,0);
});
test('Native band depth bias is bounded and retains normal occlusion',async()=>{
 for(const p of a.parts.filter(nativeColonicTaenia)){
  const finish=colonicTaeniaFinish(p);assert.equal(finish.polygonOffset,true);assert.equal(finish.polygonOffsetFactor,-1);assert.equal(finish.polygonOffsetUnits,-1);assert.notEqual(finish.depthTest,false);
 }
 const scene=await fs.readFile('app/scene.tsx','utf8');assert.match(scene,/polygonOffset:source\?\.polygonOffset\?\?false/);assert.match(scene,/polygonOffsetFactor:source\?\.polygonOffsetFactor\?\?0/);
});
