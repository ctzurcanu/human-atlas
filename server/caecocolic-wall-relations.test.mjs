import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/caecocolic-wall-relations';export {anatomicalRelations} from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {caecocolicWallRelations,anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
for(const sex of ['female','male']){
 const atlas=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),parts=atlas.parts;
 test(sex+': all five regional interfaces connect reciprocally through eight stored coats',()=>{
  const graph=caecocolicWallRelations(parts);const colonicLinks=[...graph.entries()].filter(([id])=>!id.includes(':rectum:')).flatMap(([,links])=>links.filter(r=>r.target.caecocolicRegion!=='rectum'));assert.equal(colonicLinks.length,atlas.caecocolicWallIntegration.taeniaPartition?130:80);
  for(const [id,links] of graph){const source=parts.find(p=>p.id===id);for(const r of links){assert.equal(r.kind,'continuous');assert.ok(graph.get(r.target.id).some(back=>back.target===source&&back.note===r.note));assert.equal(source.caecocolicLayer??id.split(':').at(-1),r.target.caecocolicLayer??r.target.id.split(':').at(-1));assert.ok(anatomicalRelations(source,parts).some(x=>x.kind===r.kind&&x.target===r.target));}}
  const epithelial=parts.find(p=>p.id===`LOCAL:${sex}:CAECOWALL:caecum:epithelium`);assert.deepEqual(graph.get(epithelial.id).map(r=>r.target.caecocolicRegion).sort(),['appendix','ascending']);
 });
 test(sex+': moved, foreign, suppressed and nonreciprocal banks cannot imply tissue continuity',()=>{
  const a=parts.find(p=>p.id===`LOCAL:${sex}:CAECOWALL:ascending:serosa`),b=parts.find(p=>p.id===`LOCAL:${sex}:CAECOWALL:transverse:serosa`);
  assert.equal(caecocolicWallRelations([a,b]).size,2);
  for(const change of [{suppressed:true},{positions:b.positions+4},{normals:b.normals+4},{chunk:b.chunk+1},{geometryOrigin:[0,0,0]},{sourceOffset:[0,0,0]},{caecocolicWallAssembly:'foreign'},{tissueDrawGroups:[]},{system:'arterial'}])assert.equal(caecocolicWallRelations([a,{...b,...change}]).size,0);
  assert.equal(caecocolicWallRelations([a]).size,0);
 });
}
test('female staged rectal coats expose only reciprocal stored sigmoid junctions',()=>{
 const atlas=JSON.parse(fs.readFileSync('.local-models/ta98-review/rectal-wall-main-atlas/female/atlas.json')),parts=atlas.parts;
 const source=parts.find(p=>p.id==='LOCAL:female:large_intestine_rectum:80');
 const rectal=parts.filter(p=>p.caecocolicRegion==='rectum');assert.equal(rectal.length,8);
 const graph=caecocolicWallRelations(parts);
 for(const coat of rectal){
  assert.equal(coat.sourceId,source.sourceId);
  const links=graph.get(coat.id);assert.ok(links?.length);
  assert.ok(links.every(r=>r.kind==='continuous'&&r.target.caecocolicRegion==='sigmoid'&&r.target.caecocolicLayer===coat.caecocolicLayer));
  for(const link of links){
   assert.ok(graph.get(link.target.id)?.some(back=>back.target===coat));
   assert.ok(anatomicalRelations(coat,parts).some(r=>r.kind==='continuous'&&r.target===link.target));
  }
  for(const change of [{positions:coat.positions+4},{sourceOffset:[0,0,0]},{suppressed:true},{tissueDrawGroups:[]}]){
   const altered=parts.map(p=>p.id===coat.id?{...p,...change}:p);assert.equal(caecocolicWallRelations(altered).has(coat.id),false);
  }
 }
});
