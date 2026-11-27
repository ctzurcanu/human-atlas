import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/ileocaecal-reference-relations'; export {anatomicalRelations} from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {ileocaecalWallInterface,ileocaecalWallRelations,anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const atlas=JSON.parse(await fs.readFile('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
const layers=['serosa','subserosa','longitudinal-muscle','circular-muscle','submucosa','muscularis-mucosae','lamina-propria','epithelium'];
test('Every installed coat connects ileum and caecum through a reciprocal annular tissue interface',()=>{
 const graph=ileocaecalWallRelations(atlas.parts);
 for(const layer of layers){
  const c=atlas.parts.find(p=>p.id===`LOCAL:female:CAECOWALL:caecum:${layer}`),i=atlas.parts.find(p=>p.id===`LOCAL:female:BOWELWALL:ileum:${layer}`);
  assert.ok(ileocaecalWallInterface(c,i));assert.ok(ileocaecalWallInterface(i,c));
  assert.ok(graph.get(c.id)?.some(r=>r.kind==='continuous'&&r.target===i));
  assert.ok(anatomicalRelations(i,atlas.parts).some(r=>r.kind==='continuous'&&r.target.id===c.id));
  for(const changed of [{...i,suppressed:true},{...i,geometryOrigin:[0,0,0]},{...i,positions:i.positions+4},{...i,chunk:i.chunk+1},{...i,ileocaecalJunctionRole:'caecum'},{...i,ileocaecalJunctionLayer:'wrong'},{...i,tissueDrawGroups:i.tissueDrawGroups.filter(g=>g.neighbor!==c.id)},{...i,sourceOffset:[0,0,0]}])assert.equal(ileocaecalWallInterface(c,changed),false);
 }
});
test('No opening plug or reference donor is substituted for the joined native wall',()=>{
 const c=atlas.parts.find(p=>p.id==='LOCAL:female:CAECOWALL:caecum:serosa');
 assert.equal(ileocaecalWallRelations([c]).size,0);
 for(const p of atlas.parts.filter(p=>p.id.startsWith('HRA:')))assert.equal(ileocaecalWallInterface(c,p),false);
 assert.equal(atlas.ileocaecalJunctionIntegration.lumenOpeningCreated,true);
 assert.equal(atlas.ileocaecalJunctionIntegration.intracaecalTailRemoved,true);
});
