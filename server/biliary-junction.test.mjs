import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/biliary-junction';export * from './app/digestive-duct-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const a=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
test('installed confluence connects cystic and both trunk directions through shared finite tissue',()=>{
 assert.ok(a.biliaryJunctionIntegration.continuousConfluence);
 const members=a.biliaryJunctionIntegration.partIds.map(id=>a.parts.find(p=>p.id===id));
 for(let i=0;i<members.length;i++)for(let j=i+1;j<members.length;j++){
  const p=members[i],q=members[j];assert.ok(api.biliaryJunctionInterface(p,q));
  assert.ok(api.digestiveDuctRelations(a.parts).get(p.id).some(r=>r.kind==='continuous'&&r.target.id===q.id));
  for(const change of [{suppressed:true},{sourceOffset:[0,0,0]},{chunk:p.chunk+1},{positions:p.positions+4},{sourceId:'foreign'},{conceptId:'foreign'},{pairedTissue:false},{sectionSolid:false},{tissueDrawGroups:[]}])assert.equal(api.biliaryJunctionInterface({...p,...change},q),false);
 }
 const proof=JSON.parse(fs.readFileSync('reports/ta98-male/biliary-junction-packed-validation.json'));
 assert.ok(proof.continuousMacroLumen);
 for(const row of proof.macroBoundaries)assert.deepEqual([row.components,row.openTerminalLoops,row.euler,row.selfIntersections],[1,2,0,0]);
 assert.equal(proof.catalogueSha256,JSON.parse(fs.readFileSync('reports/ta98-male/biliary-junction-installed.json')).catalogueSha256);
 assert.ok(proof.literalReciprocalContacts.filter(r=>members.some(p=>p.id===r.a)&&members.some(p=>p.id===r.b)).length===3);
});
