import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/ileocaecal-reference-relations';export {anatomicalRelations} from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {ileocaecalReferenceRelations,ileocaecalReferenceLimitation,ileocaecalReferenceDescription,anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const atlas=JSON.parse(await fs.readFile('.local-models/ta98-female.json'));const valve=atlas.parts.find(p=>p.id==='HRA:VH_F_ileocecal_valve');const ids=['HRA:VH_F_caecum','HRA:VH_F_ileum_terminal','HRA:VH_F_vermiform_appendix'];
test('Original HRA junction sources have reciprocal associations without claiming a tissue join',()=>{
 const graph=ileocaecalReferenceRelations(atlas.parts);assert.deepEqual(graph.get(valve.id).map(r=>r.target.id).sort(),ids.toSorted());
 for(const r of graph.get(valve.id)){assert.equal(r.kind,'adjacent');assert.ok(graph.get(r.target.id).some(x=>x.target===valve&&x.kind===r.kind&&x.note===r.note));assert.ok(anatomicalRelations(valve,atlas.parts).some(x=>x.target.id===r.target.id&&x.kind==='adjacent'&&x.note===r.note));}
 assert.match(ileocaecalReferenceDescription(valve).text,/muscular and mucosal tissues/);assert.equal(graph.get(valve.id).find(r=>r.target.id===ids[2]).viaModeled,true);assert.match(ileocaecalReferenceLimitation(valve),/open annular sheet/);assert.match(ileocaecalReferenceLimitation(valve),/lips, frenula or wall layers/);
});
test('Suppressed, independently moved, native and unqualified sources do not receive donor associations',()=>{
 const caecum=atlas.parts.find(p=>p.id===ids[0]);
 for(const p of [{...caecum,suppressed:true},{...caecum,reviewSource:'other'},{...caecum,sourceId:'other'},{...caecum,sourceOffset:[0,0,0]},{...caecum,geometryOrigin:[0,0,0]},{...caecum,id:'LOCAL:female:_large_intestine_cecum:77'}])assert.equal(ileocaecalReferenceRelations([valve,p]).size,0);
 for(const p of [{...valve,suppressed:true},{...valve,reviewSource:'other'},{...valve,sourceId:'other'},{...valve,system:'arterial'}]){assert.equal(ileocaecalReferenceRelations([p,caecum]).size,0);assert.equal(ileocaecalReferenceLimitation(p),undefined);}
});
test('Missing caecum remains an explicit unmodeled neighborhood intermediary',()=>{const appendix=atlas.parts.find(p=>p.id===ids[2]);const r=ileocaecalReferenceRelations([valve,appendix]).get(valve.id)[0];assert.deepEqual(r.via,['Caecum']);assert.equal(r.viaModeled,false);});
test('Complete snapshot cache does not retain removed source neighbors',()=>{
 const one=[valve,atlas.parts.find(p=>p.id===ids[0])],two=[valve,atlas.parts.find(p=>p.id===ids[1])];assert.equal(ileocaecalReferenceRelations(one).get(valve.id)[0].target.id,ids[0]);assert.equal(ileocaecalReferenceRelations(two).get(valve.id)[0].target.id,ids[1]);assert.equal(ileocaecalReferenceRelations(one),ileocaecalReferenceRelations(one));
});
test('Valve reference receives only broader caecal ancestry; no finished TA98 valve tissue is invented',async()=>{
 const metadata=JSON.parse(await fs.readFile('app/data/ta98-metadata.json'));assert.deepEqual(metadata.reviewedCatalogueIdentities[valve.id],{term:'A05.7.02.001',kind:'parent'});const rules=JSON.parse(await fs.readFile('scripts/ta98-reviewed-mappings.json'));const rule=rules.find(r=>r.conceptIds?.includes(valve.id));assert.equal(rule.kind,'parent');assert.match(rule.note,/no exact papilla/);
});
