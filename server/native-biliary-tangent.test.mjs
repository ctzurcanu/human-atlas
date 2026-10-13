import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
import {build} from 'esbuild';
const result=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/biliary-junction';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(result.outputFiles[0].contents).toString('base64'));
const path='.local-models/ta98-review/female-native-biliary-tangent/atlas.json',payload=fs.readFileSync(path),atlas=JSON.parse(payload),hash=crypto.createHash('sha256').update(payload).digest('hex');
test('whole-duct tangent review remains private and retains the original anatomical selectors',()=>{
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('biliary=native-tangent'),true),'/local-models/ta98-review/female-native-biliary-tangent/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('biliary=native-tangent'),false),undefined);
 assert.ok(!api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('biliary=native-tangent'),true).includes('tangent'));
 assert.deepEqual(atlas.concepts.find(c=>c.id==='LOCAL:female:bile_duct').elements,['LOCAL:female:bile_duct:1031']);assert.equal(atlas.commonHepaticReview,undefined);
});
test('the cystic communication remains literal and the patch does not invent new TA98 coverage',()=>{
 const duct=atlas.parts.find(p=>p.id==='LOCAL:female:bile_duct:1031'),cy=atlas.parts.find(p=>p.id==='LOCAL:female:cystic_duct:authored');assert.ok(api.biliaryJunctionInterface(cy,duct));
 const replay=JSON.parse(fs.readFileSync('reports/ta98-female/native-biliary-tangent-atlas-validation.json'));assert.equal(replay.catalogueSha256,hash);assert.equal(replay.retainedGaCells,16);assert.equal(replay.namedTermsAdded,0);assert.equal(replay.unrelatedPartsExact,true);assert.equal(replay.conceptsMaterialsVisibilityExact,true);assert.equal(replay.cohortTrianglesAfter-replay.cohortTrianglesBefore,448);
});
test('these exact packed bytes retain every communicating seam and add no primary body contacts',()=>{
 const macro=JSON.parse(fs.readFileSync('reports/ta98-female/native-biliary-tangent-packed-validation.json')),packing=JSON.parse(fs.readFileSync('reports/ta98-female/native-biliary-tangent-packing.json'));
 assert.equal(macro.catalogueSha256,hash);assert.equal(packing.catalogueSha256,hash);assert.equal(macro.literalReciprocalContacts.length,30);assert.equal(macro.opposedSignedNormalsExact,true);assert.equal(macro.continuousMacroLumen,true);assert.equal(packing.atlasPacked,true);assert.equal(packing.zeroNewContacts,true);
 for(const cell of macro.finiteCells){assert.equal(cell.topology.openEdges,0);assert.equal(cell.selfIntersections.properOrCoplanarIntersections,0);}
});
