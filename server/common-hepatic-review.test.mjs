import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
import {build} from 'esbuild';
const result=await build({stdin:{contents:"export * from './app/digestive-duct-relations';export * from './app/biliary-junction';export * from './app/local-ta98-layer-catalogue';export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(result.outputFiles[0].contents).toString('base64'));
const path='.local-models/ta98-review/female-common-hepatic-wall/atlas.json',atlas=JSON.parse(fs.readFileSync(path)),hash=crypto.createHash('sha256').update(fs.readFileSync(path)).digest('hex');
const ids=atlas.commonHepaticReview.partIds,[distal,common,upper]=ids.map(id=>atlas.parts.find(p=>p.id===id));
test('TA98 row selects the finite common hepatic candidate instead of the older donor',()=>{
 const tree=JSON.parse(fs.readFileSync('public/assets/ta98.json'));
 const find=(nodes,id)=>nodes.flatMap(n=>n.id===id?[n]:find(n.children,id));
 const row=find(api.resolveGuestHierarchy(atlas,tree),'A05.8.01.061')[0];
 assert.deepEqual(row.directParts.map(p=>p.id),[common.id]);assert.ok(row.description.includes(common.coverageLimitation));assert.ok(!row.description.includes('2 exact concept(s)'));
 const displaced={...atlas,parts:atlas.parts.map(p=>p.id===common.id?{...p,sourceOffset:[.01,0,0]}:p)};
 assert.ok(!find(api.resolveGuestHierarchy(displaced,tree),'A05.8.01.061')[0].directParts.some(p=>p.id===common.id),'a displaced candidate cannot replace the donor selector');
 const ordinary=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
 assert.ok(!find(api.resolveGuestHierarchy(ordinary,tree),'A05.8.01.061')[0].directParts.some(p=>p.id===common.id));
});
test('review exposes one named finite common hepatic wall without claiming named tributaries',()=>{
 const c=atlas.concepts.find(c=>c.id===common.conceptId);assert.equal(c.ta98Term,'A05.8.01.061');assert.equal(c.ta98Kind,'exact');assert.deepEqual(c.elements,[common.id]);
 assert.ok(atlas.parts.find(p=>p.id==='LOCAL:female:bile_duct:1031').suppressed);assert.equal(atlas.commonHepaticReview.anatomicallyAccepted,false);
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('biliary=hepatic-wall'),true),'/local-models/ta98-review/female-common-hepatic-wall/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('biliary=hepatic-wall'),false),undefined);
 assert.ok(!api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('biliary=hepatic-wall'),true).includes('common-hepatic'));
});
test('actual shared walls expose both common hepatic neighbors and the one cystic confluence',()=>{
 const graph=api.digestiveDuctRelations(atlas.parts),relations=graph.get(common.id);
 for(const [kind,target] of [['before',upper],['after',distal],['continuous',upper],['continuous',distal]])assert.ok(relations.some(r=>r.kind===kind&&r.target.id===target.id),`${kind} ${target.name}`);
 const cystic=atlas.parts.find(p=>p.id==='LOCAL:female:cystic_duct:authored');assert.ok(api.biliaryJunctionInterface(cystic,distal));assert.equal(api.biliaryJunctionInterface(cystic,common),false);
 assert.equal(api.biliaryJunctionInterface({...common,sourceOffset:[.01,0,0]},distal),false);
 assert.equal(api.biliaryJunctionInterface({...common,sourceId:'foreign'},distal),false);
});
test('packed evidence covers these exact bytes and every changed finite domain',()=>{
 const proof=JSON.parse(fs.readFileSync('reports/ta98-female/common-hepatic-packed-validation.json'));assert.equal(proof.catalogueSha256,hash);assert.equal(proof.opposedSignedNormalsExact,true);assert.equal(proof.literalReciprocalContacts.length,32);assert.deepEqual(proof.finiteCells.map(c=>c.id),ids);
 for(const cell of proof.finiteCells){assert.equal(cell.topology.openEdges,0);assert.equal(cell.selfIntersections.properOrCoplanarIntersections,0);}
});
