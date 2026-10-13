import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const read=p=>JSON.parse(readFileSync(p));
const report=read('reports/ta98-female/installed-biliary-source-scope.json');
const pid='LOCAL:female:bile_duct:1031',cid='LOCAL:female:bile_duct';
const built=await build({stdin:{contents:"export * from './app/anatomical-relations'; export * from './app/guest-hierarchy'; export * from './app/reviewed-catalogue-identities'; export * from './app/viewer-state'; export * from './app/anatomical-terminology'; export * from './app/digestive-duct-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const hierarchy=read('public/assets/ta98.json');

test('all four catalogues correct scope without changing any geometry resource or unrelated part',()=>{
 assert.equal(report.catalogues.length,4);
 for(const row of report.catalogues){
  const raw=readFileSync(row.path),a=JSON.parse(raw),before=JSON.parse(gunzipSync(readFileSync(row.snapshot)));
  assert.equal(createHash('sha256').update(raw).digest('hex'),row.afterSha256);
  assert.deepEqual(gunzipSync(readFileSync(row.path+'.gz')),raw);
  assert.deepEqual(a.chunks,before.chunks);assert.deepEqual(a.materials,before.materials);
  assert.deepEqual(a.parts.map(p=>p.id===pid?{...p,name:'Bile duct'}:p),before.parts);
  assert.deepEqual(a.concepts.filter(c=>![cid,'TA98FEMALE:A05.8.02.013'].includes(c.id)),before.concepts.filter(c=>![cid,'TA98FEMALE:A05.8.02.013'].includes(c.id)));
  assert.equal(a.concepts.find(c=>c.id===cid).ta98Kind,'parent');
  assert.equal(a.concepts.find(c=>c.id===cid).ta98Term,'A05.0.00.000');
  assert.deepEqual(api.resolveSelection(a,[cid]),[pid]);
 }
});

test('combined source is accessible under Alimentary system without a false common-bile-duct or pancreatic identity',()=>{
 const a=api.withReviewedCatalogueIdentities(read('.local-models/ta98-female-runtime.json'));
 const roots=api.resolveGuestHierarchy(a,hierarchy),nodes=new Map();
 const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};roots.forEach(walk);
 assert.ok(nodes.get('A05.0.00.000').parts.some(p=>p.id===pid));
 for(const code of ['A05.8.02.013','A05.9.01.015','A05.8.02.001'])assert.ok(!nodes.get(code).parts.some(p=>p.id===pid),code);
 assert.equal(api.taExactEntityForConcept(cid),undefined);
 assert.equal(api.terminologyForConcept(cid,'A05.0.00.000','parent').ta98,null);
 // Stale embedded source tags cannot reinstate the disproven exact identity.
 const stale={...a,concepts:a.concepts.map(c=>c.id===cid?{...c,ta98Term:'A05.8.02.013',ta98Kind:'exact'}:c)};
 assert.equal(api.withReviewedCatalogueIdentities(stale).concepts.find(c=>c.id===cid).ta98Kind,'parent');
});

test('source-qualified hepatic and pancreatic routes are reciprocal and preserve separate organ identity',()=>{
 const a=read('.local-models/ta98-female-runtime.json'),duct=a.parts.find(p=>p.id===pid);
 for(const id of ['LOCAL:female:liver:1032','LOCAL:female:pancreas:1039']){
  const organ=a.parts.find(p=>p.id===id);
  assert.ok(api.anatomicalRelations(organ,a.parts).some(r=>r.kind==='after'&&r.target.id===pid));
  assert.ok(api.anatomicalRelations(duct,a.parts).some(r=>r.kind==='before'&&r.target.id===id));
  assert.ok(!api.anatomicalRelations(organ,a.parts).some(r=>r.kind==='contains'&&r.target.id===pid));
  assert.equal(api.digestiveDuctRelations([duct,{...organ,sourceId:'donor'}]).size,0);
 }
 assert.ok(api.anatomicalRelations(duct,a.parts).some(r=>r.kind==='after'&&r.target.id==='LOCAL:female:intestine_duodenum:183'));
 assert.match(api.digestiveDuctSourceLimitation(duct),/combined native source/i);
});
