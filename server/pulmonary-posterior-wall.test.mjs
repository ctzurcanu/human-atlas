import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/pulmonary-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {pulmonaryRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const fullBundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(fullBundle.outputFiles[0].contents).toString('base64'));
const atlas=JSON.parse(await readFile('.local-models/organ-models/spinal-cord/male-body.json','utf8'));
const id='TA98PART:ZA:Posterior segmental artery of right lung:A12.2.01.108';
test('Installed finite posterior branch retains reciprocal flow, territory and airway associations',()=>{
 const p=atlas.parts.find(p=>p.id===id);assert.ok(p&&!p.suppressed&&p.sectionSolid&&p.pairedTissue);
 const graph=pulmonaryRelations(atlas.parts),links=graph.get(id);
 for(const [target,kind,reverse] of [['ZA:Right pulmonary artery','before','after'],['TA98PULMONARY:male:right-superior:posterior','supplies','arterial'],['LOCAL:male:pulmonary:bronchus:right_posterior_bronchus','adjacent','adjacent']]){
  assert.ok(links.some(r=>r.kind===kind&&r.target.id===target),target);
  assert.ok(graph.get(target).some(r=>r.kind===reverse&&r.target.id===id),target);
 }
 assert.match(links.find(r=>r.kind==='before').note,/junction remains unverified/);
 assert.deepEqual(anatomicalRelations(p,atlas.parts).filter(r=>r.kind==='before').map(r=>r.target.id),['ZA:Right pulmonary artery']);
});
test('Suppressed sources, opposite bodies and unqualified lookalikes cannot inherit the branch links',()=>{
 const p=atlas.parts.find(p=>p.id===id);
 for(const bad of [{...p,suppressed:true},{...p,id:'LOCAL:female:lookalike'}, {...p,sourceId:'Other artery'}, {...p,tissueWallAssembly:'other'}]){
  const parts=atlas.parts.map(q=>q.id===id?bad:q);
  assert.equal(pulmonaryRelations(parts).get(bad.id)?.some(r=>r.kind==='supplies')??false,false);
 }
});
test('Installed anterior descending branch keeps its pulmonary root and territory, excluding the hepatic namesake',()=>{
 const anteriorId='TA98PART:ZA:Anterior segmental artery of right lung:A12.2.01.106';
 const p=atlas.parts.find(p=>p.id===anteriorId);assert.ok(p&&!p.suppressed&&p.pairedTissue&&p.sectionSolid);
 const graph=pulmonaryRelations(atlas.parts),links=anatomicalRelations(p,atlas.parts);
 assert.deepEqual(links.filter(r=>r.kind==='before').map(r=>r.target.id),['ZA:Right pulmonary artery']);
 assert.ok(links.some(r=>r.kind==='supplies'&&r.target.id==='TA98PULMONARY:male:right-superior:anterior'));
 assert.ok(graph.get('LOCAL:male:pulmonary:bronchus:right_anterior_bronchus').some(r=>r.kind==='adjacent'&&r.target.id===anteriorId));
 assert.match(links.find(r=>r.kind==='before').note,/proximal conduit.*not constructed/);
 const changed=atlas.parts.map(q=>q.id===anteriorId?{...q,sourceId:'Unregistered donor'}:q);
 assert.equal(pulmonaryRelations(changed).get(anteriorId)?.some(r=>r.kind==='supplies')??false,false);
});
