import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from './atlas-test-build.mjs';

const built=await build({stdin:{contents:"export * from './app/intercostal-mouth-relations'; export * from './app/local-ta98-layer-catalogue'; export {anatomicalRelations} from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,write:false,platform:'node',format:'esm'});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const route=api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('blood=biodigital-review'),true);
const atlas=JSON.parse(fs.readFileSync('.local-models'+route.slice('/local-models'.length)));
const prefix='BDM:male:human_20_male_cardiovascular_system-';

test('Installed ninth intercostal branch mouths are reciprocal and side-specific',()=>{
 const result=api.withVerifiedIntercostalMouths(atlas);
 for(const side of ['left','right']){
  const parent=result.parts.find(p=>p.id===prefix+side+'_posterior_intercostal_artery_9_00_ID');
  for(const [kind,suffix] of [['dorsal','9_0_ID'],['collateral','9_00_ID']]){
   const branch=result.parts.find(p=>p.id===prefix+kind+'_branch_of_'+side+'_posterior_intercostal_artery_'+suffix);
   assert.ok(parent&&branch);
   assert.ok(api.intercostalMouthRelations(parent,result.parts).some(r=>r.target.id===branch.id));
   assert.ok(api.intercostalMouthRelations(branch,result.parts).some(r=>r.target.id===parent.id));
   assert.ok(api.anatomicalRelations(branch,result.parts).some(r=>r.kind==='continuous'&&r.target.id===parent.id));
  }
  assert.equal(api.intercostalMouthRelations(parent,result.parts).length,2);
 }
});

test('All installed third through eleventh intercostal openings expose reciprocal connections',()=>{
 const result=api.withVerifiedIntercostalMouths(atlas);
 for(const side of ['left','right'])for(let level=3;level<=11;level++){
  const parent=result.parts.find(p=>p.id===prefix+side+'_posterior_intercostal_artery_'+level+'_00_ID');
  assert.ok(parent);
  for(const [kind,suffix] of [['dorsal','0_ID'],['collateral','00_ID']]){
   const branch=result.parts.find(p=>p.id===prefix+kind+'_branch_of_'+side+'_posterior_intercostal_artery_'+level+'_'+suffix);
   assert.ok(branch);
   assert.ok(api.intercostalMouthRelations(parent,result.parts).some(r=>r.target.id===branch.id));
   assert.ok(api.intercostalMouthRelations(branch,result.parts).some(r=>r.target.id===parent.id));
  }
 }
});

test('Installed dorsal intercostal openings link to the spinal branches on the same side',()=>{
 const result=api.withVerifiedIntercostalMouths(atlas);
 for(const side of ['left','right'])for(let level=3;level<=11;level++){
  const branch=result.parts.find(p=>p.id===prefix+'dorsal_branch_of_'+side+'_posterior_intercostal_artery_'+level+'_0_ID');
  const spinal=result.parts.find(p=>p.id===prefix+side+'_spinal_arterial_branches_ID');
  const opposite=result.parts.find(p=>p.id===prefix+(side==='left'?'right':'left')+'_spinal_arterial_branches_ID');
  assert.ok(branch&&spinal&&opposite);
  assert.ok(api.intercostalMouthRelations(branch,result.parts).some(r=>r.target.id===spinal.id));
  assert.ok(api.intercostalMouthRelations(spinal,result.parts).some(r=>r.target.id===branch.id));
  assert.ok(!api.intercostalMouthRelations(branch,result.parts).some(r=>r.target.id===opposite.id));
  assert.ok(api.anatomicalRelations(branch,result.parts).some(r=>r.kind==='continuous'&&r.target.id===spinal.id));
 }
});

test('Changing a bank or placement revokes the corresponding connection',()=>{
 const id=prefix+'collateral_branch_of_left_posterior_intercostal_artery_9_00_ID';
 const branch=atlas.parts.find(p=>p.id===id);
 const variants=[{...atlas,chunks:atlas.chunks.map((c,i)=>i===branch.chunk?{...c,sha256:'changed'}:c)}, {...atlas,parts:atlas.parts.map(p=>p.id===id?{...p,geometryOrigin:[0,0,0]}:p)}];
 for(const variant of variants){const result=api.withVerifiedIntercostalMouths(variant);const item=result.parts.find(p=>p.id===id);assert.equal(api.intercostalMouthRelations(item,result.parts).length,0);}
});
