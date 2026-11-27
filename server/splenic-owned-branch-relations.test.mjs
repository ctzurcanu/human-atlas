import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/foregut-arterial-relations';export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const male=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
const trunkId='TA98PART:ZA:Splenic artery:trunk',branchId='TA98PART:ZA:Splenic artery:A12.2.12.051';
function fixture(){
 const base=male.parts.find(p=>p.id==='ZA:Splenic artery');
 const cell=(role,id,term,name,neighbor)=>({...base,id,conceptId:id,name,ta98Term:term,suppressed:false,sourceId:'Splenic artery',chunk:99,geometryOrigin:[0,1.25,0],pairedTissue:true,sectionSolid:true,indexCount:300,tissueWallAssembly:'splenic-artery:male:owned-branches',splenicBranchOwnership:{schema:'source-bound-hilar-tree/v1',role,boundaryProvisional:true,organInterfaceVerified:false},tissueDrawGroups:[{start:0,count:228},{start:228,count:72,neighbor}]});
 const trunk=cell('trunk',trunkId,'A12.2.12.040','Splenic arterial trunk',branchId),branches=cell('branches',branchId,'A12.2.12.051','Splenic branches',trunkId);
 const retired={...base,suppressed:true,replacementParts:[trunkId,branchId]},coeliac=male.parts.find(p=>p.id==='ZA:Coeliac trunk');
 const tissues=['serosa','capsule','pulp','trabeculae'].map(role=>male.parts.find(p=>p.id===`TA98PART:male:Spleen:${role}`));
 return {trunk,branches,retired,coeliac,tissues,parts:[retired,trunk,branches,coeliac,...tissues]};
}
test('owned cells expose reciprocal shared-annulus and qualified coeliac parent links',()=>{
 const {trunk,branches,coeliac,parts}=fixture(),links=p=>api.anatomicalRelations(p,parts);
 assert.ok(links(trunk).some(r=>r.kind==='after'&&r.target===branches&&/reciprocal annular/.test(r.note)));
 assert.ok(links(branches).some(r=>r.kind==='before'&&r.target===trunk));
 assert.ok(links(trunk).some(r=>r.kind==='before'&&r.target===coeliac&&/junction remains unverified/.test(r.note)));
 assert.ok(links(coeliac).some(r=>r.kind==='after'&&r.target===trunk));
 assert.ok(!links(coeliac).some(r=>r.target===branches||r.target.id==='ZA:Splenic artery'));
 assert.ok(!links(trunk).some(r=>r.kind==='supplies'));
});
test('supply targets the exact whole male spleen and preserves recipient/perfusion qualification',()=>{
 const {branches,tissues,parts}=fixture(),whole={id:'ZA:Spleen',name:'Spleen',elements:tissues.map(p=>p.id)};
 const supply=api.anatomicalRelations(branches,parts).filter(r=>r.kind==='supplies');
 assert.equal(supply.length,1);assert.equal(supply[0].target,tissues[0]);assert.deepEqual(supply[0].targetConcept,whole);
 assert.match(supply[0].note,/Recipient tissue entry, biological perfusion.*unverified/);
 assert.deepEqual(api.anatomicalConceptFunctionalRelations(whole,tissues,parts).map(r=>[r.kind,r.target.id]),[['arterial',branchId]]);
 assert.ok(api.anatomicalRelations(tissues[0],parts).some(r=>r.kind==='arterial'&&r.target===branches));
 assert.deepEqual(api.foregutWholeFunctionalRelations({...whole,id:'selection-set'},tissues,parts),undefined);
 assert.deepEqual(api.foregutWholeFunctionalRelations(whole,tissues.slice(1),parts),[]);
});
test('exact IDs, schema, role, pose, source and reciprocal rim metadata are mandatory',()=>{
 for(const patch of [{id:'HRA:donor:splenic-branches'},{sourceId:'donor'},{suppressed:true},{pairedTissue:false},{sectionSolid:false},{normalSigns:undefined},{tissueWallAssembly:'other'},{geometryOrigin:[0,0,0]},{sourceOffset:[0.001,0,0]},{geometryScale:2},{geometryInstance:{rotation:[0,0,0,1],scale:1}},{system:'venous'},{license:'unknown'},{splenicBranchOwnership:{schema:'unreviewed',role:'branches'}},{splenicBranchOwnership:{schema:'source-bound-hilar-tree/v1',role:'trunk'}},{tissueDrawGroups:[]}]){
  const f=fixture(),bad={...f.branches,...patch},parts=f.parts.map(p=>p===f.branches?bad:p);
  assert.ok(!api.anatomicalRelations(f.trunk,parts).some(r=>r.target===bad),JSON.stringify(patch));
  assert.ok(!api.anatomicalRelations(bad,parts).some(r=>r.target===f.trunk||r.kind==='supplies'),JSON.stringify(patch));
 }
});
test('generic, retired, female and donor peers cannot replace registered recipient or parent cells',()=>{
 const f=fixture();
 for(const bad of [{...f.coeliac,id:'HRA:donor:coeliac'},{...f.coeliac,suppressed:true},{...f.coeliac,sourceId:'unrelated'},{...f.coeliac,sourceOffset:[0.01,0,0]}]){
  const parts=f.parts.map(p=>p===f.coeliac?bad:p);
  assert.ok(!api.anatomicalRelations(f.trunk,parts).some(r=>r.target===bad));
 }
 for(const patch of [{id:'LOCAL:female:spleen:serosa'},{suppressed:true},{sourceId:'donor'},{geometryOrigin:[0,0,0]},{tissueWallAssembly:'splenic:female'}]){
  const bad={...f.tissues[0],...patch},parts=f.parts.map(p=>p===f.tissues[0]?bad:p);
  assert.ok(!api.anatomicalRelations(f.branches,parts).some(r=>r.kind==='supplies'));
 }
});
