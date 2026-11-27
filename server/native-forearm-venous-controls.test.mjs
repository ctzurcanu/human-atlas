import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {readFileSync,writeFileSync,existsSync} from 'node:fs';
import {build} from 'esbuild';
import test from 'node:test';

const base='.local-models/geometry-workbench/circulatory-completion/male-native-venous-regions-v1/literal-forearm-pack-v1';
const candidatePath=`${base}/proposed-atlas.json`;
const sourcePath='.local-models/geometry-workbench/circulatory-completion/male-cystic-artery-correction-proposal-v43/proposed-atlas.json';
const currentPath='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
const bodyPath='.local-models/organ-models/spinal-cord/male-body.json';
const sourceSha='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69';
const candidateSha='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de';
const hash=b=>createHash('sha256').update(b).digest('hex');
if(!existsSync(sourcePath)||!existsSync(candidatePath)){
 test('Native forearm venous controls require private source geometry',{skip:'Local licensed source and candidate catalogues are intentionally absent'},()=>{});
}else test('Native forearm venous candidate and installed ordinary-model controls',async()=>{
const sourceBytes=readFileSync(sourcePath),candidateBytes=readFileSync(candidatePath),atlas=JSON.parse(candidateBytes),source=JSON.parse(sourceBytes);
assert.equal(hash(sourceBytes),sourceSha,'Immutable pre-install source must remain bound');
assert.equal(hash(candidateBytes),candidateSha,'Candidate must match the reviewed installed revision');
const built=await build({stdin:{contents:`
 export {registeredNativeForearmVenousRegions,withNativeForearmVenousEntries} from './app/native-forearm-venous-controls.ts';
 export {anatomicalComponentIds} from './app/anatomical-components.ts';
 export {hierarchyEntries,systemPathFor} from './app/anatomy-hierarchy.ts';
 export {buildSystemTree} from './app/system-tree.tsx';
`,resolveDir:process.cwd(),sourcefile:'native-forearm-control-check.ts'},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import(`data:text/javascript;base64,${Buffer.from(built.outputFiles[0].contents).toString('base64')}`);
const registered=api.registeredNativeForearmVenousRegions(atlas);assert.equal(registered?.length,4);
const parts=new Map(atlas.parts.map(p=>[p.id,p])),owned=new Set(registered.flatMap(r=>[r.remainder.id,r.region.id]));
for(const r of registered){
 assert.equal(r.regionConcept.partOf,r.whole.id);
 assert.deepEqual(new Set(api.anatomicalComponentIds(r.whole,[],parts,atlas.concepts)),new Set([r.remainder.id,r.region.id]),'Direct Show components retains both disjoint literal source portions');
}
const entries=api.hierarchyEntries(atlas),scoped=api.withNativeForearmVenousEntries(atlas,entries);
assert.notEqual(entries,scoped);assert.equal(entries.length,scoped.length);
for(let i=0;i<entries.length;i++)if(!entries[i].parts.some(p=>owned.has(p.id)))assert.equal(entries[i],scoped[i],'Unrelated Systems entries remain the same objects');
const tree=api.buildSystemTree(entries,'systems',atlas),circulatory=tree.find(b=>b.id==='circulatory');assert.ok(circulatory);
const groups=[];const leaves=[];
const visit=nodes=>{for(const n of nodes){if(n.kind==='group'){if(['Cephalic vein','Basilic vein'].includes(n.name)&&n.parts.every(p=>owned.has(p.id)))groups.push(n);visit(n.nodes);}else if(n.kind==='bilateral')leaves.push(...n.entries);else leaves.push(n.entry);}};
visit(circulatory.nodes);
assert.equal(groups.length,2,'One whole assembly per named vein, each including both sides');
for(const name of ['Cephalic vein','Basilic vein']){
 const g=groups.find(g=>g.name===name);assert.equal(g.parts.length,4);assert.equal(new Set(g.parts.map(p=>p.id)).size,4);
 assert.ok(g.parts.every(p=>owned.has(p.id)));assert.equal(g.nodes.length,2,'Forearm region and remaining source course are the two child controls');
 assert.equal(g.nodes.filter(n=>n.name==='Other source course').length,1);
}
const selected=leaves.filter(e=>e.parts.some(p=>owned.has(p.id)));
assert.equal(selected.length,8);assert.deepEqual(selected.flatMap(e=>e.parts.map(p=>p.id)).sort(),[...owned].sort(),'Every physical portion occurs exactly once as a leaf');
for(const term of ['A12.3.08.021','A12.3.08.022'])assert.equal(selected.filter(e=>e.ta98Term===term).length,2,'Both exact side-qualified forearm references are discoverable');
assert.equal(selected.filter(e=>e.name.startsWith('Other source course')).length,4);
const noMarker={...atlas};delete noMarker.nativeVenousRegionIntegration;
assert.equal(api.withNativeForearmVenousEntries(noMarker,entries),entries,'Marker absence preserves legacy behavior');
const without=api.buildSystemTree(entries,'systems',noMarker),unrelated=nodes=>nodes.flatMap(n=>{
 const ids=n.parts.filter(p=>!owned.has(p.id)).map(p=>p.id);if(!ids.length)return [];
 return [{id:n.id,name:n.name,kind:n.kind,parts:ids,...(n.kind==='group'?{nodes:unrelated(n.nodes)}:{})}];
});
assert.deepEqual(unrelated(circulatory.nodes),unrelated(without.find(b=>b.id==='circulatory').nodes),'Existing accessory-cephalic and other unrelated ancestry groups remain literal');
assert.deepEqual(api.buildSystemTree(entries,'regions',atlas),api.buildSystemTree(entries,'regions',noMarker),'Regional hierarchy is outside this Systems-only presentation correction');
let invalid=0;
for(const mutate of [
 a=>a.parts.find(p=>p.id===registered[0].region.id).sourceOffset=[0,0,0],
 a=>a.parts.find(p=>p.id===registered[0].region.id).ta98Term='A12.3.08.022',
 a=>a.parts.find(p=>p.id===registered[0].region.id).sourceOuterSurfaceOnly=false,
 a=>a.parts.find(p=>p.id===registered[0].region.id).pairedTissue=true,
 a=>a.concepts.find(c=>c.id===registered[0].regionConcept.id).partOf=registered[1].whole.id,
 a=>a.concepts.find(c=>c.id===registered[0].whole.id).elements=[registered[0].remainder.id],
 a=>a.nativeVenousRegionIntegration.groups.pop(),
 a=>a.nativeVenousRegionIntegration.groups[0].sharedInterfaces[0].added_cap_triangles=1,
 a=>a.nativeVenousRegionIntegration.bankSha256='0'.repeat(64),
]){const a=structuredClone(atlas);mutate(a);assert.equal(api.registeredNativeForearmVenousRegions(a),undefined);invalid++;}
assert.equal(api.withNativeForearmVenousEntries(source,api.hierarchyEntries(source)).length,api.hierarchyEntries(source).length);
assert.equal(hash(readFileSync(sourcePath)),hash(sourceBytes));assert.equal(hash(readFileSync(candidatePath)),hash(candidateBytes));
const output={schema:'native-forearm-venous-controls-proof/v1',passed:true,source_catalogue_sha256:hash(sourceBytes),candidate_sha256:hash(candidateBytes),
 exact_new_terms:2,side_qualified_forearm_controls:4,whole_systems_groups:2,physical_leaf_controls:8,components_per_whole_side:2,
 unrelated_entry_objects_preserved:true,marker_absence_legacy_behavior_preserved:true,invalid_candidate_cases_rejected:invalid,
 source_catalogue_unchanged:true,installed:false,anatomical_acceptance:false};
const out=`${base}/controls-proof.json`;
if(existsSync(out))assert.deepEqual(JSON.parse(readFileSync(out)),output,'A repeated check must reproduce the sealed controls receipt');
let installed;
if(existsSync(currentPath)&&hash(readFileSync(currentPath))===candidateSha){
 assert.ok(existsSync(bodyPath),'The installed ordinary body catalogue is required');
 const bodyBytes=readFileSync(bodyPath),body=JSON.parse(bodyBytes),actual=api.registeredNativeForearmVenousRegions(body);
 assert.equal(actual?.length,4,'All four forearm assemblies must survive ordinary-model publication');
 assert.equal(body.sourceChunkCount,atlas.chunks.length);
 assert.equal(body.sourceChunkIndices.length,body.chunks.length);
 assert.equal(new Set(body.sourceChunkIndices).size,body.chunks.length);
 const ordinaryParts=new Map(body.parts.map(p=>[p.id,p]));
 for(const r of actual){
  for(const p of [r.remainder,r.region]){
   assert.deepEqual({...p,chunk:body.sourceChunkIndices[p.chunk]},parts.get(p.id),'Actual ordinary records must match the installed source, including array offsets');
   assert.deepEqual(body.chunks[p.chunk],atlas.chunks[body.sourceChunkIndices[p.chunk]]);
  }
  assert.deepEqual(new Set(api.anatomicalComponentIds(r.whole,[],ordinaryParts,body.concepts)),new Set([r.remainder.id,r.region.id]),'Ordinary-model Show components retains both disjoint portions');
 }
 const ordinaryTree=api.buildSystemTree(api.hierarchyEntries(body),'systems',body).find(b=>b.id==='circulatory');
 const ordinaryGroups=[],ordinaryLeaves=[];
 const scan=nodes=>{for(const n of nodes){if(n.kind==='group'){
  if(['Cephalic vein','Basilic vein'].includes(n.name)&&n.parts.every(p=>owned.has(p.id)))ordinaryGroups.push(n);
  scan(n.nodes);
 }else if(n.kind==='bilateral')ordinaryLeaves.push(...n.entries);else ordinaryLeaves.push(n.entry);}};
 scan(ordinaryTree.nodes);
 assert.equal(ordinaryGroups.length,2);
 assert.deepEqual(ordinaryLeaves.filter(e=>e.parts.some(p=>owned.has(p.id))).flatMap(e=>e.parts.map(p=>p.id)).sort(),[...owned].sort(),'Installed leaf controls select each actual source portion once');
 for(const term of ['A12.3.08.021','A12.3.08.022'])assert.equal(ordinaryLeaves.filter(e=>e.ta98Term===term).length,2);
 const bankUrls=new Set(actual.flatMap(r=>[r.remainder,r.region]).map(p=>body.chunks[p.chunk].url));assert.equal(bankUrls.size,1);
 const url=[...bankUrls][0];assert.ok(url.startsWith('/local-models/'));
 const bankSha=hash(readFileSync('.local-models/'+url.slice('/local-models/'.length)));
 assert.equal(bankSha,body.nativeVenousRegionIntegration.bankSha256);
 installed={schema:'installed-native-forearm-venous-controls-proof/v1',passed:true,immutable_source_sha256:sourceSha,installed_catalogue_sha256:candidateSha,
  ordinary_body_sha256:hash(bodyBytes),actual_bank_sha256:bankSha,registered_assemblies:4,whole_systems_groups:2,physical_leaf_controls:8,
  exact_new_terms:['A12.3.08.021','A12.3.08.022'],source_records_and_chunk_remap_literal:true,components_per_whole_side:2,
  installed:true,anatomical_acceptance:false};
 const installedOut=`${base}/installed-controls-proof.json`;
 if(existsSync(installedOut))assert.deepEqual(JSON.parse(readFileSync(installedOut)),installed,'Installed receipt must replay unchanged');
 else writeFileSync(installedOut,JSON.stringify(installed,null,2)+'\n');
 assert.equal(hash(readFileSync(bodyPath)),hash(bodyBytes),'Ordinary model must remain unchanged throughout the check');
 assert.equal(hash(readFileSync(currentPath)),candidateSha);
}
assert.equal(hash(readFileSync(sourcePath)),sourceSha);assert.equal(hash(readFileSync(candidatePath)),candidateSha);
console.log(JSON.stringify({candidate_controls:output,installed_controls:installed??'Current male source is a different revision; immutable candidate regression still passed'}));
});
