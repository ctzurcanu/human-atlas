import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import vm from 'node:vm';
import {gunzipSync} from 'node:zlib';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/foregut-arterial-relations'; export * from './app/anatomical-relations'; export * from './app/guest-hierarchy'; export * from './app/anatomical-components'; export * from './app/ta98-modeled-parent';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const read=path=>JSON.parse(fs.readFileSync(path));
const female=read('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json');
test('installation preserves every original source attribute and compressed catalogue, adding no wall claims',()=>{
 const report=read('reports/ta98-female/celiac-native-source-identity.json');
 assert.equal(report.status,'installed-metadata-only');
 for(const row of report.records){
  assert.deepEqual(gunzipSync(fs.readFileSync(row.path+'.gz')),fs.readFileSync(row.path));
  const a=read(row.path),p=a.parts.find(p=>p.id===report.partId);
  assert.deepEqual({...p,system:row.originalAffectedRecords.part.system},row.originalAffectedRecords.part);
  assert.ok(!p.pairedTissue&&!p.sectionSolid);
 }
 assert.equal(report.addedTriangles,0);assert.equal(report.addedTextures,0);assert.equal(report.addedDrawGroups,0);
 assert.deepEqual(report.exactBranchInterfaces.map(x=>x.count),[6,5,7]);
});
test('all four installed female catalogues use the original celiac source under the arterial TA98 identity',()=>{
 for(const path of ['.local-models/female.json','.local-models/ta98-female.json','.local-models/ta98-female-runtime.json','.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json']){
  const a=read(path),p=a.parts.find(p=>p.id==='LOCAL:female:celiac_trunk:1756'),c=a.concepts.find(c=>c.id===p.conceptId);
  assert.equal(p.system,'arterial');assert.equal(p.sourceId,'celiac_trunk');assert.equal(p.indexCount,399);
  assert.equal(c.ta98Term,'A12.2.12.012');assert.equal(c.ta98Kind,'exact');assert.deepEqual(c.elements,[p.id]);
 }
 const tree=api.resolveGuestHierarchy(female,api.parseGuestHierarchy(read('public/assets/ta98.json')));let node;
 const walk=n=>{if(n.id.includes('A12.2.12.012'))node=n;n.children.forEach(walk);};tree.forEach(walk);
 assert.ok(node);assert.ok(node.parts.some(p=>p.id==='LOCAL:female:celiac_trunk:1756'));
});
test('the importer preserves arterial identity despite the misleading source ancestry',()=>{
 const text=fs.readFileSync('scripts/import-local-zygote.mjs','utf8'),start=text.indexOf('function systemFor('),end=text.indexOf('\nfunction decode(',start);
 const fn=vm.runInNewContext('('+text.slice(start,end)+')',{systems:{7:'arterial'}});
 assert.equal(fn('celiac_trunk',{layer:7,route:['Circulatory system','Veins']}),'arterial');
 assert.equal(fn('portal_vein',{layer:7,route:['Veins']}),'venous');
});
test('both native foregut trees have reciprocal branches without selecting independent donors',()=>{
 for(const sex of ['male','female']){
  const atlas=read(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`),p=atlas.parts.find(p=>p.id===(sex==='male'?'ZA:Coeliac trunk':'LOCAL:female:celiac_trunk:1756'));
  const links=api.anatomicalRelations(p,atlas.parts).filter(r=>r.kind==='after');
  for(const name of sex==='male'?['Left gastric artery','Common hepatic artery','Splenic artery']:['Gastric artery (left)','Common hepatic artery','Splenic artery']){
   const native=atlas.parts.find(q=>q.name===name&&q.id.startsWith(sex==='male'?'ZA:':'LOCAL:female:'));
   assert.ok(links.some(r=>r.target.id===native.id),name);
   const back=api.anatomicalRelations(native,atlas.parts).filter(r=>r.kind==='before');
   assert.ok(back.some(r=>r.target.id===p.id));assert.ok(back.every(r=>r.target.id.startsWith(sex==='male'?'ZA:':'LOCAL:female:')));
  }
  assert.ok(links.every(r=>r.target.id.startsWith(sex==='male'?'ZA:':'LOCAL:female:')));
  const choice=atlas.concepts.find(c=>c.id===p.conceptId),parts=new Map(atlas.parts.map(p=>[p.id,p]));
  assert.equal(api.ta98ModeledParent(choice,atlas.concepts,parts),undefined,'aortic ancestry is not tissue containment');
  assert.deepEqual(api.anatomicalComponentIds(choice,api.anatomicalRelations(p,atlas.parts),parts,atlas.concepts),[],'branch arteries are related anatomy, not trunk wall components');
 }
});
test('foregut branches reject wrong systems, other-sex sources and retired references in either direction',()=>{
 const a=female.parts.find(p=>p.id==='LOCAL:female:celiac_trunk:1756'),b=female.parts.find(p=>p.id==='LOCAL:female:common_hepatic_artery:1722');
 for(const bad of [{...b,id:'HRA:VH_F_common_hepatic_artery'},{...b,id:'LOCAL:male:common_hepatic_artery:1722'},{...b,suppressed:true},{...b,system:'venous'},{...b,sourceId:'unrelated'}]){
  assert.equal(api.foregutArterialFrameCompatible('after',a,bad),false);
  assert.equal(api.foregutArterialFrameCompatible('before',bad,a),false);
  assert.ok(!api.anatomicalRelations(a,[a,bad]).some(r=>r.kind==='after'&&r.target.id===bad.id));
  assert.ok(!api.anatomicalRelations(bad,[a,bad]).some(r=>r.kind==='before'&&r.target.id===a.id));
 }
});
