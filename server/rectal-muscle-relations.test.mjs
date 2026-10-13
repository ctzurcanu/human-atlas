import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations,anatomicalConceptFunctionalRelations}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
// Exercise the active native-source contract independently of later wall
// installation. Actual replacement catalogs are checked below.
const female=JSON.parse(readFileSync('.local-models/ta98-female-runtime.json')).parts.map(p=>p.id==='LOCAL:female:large_intestine_rectum:80'?{...p,suppressed:false}:p);
const rectum=female.find(p=>p.id==='LOCAL:female:large_intestine_rectum:80');
const muscle=female.find(p=>p.id==='LOCAL:female:puborectalis:236');
assert.ok(rectum&&muscle);
test('actual female rectum and puborectalis expose reciprocal topography without containment',()=>{
 for(const [a,b] of [[rectum,muscle],[muscle,rectum]]){
  const rows=anatomicalRelations(a,female).filter(r=>r.target.id===b.id);
  assert.equal(rows.filter(r=>r.kind==='adjacent').length,1);
  assert.ok(!rows.some(r=>['partOf','contains','continuous'].includes(r.kind)));
  assert.match(rows.find(r=>r.kind==='adjacent').note,/anorectal junction/);
 }
});
test('native male framing works without borrowing female, donor, hidden or wrong-system muscles',()=>{
 const r={...rectum,id:'LOCAL:male:rectum',name:'Rectum'};
 const left={...muscle,id:'LOCAL:male:puborectalis-left',name:'Left puborectalis'};
 const right={...muscle,id:'LOCAL:male:puborectalis-right',name:'Puborectalis (right)'};
 const invalid=[muscle,{...left,id:'HRA:foreign-puborectalis'},{...right,id:'LOCAL:male:hidden',suppressed:true},{...left,id:'LOCAL:male:nerve',system:'nervous'}];
 const parts=[r,left,right,...invalid];
 assert.deepEqual(anatomicalRelations(r,parts).filter(x=>x.kind==='adjacent').map(x=>x.target.id).sort(),[left.id,right.id].sort());
 const donor={...r,id:'FJ2571'};
 assert.ok(!anatomicalRelations(donor,[donor,left,right]).some(x=>x.kind==='adjacent'));
});
test('new catalogue snapshots remove retired relationships',()=>{
 assert.ok(anatomicalRelations(rectum,[rectum,muscle]).some(x=>x.kind==='adjacent'));
 assert.ok(!anatomicalRelations(rectum,[rectum,{...muscle,suppressed:true}]).some(x=>x.kind==='adjacent'));
});
test('native female rectal continuation excludes the separately registered donor sigmoid',()=>{
 const sigmoid=female.find(p=>p.id==='LOCAL:female:sigmoid_colon:74');
 const donor=female.find(p=>p.id==='HRA:VH_F_sigmoid_colon');
 assert.ok(sigmoid&&donor);
 assert.deepEqual(anatomicalRelations(rectum,female).filter(r=>r.kind==='before').map(r=>r.target.id),[sigmoid.id]);
 assert.ok(!anatomicalRelations(donor,female).some(r=>r.kind==='after'&&r.target.id===rectum.id));
});
test('foreign collectors cannot become native rectal supply or drainage by matching a name',()=>{
 const artery=female.find(p=>p.name==='Superior rectal artery');
 const vein=female.find(p=>p.name==='Superior rectal vein');
 assert.ok(artery&&vein);
 const parts=[rectum,artery,vein,...[artery,vein].flatMap(p=>[{...p,id:`HRA:${p.id}`},{...p,id:p.id.replace('LOCAL:female:','LOCAL:male:')},{...p,id:p.id+':retired',suppressed:true}])];
 for(const kind of ['arterial','venous']){
  const links=anatomicalRelations(rectum,parts).filter(r=>r.kind===kind);
  assert.ok(links.length);assert.ok(links.every(r=>[artery.id,vein.id].includes(r.target.id)));
 }
});
test('actual installed female covering carries rectal route and sling links while retired source stays excluded',()=>{
 const atlas=JSON.parse(readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json')).parts;
 const retired=atlas.find(p=>p.id===rectum.id),source=atlas.find(p=>p.id==='LOCAL:female:CAECOWALL:rectum:serosa'),cover=atlas.find(p=>p.id==='LOCAL:female:CAECOWALL:sigmoid:serosa');assert.ok(retired?.suppressed&&source&&cover);
 assert.ok(!anatomicalRelations(retired,atlas).some(r=>['before','adjacent'].includes(r.kind)));
 const before=anatomicalRelations(source,atlas).filter(r=>r.kind==='before');assert.deepEqual(before.map(r=>r.target.id),[cover.id]);assert.match(before[0].note,/whole named structure/);
 assert.ok(anatomicalRelations(cover,atlas).some(r=>r.kind==='after'&&r.target.id===source.id));
 assert.ok(anatomicalRelations(source,atlas).some(r=>r.kind==='adjacent'&&r.target.id===muscle.id));
 assert.ok(anatomicalRelations(atlas.find(p=>p.id===muscle.id),atlas).some(r=>r.kind==='adjacent'&&r.target.id===source.id));
 for(const invalid of [{...cover,suppressed:true},{...cover,sourceOffset:[0,.01,0]},{...cover,bowelWallSource:'HRA:foreign-source'}]){
  const snapshot=atlas.map(p=>p.id===cover.id?invalid:p);assert.ok(!anatomicalRelations(source,snapshot).some(r=>r.kind==='before'&&r.target.id===cover.id));
 }
});
test('whole installed rectum has organ-level supply links; individual coats and incomplete assemblies do not',()=>{
 const atlas=JSON.parse(readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
 const ids=atlas.rectalWallIntegration.partIds,selected=ids.map(id=>atlas.parts.find(p=>p.id===id));
 const concept=atlas.concepts.find(c=>c.id==='LOCAL:female:large_intestine_rectum');assert.ok(concept);
 const links=anatomicalConceptFunctionalRelations(concept,selected,atlas.parts);
 assert.ok(links.some(r=>r.kind==='arterial'&&r.target.id==='LOCAL:female:superior_rectal_artery:1687'));
 assert.ok(links.some(r=>r.kind==='venous'&&r.target.id==='LOCAL:female:superior_rectal_vein:1789'));
 const hierarchyChoice={id:'hierarchy:guest:ta98:A05.7.04.001',name:'rectum',elements:concept.elements};
 assert.deepEqual(anatomicalConceptFunctionalRelations(hierarchyChoice,selected,atlas.parts).map(r=>[r.kind,r.target.id]),links.map(r=>[r.kind,r.target.id]));
 assert.ok(links.every(r=>!r.target.suppressed&&!r.target.id.startsWith('HRA:')));
 for(const subset of [selected.slice(1),[selected[0]]]){
  assert.ok(!anatomicalConceptFunctionalRelations({...concept,elements:subset.map(p=>p.id)},subset,atlas.parts).some(r=>['arterial','venous'].includes(r.kind)));
  assert.deepEqual(anatomicalConceptFunctionalRelations({...concept,name:'Rectum',elements:subset.map(p=>p.id)},subset,atlas.parts),[]);
 }
 for(const bad of [{...selected[1],sourceOffset:[0,.01,0]},{...selected[1],bowelWallSource:'HRA:foreign'}]){
  const parts=atlas.parts.map(p=>p.id===bad.id?bad:p),selection=selected.map(p=>p.id===bad.id?bad:p);
  assert.ok(!anatomicalConceptFunctionalRelations(concept,selection,parts).some(r=>['arterial','venous'].includes(r.kind)));
 }
});
test('male registered rectum routes to the registered sigmoid without borrowing differently framed vessels',()=>{
 const atlas=JSON.parse(readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
 const selected=atlas.rectalWallIntegration.partIds.map(id=>atlas.parts.find(p=>p.id===id));
 const concept=atlas.concepts.find(c=>c.id==='FMA14544');assert.ok(concept);
 const links=anatomicalConceptFunctionalRelations(concept,selected,atlas.parts);
 assert.deepEqual(links.filter(r=>r.kind==='before').map(r=>r.target.id),['LOCAL:male:CAECOWALL:sigmoid:serosa']);
 assert.deepEqual(anatomicalConceptFunctionalRelations({id:'hierarchy:guest:ta98:A05.7.04.001',name:'rectum',elements:concept.elements},selected,atlas.parts).map(r=>[r.kind,r.target.id]),links.map(r=>[r.kind,r.target.id]));
 assert.ok(!links.some(r=>r.target.suppressed||r.target.id.startsWith('HRA:')||r.kind==='arterial'&&r.target.id.startsWith('ZA:')));
 const retired=atlas.parts.find(p=>p.id===atlas.rectalWallIntegration.sourceRectum);
 const foreign={...retired,id:'HRA:unregistered-rectum'};
 const selection=selected.map(p=>({...p,bowelWallSource:foreign.id}));
 const parts=atlas.parts.filter(p=>p.id!==retired.id&&!selected.some(s=>s.id===p.id)).concat(foreign,selection);
 assert.deepEqual(anatomicalConceptFunctionalRelations(concept,selection,parts),[]);
});
