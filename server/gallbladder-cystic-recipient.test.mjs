import assert from 'node:assert/strict';
import test from 'node:test';
import {existsSync} from 'node:fs';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/gallbladder-relations.ts','app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/gallbladder-cystic-recipient-tests'});
const load=async name=>import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles.find(f=>f.path.endsWith(name+'.js')).contents).toString('base64'));
const {gallbladderRelations}=await load('gallbladder-relations');
const {anatomicalRelations,anatomicalConceptFunctionalRelations}=await load('anatomical-relations');
const paths=['.local-models/organ-models/spinal-cord/male-body.json','.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json','.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'];
if(!paths.some(existsSync)){
 test('Cystic recipient regression requires the private installed anatomy catalogues',{skip:'Local anatomy sources are intentionally absent'},()=>{});
}else{
// A partially present or corrupt installation must fail; it is not a skip.
const raw=await Promise.all(paths.map(path=>readFile(path)));
const [male,master,female]=raw.map(bytes=>JSON.parse(bytes));
const arteryId='TA98CYSTICBLOOD:male:artery',parentId='TA98HEPATIC:male:artery';
const artery=male.parts.find(p=>p.id===arteryId),parent=male.parts.find(p=>p.id===parentId);
const veinId='TA98CYSTICBLOOD:male:vein',portalId='TA98HEPATIC:male:portal';
const vein=male.parts.find(p=>p.id===veinId),portal=male.parts.find(p=>p.id===portalId);
const whole=male.concepts.find(c=>c.id==='ZA:Gallbladder');
const selected=whole.elements.map(id=>male.parts.find(p=>p.id===id));
const covering=selected.filter(p=>p.tissueWallLayer==='covering');
const withoutArtery=male.parts.filter(p=>p.id!==arteryId);
const replace=(original,bad)=>male.parts.map(p=>p===original?bad:p);
const links=(graph,id,kind)=>graph.get(id)?.filter(row=>row.kind===kind&&row.target.id===arteryId)??[];

test('Installed ordinary-body gallbladder and the cystic source use the reviewed live assembly',async()=>{
 assert.ok(artery&&parent);assert.equal(artery.sourceId,'HRA:VH_F_cystic_artery');
 assert.equal(male.parts.find(p=>p.id==='ZA:Gallbladder').suppressed,true);
 assert.equal(selected.length,16);assert.ok(selected.every(p=>!p.suppressed&&p.tissueWallAssembly==='gallbladder:male'&&p.sourceId==='Gallbladder'));
 const masterArtery=master.parts.find(p=>p.id===arteryId);
 assert.equal(covering.length,4);assert.deepEqual({...artery,chunk:undefined},{...masterArtery,chunk:undefined});
 assert.deepEqual(male.chunks[artery.chunk],master.chunks[masterArtery.chunk]);
 const path=male.chunks[artery.chunk].url.replace(/^\/local-models\//,'.local-models/');
 const bank=await readFile(path);assert.equal(createHash('sha256').update(bank).digest('hex'),'8fa533756b17f71c194ec5a9e51ea4e24229a01ab5342fa93106575091bf4947');
});

test('Actual active covering cells and the whole concept receive reciprocal cystic arterial supply',()=>{
 const graph=gallbladderRelations(male.parts),back=graph.get(arteryId);
 for(const cell of covering){
  const row=links(graph,cell.id,'arterial')[0];assert.ok(row,cell.id);
  assert.ok(back.some(r=>r.kind==='supplies'&&r.target===cell&&r.note===row.note));
  assert.ok(anatomicalRelations(cell,male.parts).some(r=>r.kind==='arterial'&&r.target.id===arteryId));
  assert.match(row.note,/entry into the gallbladder fundus/);assert.match(row.note,/Microscopic perfusion.*thickness.*provisional/);assert.match(row.note,/upstream.*unverified/);
  if(cell.id.endsWith(':A05.8.02.002'))assert.doesNotMatch(row.note,/no local entry/);
  else assert.match(row.note,/Whole-gallbladder arterial reference only; no local entry/);
 }
 for(const inner of selected.filter(p=>p.tissueWallLayer!=='covering'))assert.equal(links(graph,inner.id,'arterial').length,0,'No direct supply is invented for an inner coat');
 const rows=anatomicalConceptFunctionalRelations(whole,selected,male.parts);
 assert.ok(rows.some(r=>r.kind==='arterial'&&r.target.id===arteryId));
 assert.ok(anatomicalRelations(artery,male.parts).some(r=>r.kind==='supplies'&&r.targetConcept?.id===whole.id),'Reciprocal whole-organ relation selects the active registered assembly');
});

test('Unrelated, shifted, opposite-frame and incomplete cystic donors cannot borrow the installed relation',()=>{
 const invalid=[{...artery,id:'HRA:VH_F_cystic_artery'},{...artery,conceptId:'unrelated'},{...artery,sourceId:'unrelated'},
  {...artery,hepaticSourceModel:'local-female-ta98'},{...artery,sourceOffset:[0,0,.001]},{...artery,geometryOrigin:[0,0,0]},
  {...artery,ta98Term:'A12.2.12.015'},{...artery,pairedTissue:false},{...artery,sectionSolid:false},{...artery,normalSigns:undefined},
  {...artery,suppressed:true},{...artery,tissueDrawGroups:artery.tissueDrawGroups.filter(g=>!g.neighbor)}];
 for(const bad of invalid){
  const parts=replace(artery,bad),graph=gallbladderRelations(parts);
  for(const cell of covering)assert.ok(!graph.get(cell.id)?.some(r=>r.kind==='arterial'&&r.target.id===bad.id),JSON.stringify(bad));
  assert.ok(!anatomicalConceptFunctionalRelations(whole,selected,parts).some(r=>r.kind==='arterial'&&r.target.id===bad.id));
  assert.ok(!anatomicalRelations(covering[0],parts).some(r=>r.kind==='arterial'&&r.target.id===bad.id));
  assert.ok(!anatomicalRelations(bad,parts).some(r=>r.kind==='supplies'&&covering.some(p=>p.id===r.target.id)));
 }
 for(const bad of [{...parent,suppressed:true},{...parent,hepaticSourceModel:'local-female-ta98'},
  {...parent,tissueDrawGroups:parent.tissueDrawGroups.filter(g=>g.neighbor!==arteryId)}]){
  const graph=gallbladderRelations(replace(parent,bad));
  for(const cell of covering)assert.equal(links(graph,cell.id,'arterial').length,0);
 }
});

test('Matching gallbladder donor names and incomplete recipient assemblies receive no covering links',()=>{
 const source=male.parts.find(p=>p.id==='ZA:Gallbladder');
 const donor={...source,id:'HRA:VH_F_gallbladder',conceptId:'HRA:VH_F_gallbladder',suppressed:false};
 assert.equal(gallbladderRelations([donor,artery,parent]).get(donor.id),undefined);
 for(const parts of [male.parts.filter(p=>p!==source),male.parts.filter(p=>p!==selected[0]),
  replace(covering[0],{...covering[0],sourceOffset:[0,0,.001]}),
  replace(covering[1],{...covering[1],geometryScale:.9}),
  replace(selected.find(p=>p.tissueWallLayer==='mucosa'),{...selected.find(p=>p.tissueWallLayer==='mucosa'),tissueWallLayer:'covering'})]){
  const graph=gallbladderRelations(parts);for(const cell of covering)assert.equal(links(graph,cell.id,'arterial').length,0);
 }
});

test('Installed cystic vein and actual covering cells have reciprocal source-specific drainage',async()=>{
 assert.ok(vein&&portal,'The reviewed vein must be present in the ordinary body');
 assert.equal(vein.sourceId,'HRA:VH_F_cystic_vein');assert.equal(vein.ta98Term,'A12.3.12.012');
 const masterVein=master.parts.find(p=>p.id===veinId);
 assert.deepEqual({...vein,chunk:undefined},{...masterVein,chunk:undefined});
 assert.deepEqual(male.chunks[vein.chunk],master.chunks[masterVein.chunk]);
 const bank=await readFile(male.chunks[vein.chunk].url.replace(/^\/local-models\//,'.local-models/'));
 assert.equal(createHash('sha256').update(bank).digest('hex'),'fb7943d2f5c20f4e7848a386649ef1b646ad725da9eadc82d54a1173ad3392ef');
 const graph=gallbladderRelations(male.parts);
 for(const cell of covering){
  const row=graph.get(cell.id)?.find(r=>r.kind==='venous'&&r.target===vein);assert.ok(row,cell.id);
  assert.ok(graph.get(veinId)?.some(r=>r.kind==='drains'&&r.target===cell&&r.note===row.note));
  assert.ok(anatomicalRelations(cell,male.parts).some(r=>r.kind==='venous'&&r.target.id===veinId));
  assert.match(row.note,/right-portal venous junction/);assert.match(row.note,/without breaching its inner surface/);
  assert.match(row.note,/one source-specific venous route/);
  if(cell.id.endsWith(':A05.8.02.002'))assert.doesNotMatch(row.note,/no local entry/);
  else assert.match(row.note,/Whole-gallbladder venous reference only; no local entry/);
 }
 for(const inner of selected.filter(p=>p.tissueWallLayer!=='covering'))assert.ok(!graph.get(inner.id)?.some(r=>r.kind==='venous'&&r.target===vein));
 assert.ok(anatomicalConceptFunctionalRelations(whole,selected,male.parts).some(r=>r.kind==='venous'&&r.target.id===veinId));
 assert.ok(anatomicalRelations(vein,male.parts).some(r=>r.kind==='drains'&&r.targetConcept?.id===whole.id));
});

test('Invalid cystic veins and portal interfaces cannot borrow the installed drainage',()=>{
 assert.ok(vein&&portal);
 for(const bad of [{...vein,id:'HRA:VH_F_cystic_vein'},{...vein,sourceId:'unrelated'},
  {...vein,ta98Term:'A12.3.12.001'},{...vein,system:'arterial'},
  {...vein,hepaticSourceModel:'local-female-ta98'},{...vein,sourceOffset:[0,0,.001]},
  {...vein,geometryScale:.9},{...vein,suppressed:true},{...vein,normalSigns:undefined},
  {...vein,tissueWallAssembly:'hepatic-connections:male:artery'},
  {...vein,tissueDrawGroups:vein.tissueDrawGroups.filter(g=>g.neighbor!==portalId)}]){
  const parts=replace(vein,bad),graph=gallbladderRelations(parts);
  for(const cell of covering)assert.ok(!graph.get(cell.id)?.some(r=>r.kind==='venous'&&r.target.id===bad.id));
  assert.ok(!anatomicalConceptFunctionalRelations(whole,selected,parts).some(r=>r.kind==='venous'&&r.target.id===bad.id));
  assert.ok(!anatomicalRelations(bad,parts).some(r=>r.kind==='drains'&&covering.some(p=>r.target.id===p.id)));
 }
 for(const bad of [{...portal,suppressed:true},{...portal,hepaticSourceModel:'local-female-ta98'},
  {...portal,tissueDrawGroups:portal.tissueDrawGroups.filter(g=>g.neighbor!==veinId)}]){
  const graph=gallbladderRelations(replace(portal,bad));
  for(const cell of covering)assert.ok(!graph.get(cell.id)?.some(r=>r.kind==='venous'&&r.target.id===veinId));
 }
 const concept=female.concepts.find(c=>c.id==='LOCAL:female:gallbladder');
 const members=concept.elements.map(id=>female.parts.find(p=>p.id===id));
 const before=anatomicalConceptFunctionalRelations(concept,members,female.parts);
 const after=anatomicalConceptFunctionalRelations(concept,members,[...female.parts,vein,portal]);
 assert.deepEqual(after.map(r=>[r.kind,r.target.id,r.note,r.via]),before.map(r=>[r.kind,r.target.id,r.note,r.via]));
});

test('Native female fallback and existing nonarterial relations stay unchanged',()=>{
 const concept=female.concepts.find(c=>c.id==='LOCAL:female:gallbladder');assert.ok(concept);
 const members=concept.elements.map(id=>female.parts.find(p=>p.id===id));
 const baseline=anatomicalConceptFunctionalRelations(concept,members,female.parts);
 assert.ok(baseline.some(r=>r.kind==='arterial'&&r.target.id.startsWith('LOCAL:female:')),'The actual female native arterial fallback remains available');
 const mixed=anatomicalConceptFunctionalRelations(concept,members,[...female.parts,artery,parent]);
 assert.deepEqual(mixed.map(r=>[r.kind,r.target.id,r.note,r.via]),baseline.map(r=>[r.kind,r.target.id,r.note,r.via]));
 const nonarterial=parts=>anatomicalConceptFunctionalRelations(whole,selected,parts).filter(r=>r.kind!=='arterial').map(r=>[r.kind,r.target.id,r.note,r.via]);
 assert.deepEqual(nonarterial(male.parts),nonarterial(withoutArtery));
});

test('The regression does not rewrite the installed geometry masters or ordinary body',async()=>{
 for(let i=0;i<paths.length;i++)assert.deepEqual(await readFile(paths[i]),raw[i]);
});
}
