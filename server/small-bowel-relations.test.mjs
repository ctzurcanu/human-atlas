import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile,writeFile} from 'node:fs/promises';
import {build} from 'esbuild';

const built=await build({entryPoints:['app/anatomical-relations.ts','app/small-bowel-relations.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/small-bowel-relations-tests'});
const load=async name=>import('data:text/javascript;base64,'+Buffer.from(built.outputFiles.find(f=>f.path.endsWith(name+'.js')).contents).toString('base64'));
const {anatomicalRelations}=await load('anatomical-relations');
const {smallBowelRelations,smallBowelSourceLimitation}=await load('small-bowel-relations');
const [male,female]=await Promise.all(['.local-models/ta98-runtime.json','.local-models/ta98-female-runtime.json'].map(async p=>JSON.parse(await readFile(p))));
const functional=new Set(['arterial','venous','innervation','supplies','drains','innervates']);
function check(parts,id,target,kind,reverse,via){
 const a=parts.find(p=>p.id===id),b=parts.find(p=>p.id===target);assert.ok(a&&b);
 const forward=anatomicalRelations(a,parts).find(r=>r.kind===kind&&r.target.id===target);
 const back=anatomicalRelations(b,parts).find(r=>r.kind===reverse&&r.target.id===id);
 assert.ok(forward&&back);assert.deepEqual(forward.via,via);assert.deepEqual(back.via,via);assert.equal(back.note,forward.note);assert.match(forward.note,/unverified/);
 return forward;
}
test('Both female bowel regions resolve actual native vagus aliases through missing autonomic intermediates',()=>{
 for(const id of ['LOCAL:female:small_intestine_jejunum:182','LOCAL:female:small_intestine_ileum:181'])for(const nerve of ['LOCAL:female:l_vagus_nerve:894','LOCAL:female:r_vagus_nerve:895'])check(female.parts,id,nerve,'innervation','innervates',['Superior mesenteric plexus','Enteric plexuses']);
});
test('Named native arterial and venous collectors retain branch intermediates in both directions',()=>{
 for(const [id,title] of [['LOCAL:female:small_intestine_jejunum:182','Jejunal'],['LOCAL:female:small_intestine_ileum:181','Ileal']]){
  check(female.parts,id,'LOCAL:female:superior_mesenteric_artery:1670','arterial','supplies',[title+' arteries']);
  check(female.parts,id,'LOCAL:female:superior_mesenteric_vein:1755','venous','drains',[title+' veins']);
 }
 check(male.parts,'ZA:Jejunum','ZA:Superior mesenteric artery','arterial','supplies',['Jejunal arteries']);
 check(male.parts,'ZA:Jejunum','ZA:Superior mesenteric vein','venous','drains',['Jejunal veins']);
 for(const nerve of ['ZA:Vagus nerve (X).l','ZA:Vagus nerve (X).r'])check(male.parts,'ZA:Jejunum',nerve,'innervation','innervates',['Superior mesenteric plexus','Enteric plexuses']);
});
test('Generic inverse links cannot restore donor or opposite-frame bowel supply',()=>{
 for(const model of [male,female]){
  const artery=model.parts.find(p=>p.name==='Superior mesenteric artery');assert.ok(artery);
  const donor={...female.parts.find(p=>p.id==='LOCAL:female:small_intestine_ileum:181'),id:'HRA:test_ileum',conceptId:'HRA:test_ileum',sourceId:'test_ileum',name:'Ileum'};
  const parts=[artery,donor];
  assert.ok(!anatomicalRelations(artery,parts).some(r=>functional.has(r.kind)&&r.target.id===donor.id));
  assert.ok(!anatomicalRelations(donor,parts).some(r=>functional.has(r.kind)&&r.target.id===artery.id));
 }
 const ileum=male.parts.find(p=>p.id==='LOCAL:male:small_intestine_illium:1263');
 assert.ok(ileum);assert.ok(!anatomicalRelations(ileum,male.parts).some(r=>functional.has(r.kind)&&r.target.id.startsWith('ZA:')));
});
test('Native bowel continuations exclude independent references without skipping an available unregistered stage',()=>{
 const j=female.parts.find(p=>p.id==='LOCAL:female:small_intestine_jejunum:182');
 const links=anatomicalRelations(j,female.parts);
 assert.deepEqual(links.filter(r=>r.kind==='before').map(r=>r.target.id),['LOCAL:female:intestine_duodenum:183']);
 assert.deepEqual(links.filter(r=>r.kind==='after').map(r=>r.target.id),['LOCAL:female:small_intestine_ileum:181']);
 const mj=male.parts.find(p=>p.id==='ZA:Jejunum');
 assert.ok(!anatomicalRelations(mj,male.parts).some(r=>r.kind==='after'),'A differently registered ileal source does not authorize a fabricated joined course');
});
test('Suppressed, wrong-system or renamed source organs do not acquire reviewed functional links',()=>{
 const source=female.parts.find(p=>p.id==='LOCAL:female:small_intestine_jejunum:182'),artery=female.parts.find(p=>p.id==='LOCAL:female:superior_mesenteric_artery:1670');
 for(const bad of [{...source,suppressed:true},{...source,system:'muscular'},{...source,name:'Ileum'},{...source,sourceId:'other'}])assert.equal(smallBowelRelations([bad,artery]).size,0);
 assert.equal(smallBowelRelations([source,{...artery,system:'venous'}]).size,0);
 assert.match(smallBowelSourceLimitation(source),/Finite wall.*mesenteric support.*unverified/);
 assert.equal(smallBowelSourceLimitation({...source,id:'HRA:donor'}),undefined);
 const graph=smallBowelRelations(female.parts);assert.equal(graph,smallBowelRelations(female.parts));
});
test('Current route evidence records frames and source limitations without geometry resources',async()=>{
 const records=[male,female].map(model=>({links:[...smallBowelRelations(model.parts)].map(([part,links])=>({part,links:links.map(r=>({kind:r.kind,target:r.target.id,via:r.via,note:r.note}))}))}));
 await writeFile('reports/ta98-completeness/small-bowel-native-functional-relations.json',JSON.stringify({records,addedTriangles:0,addedTextures:0,addedDrawGroups:0,geometryModified:false,accepted:false},null,2)+'\n');
});
