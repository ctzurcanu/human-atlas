import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations,anatomicalInnervationNames}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
const male=JSON.parse(readFileSync('.local-models/ta98-runtime.json','utf8')).parts;
const left=male.find(p=>p.id==='ZA:External anal sphincter.l');
const right=male.find(p=>p.id==='ZA:External anal sphincter.r');
assert.ok(left&&right);

test('actual male halves resolve ipsilateral pudendal supply through explicitly missing inferior anal nerves',()=>{
 for(const muscle of [left,right]){
 const which=muscle.name.includes('(left)')?'left':'right';
  assert.ok(anatomicalInnervationNames(muscle).includes('Pudendal nerve'));
  const links=anatomicalRelations(muscle,male).filter(r=>r.kind==='innervation');
  const route=links.find(r=>r.target.name===`Pudendal nerve (${which})`);
  assert.ok(route);assert.deepEqual(route.via,['Inferior anal nerves']);assert.equal(route.viaModeled,false);
  assert.ok(!links.some(r=>r.target.name===`Pudendal nerve (${which==='left'?'right':'left'})`));
  const reverse=anatomicalRelations(route.target,male).find(r=>r.kind==='innervates'&&r.target.id===muscle.id);
  assert.ok(reverse);assert.deepEqual(reverse.via,route.via);
  assert.match(route.note,/variants/);assert.match(route.note,/does not certify/);
 }
});

test('native female route support does not borrow male, suppressed or donor anatomy',()=>{
 const muscle={...left,id:'LOCAL:female:external-sphincter',conceptId:'test:female:sphincter'};
 const pudendal={...male.find(p=>p.name==='Pudendal nerve (left)'),id:'LOCAL:female:pudendal'};
 const artery={...pudendal,id:'LOCAL:female:wrong-system',system:'arterial'};
 const donor={...pudendal,id:'HRA:independent-pudendal'};
 const suppressed={...pudendal,id:'LOCAL:female:hidden-pudendal',suppressed:true};
 const wrongSex={...pudendal,id:'LOCAL:male:foreign-pudendal'};
 const parts=[muscle,pudendal,artery,donor,suppressed,wrongSex];
 const links=anatomicalRelations(muscle,parts).filter(r=>r.kind==='innervation');
 assert.deepEqual(links.map(r=>r.target.id),[pudendal.id]);
 const inferior={...pudendal,id:'LOCAL:female:inferior-anal',name:'Inferior anal nerves (left)'};
 const modeled=[...parts,inferior];
 const direct=anatomicalRelations(muscle,modeled).filter(r=>r.kind==='innervation');
 assert.deepEqual(direct.map(r=>r.target.id),[inferior.id]);assert.equal(direct[0].via,undefined);
 assert.ok(anatomicalRelations(inferior,modeled).some(r=>r.kind==='innervates'&&r.target.id===muscle.id));
 for(const invalid of [{...muscle,suppressed:true},{...muscle,name:'Internal anal sphincter'},{...muscle,name:'Subcutaneous part of external anal sphincter'},{...muscle,system:'digestive'}]){
  assert.ok(!anatomicalRelations(invalid,[invalid,pudendal]).some(r=>r.kind==='innervation'));
 }
});

test('the current female catalogue does not acquire invented sphincter geometry or nerve targets',()=>{
 const female=JSON.parse(readFileSync('.local-models/ta98-female-runtime.json','utf8')).parts;
 assert.equal(female.filter(p=>p.name.toLowerCase().startsWith('external anal sphincter')).length,0);
 for(const pudendal of female.filter(p=>p.name.startsWith('Pudendal nerve'))){
  assert.ok(!anatomicalRelations(pudendal,female).some(r=>r.kind==='innervates'&&r.target.name.toLowerCase().includes('external anal sphincter')));
 }
});
