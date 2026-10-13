import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations,anatomicalConceptFunctionalRelations}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
const parts=JSON.parse(readFileSync('.local-models/ta98-female-runtime.json')).parts;
const organs=parts.filter(p=>!p.suppressed&&p.id.startsWith('LOCAL:female:')&&['Uterus','Cervix','Vagina','Ovary (left)','Ovary (right)','Fallopian tube (left)','Fallopian tube (right)'].includes(p.name));
assert.equal(organs.length,7);
const functional=new Set(['arterial','venous','innervation','supplies','drains','innervates']);
const signature=rows=>rows.filter(r=>functional.has(r.kind)).map(r=>`${r.kind}:${r.target.id}`).sort();
test('all seven native organs retain their available functional links without borrowing matching foreign or retired targets',()=>{
 for(const organ of organs){
  const baseline=anatomicalRelations(organ,parts).filter(r=>functional.has(r.kind));
  assert.ok(baseline.length,organ.name);
  const invalid=baseline.flatMap(({target})=>[
   {...target,id:`HRA:foreign:${target.id}`},
   {...target,id:target.id.replace('LOCAL:female:','LOCAL:male:')},
   {...target,id:target.id+':retired',suppressed:true},
  ]);
  const snapshot=[...parts,...invalid];
  assert.deepEqual(signature(anatomicalRelations(organ,snapshot)),signature(baseline),organ.name);
  for(const target of invalid)assert.ok(!anatomicalRelations(target,snapshot).some(r=>functional.has(r.kind)&&r.target.id===organ.id),`No reverse link from ${target.id}`);
 }
});
test('foreign reproductive organs cannot acquire native vessels through reverse supply lookup',()=>{
 const uterus=organs.find(p=>p.name==='Uterus');
 const supply=anatomicalRelations(uterus,parts).filter(r=>['arterial','venous'].includes(r.kind));
 const invalid=[{...uterus,id:'HRA:foreign-uterus'},{...uterus,id:'LOCAL:male:foreign-uterus'},{...uterus,id:uterus.id+':retired',suppressed:true}];
 const snapshot=[...parts,...invalid];
 for(const organ of invalid){
  assert.deepEqual(signature(anatomicalRelations(organ,snapshot)),[]);
  for(const {target} of supply)assert.ok(!anatomicalRelations(target,snapshot).some(r=>functional.has(r.kind)&&r.target.id===organ.id));
 }
});
test('whole-organ functional selection keeps the native frame and rejects mixed or retired assemblies',()=>{
 const uterus=organs.find(p=>p.name==='Uterus');
 const leaves=['Body of uterus','Fundus of uterus'].map((name,i)=>({...uterus,id:`LOCAL:female:test-region:${i}`,name}));
 const concept={id:'test:whole-uterus',name:'Uterus',elements:leaves.map(p=>p.id)};
 assert.ok(anatomicalConceptFunctionalRelations(concept,leaves,[...parts,...leaves]).some(r=>r.kind==='arterial'));
 for(const replacement of [{...leaves[1],id:'HRA:foreign-fundus'},{...leaves[1],id:'LOCAL:male:foreign-fundus'},{...leaves[1],suppressed:true}]){
  const selection=[leaves[0],replacement],mixed={...concept,elements:selection.map(p=>p.id)};
  assert.deepEqual(anatomicalConceptFunctionalRelations(mixed,selection,[...parts,...selection]),[]);
 }
});
