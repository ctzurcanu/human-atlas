import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/anatomical-components';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const ids=['native-remaining-surface','left-repair','right-repair'];
function fixture(){
 const parts=new Map(ids.map((id,i)=>[id,{id,conceptId:i?'repair:'+id:'whole',name:id,system:'respiratory',bounds:[[0,0,0],[1,1,1]]}]));
 const whole={id:'whole',name:'Whole cartilage',elements:[...ids],componentElements:[...ids]};
 return {parts,whole,concepts:[whole]};
}
test('Show components retains the original remaining surface and both new patches',()=>{
 const f=fixture();assert.deepEqual(api.anatomicalComponentIds(f.whole,[],f.parts,f.concepts),ids);
 const contains=api.anatomicalComponentRelations(f.whole,[],f.parts,f.concepts).filter(r=>r.kind==='contains');
 assert.deepEqual(contains.map(r=>r.target.id),ids);
});
test('a malformed declaration cannot introduce an unrelated or absent piece',()=>{
 for(const declaration of [[...ids,'unrelated'],[ids[0],ids[0]],['absent',ids[1]],null]){
  const f=fixture();f.whole.componentElements=declaration;
  assert.ok(!api.anatomicalComponentIds(f.whole,[],f.parts,f.concepts).includes(ids[0]));
 }
 const f=fixture();f.parts.get(ids[0]).suppressed=true;
 assert.ok(!api.anatomicalComponentIds(f.whole,[],f.parts,f.concepts).includes(ids[0]));
});
test('a partial selection cannot acquire the missing original surface through declaration',()=>{
 const f=fixture(),partial={...f.whole,elements:ids.slice(1)};
 assert.deepEqual(api.anatomicalComponentIds(partial,[],f.parts,f.concepts),ids.slice(1));
});
test('legacy whole-source identity remains excluded without explicit component membership',()=>{
 const f=fixture();delete f.whole.componentElements;
 assert.deepEqual(api.anatomicalComponentIds(f.whole,[],f.parts,f.concepts),ids.slice(1));
});
