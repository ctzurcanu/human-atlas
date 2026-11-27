import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/suprarenal-wall';export * from './app/viewer-state';export * from './app/bowel-wall-peeling';export * from './app/guest-hierarchy';export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const base={selected:[],hidden:[],visible:['endocrine'],depthHidden:[],view:'front',region:'all',contextOpacity:0,skinOpacity:0,explode:0,isolate:false,rotate:false,reset:0,focus:0};
function fixture(){
 const a={sex:'female',parts:[],concepts:[],defaultHidden:[],finiteSuprarenalIntegration:{schema:'suprarenal-local-tissues-v223',sex:'female',glands:[],partIds:[]}};
 for(const side of ['left','right']){
  const source='source:'+side,whole='gland:'+side,roles=['capsule','cortex','medulla'],ids=roles.map(role=>`TA98PART:female:Suprarenal:${side}:${role}`);
  a.finiteSuprarenalIntegration.glands.push({side,sourcePartId:source,wholeConceptId:whole,partIds:ids});a.finiteSuprarenalIntegration.partIds.push(...ids);a.parts.push({id:source,name:side,conceptId:whole,system:'endocrine',suppressed:true,replacementParts:ids,bounds:[[0,1,0],[.1,1.1,.1]]});a.defaultHidden.push(source);a.concepts.push({id:whole,name:side,elements:ids});
  for(const [i,role]of roles.entries())a.parts.push({id:ids[i],name:role+' '+side,conceptId:ids[i],system:'endocrine',pairedTissue:true,sectionSolid:true,tissueWallAssembly:'suprarenal:female:'+side,tissueWallLayer:role,normalSigns:200+i,geometryOrigin:[0,1.25,0],positions:0,normals:100,chunk:1,vertexCount:20,indexCount:role==='medulla'?6:12,tissueAppearance:{color:[.5,.5,.5],pattern:'solid'},bounds:[[0,1,0],[.1,1.1,.1]],tissueDrawGroups:role==='capsule'?[{start:0,count:6},{start:6,count:6,neighbor:ids[1]}]:role==='cortex'?[{start:0,count:6,neighbor:ids[0]},{start:6,count:6,neighbor:ids[2]}]:[{start:0,count:6,neighbor:ids[1]}]});
 }
 return a;
}
test('Complete single and paired gland concepts peel capsule, cortex, medulla',()=>{
 const a=fixture(),g=a.concepts[0];assert.deepEqual(api.suprarenalTissuePeelSteps(a,g).map(x=>x.ids),g.elements.map(id=>[id]));
 const choice={id:'all-glands',elements:a.finiteSuprarenalIntegration.partIds};const steps=api.suprarenalTissuePeelSteps(a,choice);assert.equal(steps.length,3);assert.ok(steps.every(s=>s.ids.length===2));assert.equal(api.suprarenalTissuePeelSteps(a,{...choice,id:'selection-set'}),null);
 assert.equal(api.suprarenalTissuePeelSteps(a,{...g,elements:g.elements.slice(1)}),null);assert.equal(api.suprarenalTissuePeelSteps(a,{...g,elements:[...g.elements,'unrelated']}),null);
});
test('Retired gland defaults do not hide new tissue, while explicit hides do',()=>{
 const a=fixture(),g=a.concepts[0],state=api.readViewUrl('?select='+encodeURIComponent(g.id),a,base),steps=api.bowelWallPeelSteps(a,g),byId=new Map(a.parts.map(p=>[p.id,p]));assert.deepEqual(state.selected,g.elements);
 for(let i=0;i<=3;i++){const hidden=steps.slice(0,i).flatMap(x=>x.ids);assert.equal(api.wallPeelDepth(steps,hidden),i);assert.deepEqual(g.elements.filter(id=>api.partVisible(byId.get(id),{...state,hidden},a.parts)),steps.slice(i).flatMap(x=>x.ids));}
 assert.ok(g.elements.every(id=>!state.hidden.includes(id)));assert.ok(g.elements.every(id=>api.readViewUrl('?hide=source%3Aleft',a,base).hidden.includes(id)));
});
for(const [name,mutate]of [
 ['different source frame',a=>a.parts.find(p=>p.tissueWallLayer==='capsule').geometryOrigin=[0,0,0]],
 ['broken reciprocal tissue interface',a=>a.parts.find(p=>p.tissueWallLayer==='cortex').tissueDrawGroups[0].count=3],
 ['wrong sex marker',a=>a.finiteSuprarenalIntegration.sex='male'],
 ['retired source still active',a=>a.parts[0].suppressed=false],
 ['missing finite tissue',a=>a.parts=a.parts.filter(p=>p.id!==a.concepts[0].elements[2])],
])test('Refuse '+name,()=>{const a=fixture();mutate(a);assert.equal(api.suprarenalTissuePeelSteps(a,a.concepts[0]),null);});
test('Installed glands expose only their actual six finite tissues',()=>{
 for(const sex of ['male','female']){const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));assert.ok(a.finiteSuprarenalIntegration,sex+' suprarenal installation');
  for(const g of a.finiteSuprarenalIntegration.glands){const c=a.concepts.find(c=>c.id===g.wholeConceptId),steps=api.bowelWallPeelSteps(a,c);assert.equal(steps.length,3);assert.deepEqual(steps.flatMap(s=>s.ids),g.partIds);}
 }
});
test('TA98 contains actual cortex/medulla controls and additional capsule controls',()=>{
 for(const sex of ['female','male']){
 const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 const hierarchy=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))),roots=api.resolveGuestHierarchy(a,hierarchy),nodes=new Map(),todo=[...roots];
 while(todo.length){const n=todo.pop();if(nodes.has(n.id))continue;nodes.set(n.id,n);todo.push(...n.children);}
 const whole=nodes.get('A11.5.00.001'),marker=a.finiteSuprarenalIntegration;
 assert.deepEqual(new Set(whole.parts.map(p=>p.id)),new Set(marker.partIds));
 assert.equal(whole.children.find(n=>n.id===`SUPRARENAL:${sex}:capsule-assembly`).parts.length,2);
 for(const [term,role]of [['A11.5.00.007','cortex'],['A11.5.00.008','medulla']])assert.deepEqual(new Set(nodes.get(term).directParts.map(p=>p.id)),new Set(marker.partIds.filter(id=>id.endsWith(':'+role))));
 assert.deepEqual(api.bowelWallPeelSteps(a,{id:'hierarchy:guest:ta98:A11.5.00.001',elements:whole.parts.map(p=>p.id)}).map(s=>s.ids.length),[2,2,2]);
 }
});
test('Complete installed glands retain the old native vessel and renal associations',()=>{
 for(const sex of ['female','male']){
  const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),m=a.finiteSuprarenalIntegration;assert.ok(m);
  const retired=new Set(m.glands.map(g=>g.sourcePartId)),finite=new Set(m.partIds);
  const before=a.parts.filter(p=>!finite.has(p.id)).map(p=>retired.has(p.id)?{...p,suppressed:false}:p);
  const keys=links=>links.filter(r=>['arterial','venous','innervation','adjacent'].includes(r.kind)).map(r=>r.kind+':'+r.target.id).sort();
  for(const g of m.glands){
   const whole=a.concepts.find(c=>c.id===g.wholeConceptId),selected=g.partIds.map(id=>a.parts.find(p=>p.id===id)),source=before.find(p=>p.id===g.sourcePartId);
   const expected=keys(api.anatomicalRelations(source,before));assert.ok(expected.some(k=>k.startsWith('adjacent:')));
   assert.deepEqual(keys(api.anatomicalConceptFunctionalRelations(whole,selected,a.parts)),expected);
   const renal=api.anatomicalConceptFunctionalRelations(whole,selected,a.parts).find(r=>r.kind==='adjacent');
   const reverse=api.anatomicalRelations(renal.target,a.parts).find(r=>r.kind==='adjacent'&&r.targetConcept?.id===whole.id);
   assert.deepEqual(new Set(reverse?.targetConcept?.elements),new Set(g.partIds));
   assert.equal(api.suprarenalWholeSource(selected.slice(1),a.parts),undefined);
   const mixed=structuredClone(selected);mixed[0].adrenalSourcePartId=m.glands.find(q=>q.side!==g.side).sourcePartId;
   assert.equal(api.suprarenalWholeSource(mixed,a.parts),undefined);
   const donor=structuredClone(selected);donor[1].sourceOffset=[0,0,0];assert.equal(api.suprarenalWholeSource(donor,a.parts),undefined);
  }
 }
});
