import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/spleen-wall';export * from './app/viewer-state';export * from './app/bowel-wall-peeling';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const base={selected:[],hidden:[],visible:['lymphatic'],depthHidden:[],view:'front',region:'all',contextOpacity:0,skinOpacity:0,explode:0,isolate:false,rotate:false,reset:0,focus:0};
for(const sex of ['male','female'])test(sex+': installed tissue controls peel only the four canonical cells',()=>{
 const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),m=a.finiteSpleenIntegration,choice=a.concepts.find(c=>c.id===m.wholeConceptId),parts=new Map(a.parts.map(p=>[p.id,p])),steps=api.bowelWallPeelSteps(a,choice);
 assert.deepEqual(steps.map(s=>s.ids[0]),['serosa','capsule','pulp','trabeculae'].map(t=>'TA98PART:'+sex+':Spleen:'+t));
 const state=api.readViewUrl('?select='+encodeURIComponent(choice.id),a,base);assert.deepEqual(new Set(state.selected),new Set(m.partIds));assert.ok(state.hidden.includes(m.sourcePartId));assert.ok(m.partIds.every(id=>!state.hidden.includes(id)));
 for(let i=0;i<=steps.length;i++){const hidden=steps.slice(0,i).flatMap(s=>s.ids),s={...state,hidden};assert.equal(api.wallPeelDepth(steps,hidden),i);assert.deepEqual(m.partIds.filter(id=>api.partVisible(parts.get(id),s,a.parts)),steps.slice(i).flatMap(s=>s.ids));}
 const explicit=api.readViewUrl('?hide='+encodeURIComponent(m.sourcePartId),a,base);assert.ok(m.partIds.every(id=>explicit.hidden.includes(id)),'An explicit user hide still follows the replaced organ');
 const bad=structuredClone(a);bad.parts.find(p=>p.id===m.partIds[0]).geometryOrigin=[0,0,0];assert.equal(api.spleenTissuePeelSteps(bad,choice),null);assert.equal(api.spleenTissuePeelSteps(a,{...choice,id:'selection-set'}),null);
});
test('retired catalogue defaults never hide their installed replacements',()=>{
 const bounds=[[0,1,0],[.1,1.1,.1]],p=id=>({id,conceptId:id,name:id,system:'lymphatic',bounds});const a={parts:[{...p('retired'),suppressed:true,replacementParts:['newA','newB']},p('newA'),p('newB')],concepts:[{id:'whole',name:'Whole',elements:['newA','newB']}],defaultHidden:['retired']};
 assert.deepEqual(api.readViewUrl('',a,base).hidden,['retired']);assert.deepEqual(api.readViewUrl('?hide=retired',a,base).hidden,['newA','newB']);
});
