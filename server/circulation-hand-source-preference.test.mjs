import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/preferred-source-representations';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
const read=path=>JSON.parse(readFileSync(base+path));

test('review source choice hides mismatched nerves while retaining every original mesh and TA98 concept',()=>{
 const before=read('atlas-cardiovascular-review-v34/atlas.json'),after=read('atlas-cardiovascular-review-v36/atlas.json');
 for(const key of Object.keys(before).filter(k=>k!=='defaultHidden'))assert.deepEqual(after[key],before[key],key);
 const old=['O3M:Median nerve Proper palmar digital nerv.l','O3M:Median nerve Proper palmar digital nerv.r'];
 assert.deepEqual(after.defaultHidden.filter(id=>!before.defaultHidden.includes(id)),old);
 assert.deepEqual(after.handSourceSelection.hiddenReviewAlternatives,old);
 const displayed=api.withHiddenReviewAlternatives(after),byId=new Map(displayed.parts.map(p=>[p.id,p]));
 const state={selected:[],visible:[...new Set(after.parts.map(p=>p.system))],view:'front',skinOpacity:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:false,explode:0,isolate:false};
 for(const id of old){
  assert.equal(api.partVisible(byId.get(id),state,displayed.parts),false,'Show All cannot restore the second body source');
  assert.ok(api.withHiddenReviewAlternatives(after,true).parts.some(p=>p.id===id&&!p.suppressed),'repair view retains the original buffer');
 }
 for(const id of after.handSourceSelection.preferred){
  assert.ok(!byId.get(id).suppressed&&!after.defaultHidden.includes(id));
 }
 const camera=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5];
 const q=new URLSearchParams({frame:camera.join(','),select:after.handSourceSelection.preferred[0]});
 assert.deepEqual(api.readViewUrl(q.toString(),displayed,state).camera,camera);
});

test('bilateral adductor source choice retains native heads, every buffer and recoverable alternatives',()=>{
 const before=read('atlas-cardiovascular-review-v36/atlas.json'),after=read('atlas-cardiovascular-review-v37/atlas.json');
 for(const key of Object.keys(before).filter(k=>k!=='defaultHidden'))assert.deepEqual(after[key],before[key],key);
 const old=['O3M:Adductor pollicis.l','O3M:Adductor pollicis.r'];
 assert.deepEqual(after.defaultHidden.filter(id=>!before.defaultHidden.includes(id)),old);
 const marker=after.handAdductorSourceSelection;
 assert.deepEqual(marker.hiddenReviewAlternatives,old);assert.equal(marker.preferred.length,8);
 const displayed=api.withHiddenReviewAlternatives(after),byId=new Map(displayed.parts.map(p=>[p.id,p]));
 const state={selected:[],visible:[...new Set(after.parts.map(p=>p.system))],view:'front',skinOpacity:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:false,explode:0,isolate:false};
 for(const id of old){
  assert.equal(api.partVisible(byId.get(id),state,displayed.parts),false);
  assert.ok(api.withHiddenReviewAlternatives(after,true).parts.some(p=>p.id===id&&!p.suppressed));
 }
 for(const id of marker.preferred)assert.ok(byId.has(id)&&!byId.get(id).suppressed&&!after.defaultHidden.includes(id));
});
