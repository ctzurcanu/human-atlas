import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const work='.local-models/geometry-workbench/biodigital-male-circulation-v1';
const fixture=work+'/user-provided-blood-view-query.txt';
const available=fs.existsSync(fixture)&&fs.existsSync(work+'/atlas-review-v9/atlas.json');

test('the existing blood-review snapshot displays the source heart without changing camera or anatomical visibility intent',
 {skip:!available},async()=>{
 const compiled=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/viewer-state';export * from './app/view-part-codes';export * from './app/biodigital-heart';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
 const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
 const read=p=>JSON.parse(fs.readFileSync(p));
 const search=fs.readFileSync(fixture,'utf8');
 const query=new URLSearchParams(search);query.set('blood','biodigital-review');
 const path=api.localTa98LayerCatalogue('local-ta98',query,true);
 const current=read('.local-models/'+path.replace('/local-models/',''));
 const previous=read(work+'/atlas-review-v9/atlas.json');
 assert.ok(api.sourceCardiacIntegration(current),'The normal blood-review route must include the guarded BioDigital heart');
 assert.ok(api.compatiblePartAtlas(current,query.get('p')),'All old packed mesh indices retain their identities');
 const base={selected:[],visible:[],view:'three-quarter'};
 const before=api.readViewUrl(search,previous,base);
 const after=api.readViewUrl(search,current,base);
 assert.deepEqual(after.camera,before.camera);
 assert.deepEqual(after.camera,query.get('frame').split(',').map(Number));
 assert.deepEqual(after.visible,before.visible);
 assert.deepEqual(after.selected,before.selected);
 assert.equal(after.isolate,before.isolate);
 assert.deepEqual(after.section,before.section);
 assert.equal(after.skinOpacity,before.skinOpacity);
 assert.equal(after.contextOpacity,before.contextOpacity);
 const sources=api.sourceCardiacIntegration(current).sourceParts.filter(p=>p.role==='material');
 assert.ok(sources.length>=22);
 const parts=new Map(current.parts.map(p=>[p.id,p]));
 for(const source of sources){
  assert.equal(api.partVisible(parts.get(source.partId),after,current.parts),true,source.partId);
 }
 const restored=api.readViewUrl(api.viewUrl('http://127.0.0.1:3016/?blood=biodigital-review','local-ta98',after,after.camera,current),current,base);
 assert.deepEqual(restored.camera,after.camera);
 assert.deepEqual(restored.selected,after.selected);
 for(const source of sources)assert.equal(api.partVisible(parts.get(source.partId),restored,current.parts),true,source.partId+' survives saved-view round trip');
 const oldValve='BDH:male:human_20_male_cardiovascular_system-pulmonary_valve_ID';
 const explicit=api.readViewUrl('?hide='+encodeURIComponent(oldValve),current,{...base,visible:['cardiac']});
 for(const id of parts.get(oldValve)?.replacementParts??[])assert.equal(api.partVisible(parts.get(id),explicit,current.parts),false,'explicit valve hide follows its cusps');
 for(const retired of previous.parts.filter(p=>p.system==='cardiac'&&!p.suppressed)){
  assert.equal(parts.get(retired.id)?.suppressed,true,retired.id+' is retained only for rollback');
  assert.equal(api.partVisible(parts.get(retired.id),after,current.parts),false,retired.id+' cannot duplicate the source heart');
 }
 assert.equal(api.localTa98LayerCatalogue('local-ta98',query,false),undefined,'Local donor assets never enter the public route');
});
