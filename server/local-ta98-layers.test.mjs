import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/viewer-state';export * from './app/view-part-codes';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {localTa98LayerCatalogue,readViewUrl,partFingerprint,packPartIds,compatiblePartAtlas}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const read=p=>JSON.parse(fs.readFileSync(p));
for(const sex of ['female','male'])test(`${sex} ordinary local atlas retains the installed intestinal coats and restores camera intent`,()=>{
 const model=sex==='male'?'local-ta98':'local-female-ta98',path=localTa98LayerCatalogue(model,new URLSearchParams(),true),atlas=read('.local-models/'+path.replace('/local-models/','')),prior=atlas,report=read(`reports/ta98-${sex}/small-intestine-local-atlas-pack.json`);
 assert.ok(atlas.bowelWallPrototype.ordinaryLocalAtlas);assert.equal(report.partIds.length,36);assert.ok(report.partIds.every(id=>atlas.bowelWallPrototype.partIds.includes(id)));
 const old=report.organs.find(o=>o.organ==='jejunum').oldSourceParts[0],camera=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5],q=new URLSearchParams({p:partFingerprint(prior),s:packPartIds([old],prior,false),frame:camera.join(','),cut:'coronal',slice:'.227'}),base={selected:[],visible:[],view:'three-quarter'};
 const restored=readViewUrl(q.toString(),atlas,base);assert.deepEqual(new Set(restored.selected),new Set(report.organs.find(o=>o.organ==='jejunum').newLayerParts));assert.deepEqual(restored.camera,camera);assert.equal(restored.section.position,.227);
 assert.equal(localTa98LayerCatalogue(model,new URLSearchParams('bowel=source'),true),undefined);assert.equal(localTa98LayerCatalogue(model,new URLSearchParams('delivery=0'),true),undefined);assert.equal(localTa98LayerCatalogue(model,new URLSearchParams(),false),undefined);
});
test('The user’s existing female URL fingerprint still resolves its exact hidden structures',()=>{
 const atlas=read('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'),prior=compatiblePartAtlas(atlas,'260.12sxffm.uizlzc');assert.ok(prior);assert.equal(partFingerprint(prior),'260.12sxffm.uizlzc');
 const q='p=260.12sxffm.uizlzc&h=jr.js.la-ms.mv-o3.o5.o7-oz.p1-p3.p5-px.pz-y8.153-15a.15c-15k.15m-161.163-169.1lq-1ls.1mb-1mw.1wk-1yh&l=136&skin=0&cut=coronal&slice=0.227&flip=1&frame=0.432952,0.603637,0.984524,0.565065,0.281339,0.505912,0.39865,0.480956,0.205037,0.5,1,0.5';
 const base={selected:[],visible:[],view:'three-quarter'},old=readViewUrl(q,prior,base),now=readViewUrl(q,atlas,base);assert.ok(old.hidden.length>0);assert.deepEqual(now.camera,old.camera);assert.deepEqual(now.visible,old.visible);for(const id of old.hidden){const p=atlas.parts.find(p=>p.id===id);for(const next of p.replacementParts??[id])assert.ok(now.hidden.includes(next));}
});
