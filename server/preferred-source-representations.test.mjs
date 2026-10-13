import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/preferred-source-representations'; export * from './app/viewer-state'; export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const read=sex=>JSON.parse(readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
const decisions=JSON.parse(readFileSync('app/data/preferred-source-representations.json'));
const tree=api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json')));

test('explicit preferred source assemblies replace old alternatives without changing geometry or camera',()=>{
 for(const sex of ['male','female']){
  const original=read(sex),a=api.withPreferredSourceRepresentations(original);
  const byId=new Map(a.parts.map(p=>[p.id,p]));
  assert.equal(a.parts.length,original.parts.length);
  assert.equal(a.chunks,original.chunks);assert.equal(a.materials,original.materials);
  for(const d of decisions[sex]){
   assert.equal(byId.get(d.id).suppressed,true);
   assert.deepEqual(byId.get(d.id).replacementParts,d.parts);
   assert.deepEqual(new Set(api.resolveSelection(a,[d.id])),new Set(d.parts));
   assert.deepEqual(new Set(api.resolveSelection(a,[d.conceptId])),new Set(d.parts));
   assert.ok(!a.defaultHidden.includes(d.id));
   assert.ok(d.parts.every(id=>!a.defaultHidden.includes(id)));
  }
  for(let k=0;k<a.parts.length;k++){
   const {suppressed,replacementParts,...now}=a.parts[k],{suppressed:_,replacementParts:__,...before}=original.parts[k];assert.deepEqual(now,before);
  }
  const nodes=api.resolveGuestHierarchy(a,tree),visible=[];const walk=n=>{visible.push(...n.directParts.map(p=>p.id));n.children.forEach(walk);};nodes.forEach(walk);
  assert.ok(decisions[sex].every(d=>!visible.includes(d.id)));
  assert.equal(api.withPreferredSourceRepresentations(a),a,'loading is idempotent');
  const base={focus:0,contextOpacity:1,skinOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:false,explode:0,visible:['muscular'],selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0};
  const camera=[.5,.6,.7,.4,.3,.2,.5,.45,.6,.5,1,.5];
  assert.deepEqual(api.readViewUrl('?frame='+camera.join(',')+'&select='+encodeURIComponent('FJ1555'),a,base).camera,camera);
 }
});

test('ambiguous generic tissues, unresolved alternatives and distinct layers are retained',()=>{
 const a=api.withPreferredSourceRepresentations(read('male'));
 for(const id of ['FJ1581','FJ1458M','FJ1549','FJ1549M'])assert.equal(a.parts.find(p=>p.id===id).suppressed,undefined);
 assert.ok(a.parts.filter(p=>p.pairedTissue).every(p=>!decisions.male.some(d=>d.id===p.id)));
 const bad=read('male'),d=decisions.male[0];bad.parts=bad.parts.filter(p=>p.id!==d.parts[0]);
 assert.ok(!api.withPreferredSourceRepresentations(bad).parts.find(p=>p.id===d.id).suppressed,'an incomplete replacement cannot retire its source');
 const publicAtlas={...read('male'),reviewOnly:false};assert.equal(api.withPreferredSourceRepresentations(publicAtlas),publicAtlas);
});

test('both demo bodies keep every pending source hidden after Show All, saved empty visibility and parent selection',()=>{
 for(const sex of ['male','female']){
  const source=read(sex),preferred=api.withPreferredSourceRepresentations(source),a=api.withHiddenReviewAlternatives(preferred),byId=new Map(a.parts.map(p=>[p.id,p]));
  const alternatives=preferred.defaultHidden.filter(id=>byId.has(id));
  assert.ok(alternatives.length>700);
  assert.ok(alternatives.every(id=>byId.get(id).suppressed));
  const base={focus:0,contextOpacity:1,skinOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:false,explode:0,visible:[...new Set(a.parts.map(p=>p.system))],selected:[],isolate:false,view:'right',rotate:false,reset:0};
  const restored=api.readViewUrl('?hide=nonexistent&skin=1',a,base);
  for(const id of alternatives){
   assert.equal(api.partVisible(byId.get(id),restored,a.parts),false);
   assert.equal(api.partVisible(byId.get(id),{...restored,selected:[id],isolate:true},a.parts),false);
   assert.equal(api.sectionPartVisible(byId.get(id),restored),false);
  }
  const roots=api.resolveGuestHierarchy(a,tree);const ids=new Set();const walk=n=>{n.parts.forEach(p=>ids.add(p.id));n.children.forEach(walk);};roots.forEach(walk);
  assert.ok(alternatives.every(id=>!ids.has(id)));
  assert.equal(api.withHiddenReviewAlternatives(a),a);
  assert.equal(api.withHiddenReviewAlternatives(preferred,true),preferred,'unfinished sources remain accessible in an explicit repair view');
  assert.equal(a.chunks,source.chunks);assert.equal(a.parts.length,source.parts.length);
  const native=source.parts.filter(p=>!p.suppressed&&!source.defaultHidden.includes(p.id));assert.ok(native.every(p=>byId.get(p.id)===p));
 }
});
