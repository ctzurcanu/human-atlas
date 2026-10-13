import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/bowel-wall-peeling';export * from './app/guest-hierarchy';export * from './app/guest-choice';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {bowelWallPeelSteps,parseGuestHierarchy,resolveGuestHierarchy,guestNodeChoice,partVisible}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const hierarchy=parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
for(const sex of ['female','male'])test(sex+' integrated caecocolic hierarchy exposes eight peelable tissues per region',()=>{
 const atlas=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 const hasBands=!!atlas.caecocolicWallIntegration.taeniaPartition;assert.equal(atlas.caecocolicWallIntegration.partIds.length,hasBands?60:48);
 const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 for(const region of ['caecum','appendix','ascending','transverse','descending','sigmoid']){
  const node=nodes.get(`CAECOWALL:${sex}:${region}:wall-assembly`);assert.ok(node);assert.equal(node.children.length,5);
  const choice=guestNodeChoice('ta98',node),steps=bowelWallPeelSteps(atlas,choice);assert.equal(steps?.length,8);assert.equal(choice.elements.length,hasBands&&!['caecum','appendix'].includes(region)?11:8);
  for(let d=0;d<=8;d++){
   const hidden=steps.slice(0,d).flatMap(s=>s.ids),state={selected:choice.elements,hidden,depthHidden:[],visible:['digestive'],isolate:true,region:'all',skinOpacity:0,contextOpacity:0};
   assert.equal(node.parts.filter(p=>partVisible(p,state,atlas.parts)).length,choice.elements.length-hidden.length);
  }
  assert.equal(bowelWallPeelSteps(atlas,{...choice,elements:choice.elements.slice(1)}),null);
 }
 const colon=atlas.concepts.find(c=>c.id===`CAECOWALL:${sex}:colon`);assert.deepEqual(bowelWallPeelSteps(atlas,colon).map(s=>s.ids.length),hasBands?[4,4,16,4,4,4,4,4]:Array(8).fill(4));
 const p=atlas.parts.find(p=>p.id===colon.elements[0]);assert.equal(atlas.parts.find(x=>x.id===p.bowelWallSource).suppressed,true);
 const copy=structuredClone(atlas);copy.parts.find(x=>x.id===colon.elements[0]).positions+=4;assert.equal(bowelWallPeelSteps(copy,colon),null);
});
