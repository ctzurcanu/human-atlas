import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/guest-choice';export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {parseGuestHierarchy,resolveGuestHierarchy,guestDirectAnatomyParts,guestNodeChoice}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const hierarchy=parseGuestHierarchy(JSON.parse(await readFile('public/assets/ta98.json','utf8')));
for(const sex of ['female','male'])test(`${sex}: wall tree retains tissue choices and missing terms without repeating inherited mesh rows`,async()=>{
 const atlas=JSON.parse(await readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`,'utf8'));
 const roots=resolveGuestHierarchy(atlas,hierarchy),nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};roots.forEach(walk);
 for(const organ of ['duodenum','jejunum','ileum']){
  for(const key of ['muscularis','mucosa']){
   const node=nodes.get(`SIW:${organ}:${key}`);assert.ok(node);
   assert.equal(guestDirectAnatomyParts('ta98',node).length,0);
   assert.deepEqual(new Set(guestNodeChoice('ta98',node).elements),new Set(node.parts.map(p=>p.id)));
  }
  const leaf=nodes.get(`SIW:${organ}:circular-muscle`);assert.equal(guestDirectAnatomyParts('ta98',leaf).length,1);
  for(const key of ['circular-folds','villi','intestinal-glands']){
   const missing=nodes.get(`SIW:${organ}:${key}`);assert.ok(missing);assert.equal(missing.parts.length,key==='circular-folds'?4:0);
   assert.ok(guestNodeChoice('ta98',nodes.get(`SIW:${organ}:mucosa`)).children.some(c=>c.id.endsWith(missing.id)));
  }
 }
});
