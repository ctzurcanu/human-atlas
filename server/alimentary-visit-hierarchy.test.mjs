import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const compiled=await build({stdin:{contents:"export {resolveGuestHierarchy,parseGuestHierarchy} from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {resolveGuestHierarchy,parseGuestHierarchy}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const hierarchy=parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));

for(const sex of ['male','female'])test(`${sex}: every alimentary visit landmark is selectable under the alimentary system`,()=>{
 const atlas=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 const route=JSON.parse(fs.readFileSync(`.local-models/ta98-review/visits/${sex}.json`));
 const all=new Map();const walk=node=>{all.set(node.id,node);node.children.forEach(walk)};
 resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 const alimentary=all.get('A05.0.00.000');assert.ok(alimentary);
 const branch=new Map();const collect=node=>{branch.set(node.id,node);node.children.forEach(collect)};collect(alimentary);
 const represented=new Set(alimentary.parts.map(p=>p.id));
 for(const segment of route.segments){
  for(const id of segment.parts??[]){
   assert.ok(atlas.parts.some(p=>p.id===id&&!p.suppressed),`${segment.name}: active landmark ${id}`);
   assert.ok(represented.has(id),`${segment.name}: ${id} belongs under alimentary system`);
   assert.ok([...branch.values()].some(n=>n.directParts.some(p=>p.id===id)),`${segment.name}: ${id} has a selectable node`);
  }
 }
 if(sex==='male'){
  const occurrence=branch.get('ALIMENTARY:anal-canal:external-sphincter');assert.ok(occurrence);
  const canonical=all.get('A04.5.04.012');assert.deepEqual(occurrence.directParts.map(p=>p.id).sort(),canonical.directParts.map(p=>p.id).sort());
  assert.equal(occurrence.directParts.length,2);
  assert.equal(all.get('A05.7.05.001').directParts.length,0,'a sphincter reference does not fabricate an anal-canal wall');
 }
});
