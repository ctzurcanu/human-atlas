import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const bundle=await build({stdin:{contents:"export * from './app/anatomy-hierarchy'; export * from './app/anatomical-components';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const path='.local-models/female-ta98-complete/atlas.json';
test('each side is a muscle assembly with its own intersections', {skip:!fs.existsSync(path)},()=>{
 const atlas=JSON.parse(fs.readFileSync(path,'utf8'));
 const parts=new Map(atlas.parts.map(p=>[p.id,p]));
 const entries=api.hierarchyEntries(atlas);
 const nodes=api.buildAnatomyNodes(entries.filter(e=>e.system==='muscular'),e=>api.systemPathFor(e,atlas.scope));
 const find=(nodes,name)=>{for(const n of nodes){if(n.kind==='group'&&n.name===name)return n;if(n.kind==='group'){const hit=find(n.nodes,name);if(hit)return hit;}}};
 for(const side of ['left','right']){
  const carrier=`ZRA:${side}:rectus_abdominis`,tendon=`ZRA:${side}:tendinous_intersections`;
  const concept=atlas.concepts.find(c=>c.id===carrier),parent=atlas.concepts.find(c=>c.id===concept.partOf);
  assert.equal(parent.id,`ZRA:group:${side}_rectus_abdominis_side`);
  assert.equal(atlas.concepts.find(c=>c.id===tendon).partOf,parent.id);
  const group=find(nodes,`Rectus abdominis (${side})`);
  assert.ok(group);assert.deepEqual(group.parts.map(p=>p.id).sort(),[carrier,tendon].sort());
  assert.ok(group.nodes.some(n=>n.kind==='entry'&&n.entry.id===tendon));
  assert.deepEqual(api.anatomicalComponentIds(concept,[],parts,atlas.concepts),[tendon]);
  assert.equal(api.anatomicalComponentRelations(concept,[],parts,atlas.concepts).find(r=>r.kind==='contains')?.target.id,tendon);
 }
});
