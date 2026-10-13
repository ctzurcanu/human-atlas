import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import crypto from 'node:crypto';
import {build} from 'esbuild';
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const built=await build({entryPoints:['app/guest-hierarchy.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {parseGuestHierarchy,resolveGuestHierarchy}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const hierarchy=parseGuestHierarchy(JSON.parse(await fs.readFile('public/assets/ta98.json','utf8')));
const nodes=roots=>{const by=new Map();const walk=n=>{by.set(n.id,n);n.children.forEach(walk);};roots.forEach(walk);return by;};
for(const sex of ['female','male'])test(`${sex} local wall catalogue exposes actual finite contextual tissues while leaving duodenal covering masks missing`,async()=>{
 const root=`.local-models/ta98-review/bowel-wall-main-prototype/${sex}`,atlas=JSON.parse(await fs.readFile(root+'/atlas.json','utf8')),proof=JSON.parse(await fs.readFile(`reports/ta98-${sex}/small-intestine-main-prototype-pack.json`,'utf8'));
 assert.equal(sha(await fs.readFile(`.local-models/ta98-${sex}.json`)),proof.sourceCatalogueSha256);
 assert.equal(sha(await fs.readFile(root+'/wall.bin')),proof.payloadSha256);
 assert.equal(new Set(atlas.parts.map(p=>p.id)).size,atlas.parts.length);
 assert.equal(new Set(atlas.concepts.map(c=>c.id)).size,atlas.concepts.length);
 const parts=atlas.parts.filter(p=>atlas.bowelWallPrototype.partIds.includes(p.id));assert.equal(parts.length,36);
 const graph=nodes(resolveGuestHierarchy(atlas,hierarchy));
 for(const organ of ['duodenum','jejunum','ileum']){
  for(const key of ['longitudinal-muscle','circular-muscle','submucosa','muscularis-mucosae','lamina-propria','epithelium']){
   const node=graph.get(`SIW:${organ}:${key}`);assert.equal(node.parts.length,['epithelium','lamina-propria','muscularis-mucosae','submucosa'].includes(key)?2:1);assert.doesNotMatch(node.description,/Missing separate representation/);assert.match(node.description,/not anatomically accepted/);
  }
  assert.equal(graph.get(`SIW:${organ}:muscularis`).parts.length,2);assert.equal(graph.get(`SIW:${organ}:mucosa`).parts.length,7);
 }
 for(const key of ['serosa','subserosa','adventitia'])assert.deepEqual(graph.get(`SIW:duodenum:${key}`).parts,[]);
 for(const organ of ['jejunum','ileum'])for(const key of ['serosa','subserosa'])assert.equal(graph.get(`SIW:${organ}:${key}`).parts.length,1);
 for(const key of ['villi','intestinal-glands'])assert.deepEqual(graph.get(`SIW:jejunum:${key}`).parts,[]);
 for(const entry of proof.organs){const cell=parts.filter(p=>entry.newLayerParts.includes(p.id));assert.equal(new Set(cell.map(p=>p.positions)).size,1);assert.equal(new Set(cell.map(p=>p.normals)).size,1);for(const id of entry.oldSourceParts)assert.equal(atlas.parts.find(p=>p.id===id).suppressed,true);}
 assert.equal(atlas.bowelWallPrototype.anatomicallyAccepted,false);
});
