import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/guest-hierarchy.ts','app/guest-choice.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/small-intestine-hierarchy-tests'});
const module=async name=>import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles.find(f=>f.path.endsWith(name+'.js')).contents).toString('base64'));
const {parseGuestHierarchy,resolveGuestHierarchy}=await module('guest-hierarchy');
const {guestNodeChoice}=await module('guest-choice');
const raw=JSON.parse(await readFile('public/assets/ta98.json','utf8')),hierarchy=parseGuestHierarchy(raw),graph=new Map(hierarchy.nodes.map(n=>[n.id,n]));
const plan=JSON.parse(await readFile('app/data/small-intestine-wall-plan.json','utf8'));
const findNodes=roots=>{const found=new Map();const visit=n=>{if(found.has(n.id))return;found.set(n.id,n);n.children.forEach(visit);};roots.forEach(visit);return found;};
test('Organ-specific mucosa, muscle and exterior structure retains canonical terminal ileum and authority links',()=>{
 for(const organ of plan.organs){
  assert.deepEqual(graph.get(organ.ta98).childIds.filter(id=>id.startsWith('SIW:')),[...plan.innerToOuter,...(organ.outerAlternatives??[])].map(k=>`SIW:${organ.id}:${k}`));
  assert.deepEqual(graph.get(`SIW:${organ.id}:muscularis`).childIds,[`SIW:${organ.id}:circular-muscle`,`SIW:${organ.id}:longitudinal-muscle`]);
  for(const key of ['epithelium','lamina-propria','muscularis-mucosae','circular-folds','villi','intestinal-glands'])assert.ok(graph.get(`SIW:${organ.id}:mucosa`).childIds.includes(`SIW:${organ.id}:${key}`));
 }
 assert.ok(graph.get('A05.6.04.001').childIds.includes('A05.6.04.002'));
 assert.ok(graph.get('SIW:duodenum:submucosa').childIds.includes('SIW:duodenum:duodenal-glands'));
 assert.ok(graph.get('SIW:ileum:mucosa').childIds.includes('SIW:ileum:aggregated-nodules'));
 assert.ok(graph.has('SIW:duodenum:serosa')&&graph.has('SIW:duodenum:adventitia'));
 assert.ok(!graph.has('SIW:jejunum:adventitia'));
 for(const key of ['epithelium','lamina-propria','adventitia'])assert.ok(!graph.get(`SIW:duodenum:${key}`).synonyms.some(x=>x.startsWith('TA98:')));
});
test('Neither main atlas falsely exposes an outer bowel carrier as a represented wall layer',async()=>{
 for(const file of ['.local-models/ta98-runtime.json','.local-models/ta98-female-runtime.json']){
  const atlas=JSON.parse(await readFile(file,'utf8')),nodes=findNodes(resolveGuestHierarchy(atlas,hierarchy));
  for(const [id,node] of nodes)if(id.startsWith('SIW:')){assert.deepEqual(node.parts,[]);assert.match(node.description,/Missing separate representation/);}
  for(const code of ['A05.6.02.001','A05.6.03.001','A05.6.04.001'])assert.ok(nodes.get(code).parts.length);
 }
});
test('A represented layer selector cannot leak between organs or into adventitial covering',async()=>{
 const atlas=JSON.parse(await readFile('.local-models/ta98-runtime.json','utf8'));
 const id='SIW:duodenum:serosa',part={...atlas.parts[0],id:'fixture:duodenal-serosa',conceptId:id,suppressed:false,system:'digestive'};
 const fixture={...atlas,parts:[part],concepts:[{id,name:'Duodenal serosa',elements:[part.id]}]};
 const nodes=findNodes(resolveGuestHierarchy(fixture,hierarchy));
 assert.deepEqual(nodes.get(id).parts.map(p=>p.id),[part.id]);
 for(const other of ['SIW:jejunum:serosa','SIW:ileum:serosa','SIW:duodenum:adventitia'])assert.deepEqual(nodes.get(other).parts,[]);
 assert.equal(guestNodeChoice('ta98',nodes.get('SIW:duodenum:epithelium')).terminology.ta98,null);
});
