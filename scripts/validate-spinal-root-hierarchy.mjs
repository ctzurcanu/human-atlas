import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';
import {SPINAL_LEVELS,ROOT_MAP,rootLevel,buildSpinalRootHierarchy} from './build-spinal-root-hierarchy.mjs';

const input=await loadChakraInputs(),raw=JSON.parse(readFileSync('public/assets/dermatomes-myotomes.json','utf8')),hierarchy=input.parseGuestHierarchy(raw);
assert.deepEqual(raw,buildSpinalRootHierarchy(input),'The deployed file must match the curated map and current catalogues');
assert.equal(hierarchy.nodes.filter(node=>/^ROOT:(?:C|T|L|S|Co)\d+$/.test(node.id)).length,31);
assert.equal(hierarchy.nodes.filter(node=>/^ROOT:(?:C|T|L|S|Co)\d+:(left|right)$/.test(node.id)).length,62);
assert.equal(hierarchy.nodes.filter(node=>node.id.startsWith('DERMATOME:')).length,62);
assert.equal(hierarchy.nodes.filter(node=>node.id.startsWith('MYOTOME:')).length,62);
const flatten=roots=>{const map=new Map(),visit=node=>{if(map.has(node.id))return;map.set(node.id,node);node.children.forEach(visit);};roots.forEach(visit);return map;};
const bundle=await build({stdin:{contents:"export * from './app/hierarchical-explosion';export * from './app/hierarchy-navigation';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
for(const {file,atlas} of input.catalogues){
 const start=performance.now(),roots=input.resolveGuestHierarchy(atlas,hierarchy),nodes=flatten(roots),available=atlas.parts.filter(part=>!part.suppressed),ids=new Set(available.map(part=>part.id));
 assert.deepEqual(new Set(nodes.get('ANATOMY').parts.map(part=>part.id)),ids,`${file}: full anatomy branch must include every available mesh`);
 assert.deepEqual(new Set(nodes.get('ANATOMY').children.flatMap(node=>node.parts.map(part=>part.id))),ids,`${file}: every available mesh must be browsable under its system`);
 assert.deepEqual(new Set(nodes.get('INNERVATION').parts.map(part=>part.id)),ids,`${file}: known or unresolved innervation must retain every available mesh`);
 assert(nodes.get('INNERVATION:UNRESOLVED')?.description.includes('No curated supply mapping'));
 for(const level of SPINAL_LEVELS)for(const side of ['left','right']){
  for(const prefix of ['ROOT','NERVE','DERMATOME','MYOTOME'])assert(nodes.has(`${prefix}:${level}:${side}`),`${file}: missing ${level} ${side} ${prefix}`);
  for(const part of nodes.get(`MYOTOME:${level}:${side}`).parts){assert.equal(part.system,'muscular');assert(!/\bTendon\b/i.test(part.name));assert(!part.name.endsWith(`(${side==='left'?'right':'left'})`),`${file}: contralateral muscle ${part.name}`);}
  for(const part of nodes.get(`NERVE:${level}:${side}`).directParts){assert.equal(part.system,'nervous');assert.equal(rootLevel(atlas.concepts.find(concept=>concept.id===part.conceptId).name),level,`${file}: incorrect root ${part.name}`);}
  // True dermatome geometry is absent; a sensory node must not map the entire skin.
  assert.equal(nodes.get(`DERMATOME:${level}:${side}`).directParts.length,0);
 }
 for(const node of nodes.values()){assert.equal(node.parts.length,new Set(node.parts.map(part=>part.id)).size);for(const part of node.parts)assert(ids.has(part.id));}
 const entries=input.hierarchyEntries(atlas);for(const entry of entries)if(entry.parts.some(part=>part.system==='nervous')&&rootLevel(entry.name)){
  const level=rootLevel(entry.name),side=entry.name.match(/\((left|right)\)$/)?.[1];if(!side)continue;
  for(const part of entry.parts.filter(part=>!part.suppressed))assert(nodes.get(`NERVE:${level}:${side}`).directParts.some(member=>member.id===part.id),`${file}: omitted named root ${entry.name}`);
 }
 assert.equal(nodes.get('CRANIAL').children.length,12,'All twelve cranial branches remain inspectable');
 assert.equal(nodes.get('AUTONOMIC').children.length,4,'Autonomic divisions and mixed plexuses remain explicit');
 for(const root of ['CRANIAL','AUTONOMIC'])for(const part of nodes.get(root).parts)assert(ids.has(part.id));
 const nav=api.hierarchyNavigation(atlas,'guest:dermatomes-myotomes',hierarchy);if(roots.some(node=>node.parts.length))assert(nav.children.length);
 const layout=api.createHierarchicalExplosionLayout(atlas,entries,ids,'guest',1.7,roots);
 for(const stage of layout.stages){const members=stage.clusters.flat();assert.equal(members.length,available.length);assert.equal(new Set(members).size,available.length,'Overlapping roots must not duplicate geometry');assert(stage.positions.every(Number.isFinite));}
 console.log(`${file}: 31 levels / 62 sides, ${available.length} meshes assigned once, ${Math.round(performance.now()-start)} ms`);
}
// All spinal levels must resolve independently when a model provides named roots.
const parts=SPINAL_LEVELS.flatMap(level=>['left','right'].map(side=>({id:`nerve:${level}:${side}`,conceptId:`nerve:${level}:${side}`,name:`${level} root (${side})`,system:'nervous',bounds:[[-.1,.8,0],[.1,.9,.1]],groups:[]})));
const atlas={...input.catalogues[0].atlas,parts,concepts:parts.map(part=>({id:part.conceptId,name:part.name,elements:[part.id]}))},nodes=flatten(input.resolveGuestHierarchy(atlas,hierarchy));
for(const level of SPINAL_LEVELS)for(const side of ['left','right'])assert.deepEqual(nodes.get(`NERVE:${level}:${side}`).parts.map(part=>part.id),[`nerve:${level}:${side}`]);
assert.match(ROOT_MAP.sensory.C1,/no cutaneous dermatome/);assert.match(ROOT_MAP.sensory.S4,/S4–5/);assert.match(ROOT_MAP.sensory.Co1,/not a separate ISNCSCI/);
assert(!hierarchy.nodes.some(node=>node.id.startsWith('MUSCLE:Co1:')),'Do not invent an isolated Co1 myotome');
assert.equal(rootLevel('Vertebra C5'),undefined);assert.equal(rootLevel('C5 segment of cervical spinal cord'),undefined);assert.equal(rootLevel('Cranial nerve V'),undefined);assert.equal(rootLevel('C9 root'),undefined);
const state={selected:[],visible:[],guestQuery:'ROOT:Co1',explode:0,isolate:false,view:'front',rotate:false,reset:0};
const url=new URL(api.viewUrl('http://localhost:3016','local-reference',state));url.searchParams.set('tree','guest:dermatomes-myotomes');assert.equal(api.readViewUrl(url.search,atlas,state).guestQuery,'ROOT:Co1');
console.log('Every spinal level maps independently and bilaterally; sensory geometry gaps and overlapping myotomes stay explicit.');
