import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
import {spinalCordTissues} from './reference-spinal-cord.mjs';

const atlas=JSON.parse(readFileSync('.local-models/reference.json'));
const bundle=await build({entryPoints:['app/anatomy-hierarchy.ts','app/depth-layers.ts','app/viewer-state.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/reference-validation'});
const load=async name=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(name)).text).toString('base64')}`);
const {hierarchyEntries}=await load('anatomy-hierarchy.js');
const {depthLayerFor}=await load('depth-layers.js');
const {partVisible}=await load('viewer-state.js');
const available=atlas.parts.filter(part=>!part.suppressed),entries=hierarchyEntries(atlas);
assert.equal(entries.flatMap(entry=>entry.parts).length,available.length,'Hierarchy must contain every renderable part once');
assert(entries.every(entry=>entry.parts.every(part=>!part.suppressed)),'Suppressed aliases cannot enter group checkboxes');
assert.equal(available.filter(part=>depthLayerFor(part)==='skin').length,702,'All derived skin and topographic skin regions must share Skin');
assert(atlas.parts.filter(part=>part.system==='regions').every(part=>depthLayerFor(part)==='skin'),'Topographic skin cannot fall under annotations');
assert(atlas.parts.filter(part=>part.sourceId==='Anterior region of thigh.l').every(part=>part.regions?.includes('lower-left')),'Thigh skin must follow its source lower-limb assignment');
const connective=available.filter(part=>part.system==='connective');
assert(connective.length>0&&connective.every(part=>!part.suppressed),'Connective group must reach fully-on state');
const suppressed=atlas.parts.find(part=>part.suppressed);
assert(suppressed);
assert.equal(partVisible(suppressed,{selected:[suppressed.id],visible:[suppressed.system],hidden:[],depthHidden:[],skinOpacity:0,explode:0,isolate:false,view:'back',rotate:false,reset:0}),false,'A selected suppressed alias must stay hidden');

function positions(name){
 const part=atlas.parts.find(item=>item.name===name);assert(part,`Missing ${name}`);
 const bytes=readFileSync(`.local-models/reference/chunk-${part.chunk}.bin`);
 return {part,vertices:Array.from({length:part.vertexCount},(_,index)=>[0,1,2].map(axis=>bytes.readFloatLE(part.positions+index*12+axis*4)))};
}
assert.equal(spinalCordTissues.size,27,'The source cord extrusion family changed');
const dura=positions('Spinal dura').vertices,bins=new Map();
for(const vertex of dura){const key=Math.round(vertex[1]*100),members=bins.get(key)??[];members.push(vertex);bins.set(key,members);}
const envelopes=new Map();
for(let key=105;key<=155;key++){
 const ring=[...(bins.get(key-1)??[]),...(bins.get(key)??[]),...(bins.get(key+1)??[])];
 if(ring.length)envelopes.set(key,[Math.min(...ring.map(point=>point[0])),Math.max(...ring.map(point=>point[0])),Math.min(...ring.map(point=>point[2])),Math.max(...ring.map(point=>point[2]))]);
}
for(const name of spinalCordTissues){
 const {part,vertices}=positions(name);
 assert(part.bounds[0][1]>1.10&&part.bounds[1][1]<1.53,`${name} extends beyond the cord's modeled canal`);
 assert(part.regions?.includes('torso'),`${name} must be available in torso sections`);
 for(const vertex of vertices){
  const envelope=envelopes.get(Math.round(vertex[1]*100));
  assert(envelope,`${name}: no dura near y=${vertex[1]}`);
  const [minX,maxX,minZ,maxZ]=envelope;
  assert(vertex[0]>=minX-.001&&vertex[0]<=maxX+.001&&vertex[2]>=minZ-.001&&vertex[2]<=maxZ+.001,`${name} escapes the modeled spinal dura`);
 }
}
assert.equal(atlas.parts.filter(part=>part.id.startsWith('REF:core:')&&part.system==='nervous'&&part.bounds[0][1]<1.4&&part.bounds[1][1]>1.7).length,0,'A long misplaced core nerve mesh remains');
console.log(`Reference import validated: ${available.length} hierarchy pieces, all ${spinalCordTissues.size} cord tissues inside the spinal dura.`);
