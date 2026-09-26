import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
import * as T from 'three';
import {mergeGeometries} from 'three/addons/utils/BufferGeometryUtils.js';
import {spinalCordTissues} from './reference-spinal-cord.mjs';

const atlas=JSON.parse(readFileSync('.local-models/reference.json'));
const bundle=await build({entryPoints:['app/anatomy-hierarchy.ts','app/depth-layers.ts','app/viewer-state.ts','app/section-cap.ts','app/section-cap-profile.ts','app/mesh-topology.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/reference-validation'});
const load=async name=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(name)).text).toString('base64')}`);
const {hierarchyEntries}=await load('anatomy-hierarchy.js');
const {depthLayerFor}=await load('depth-layers.js');
const {partVisible}=await load('viewer-state.js');
const {sectionCapGeometry}=await load('section-cap.js');
const {sectionCapProfile}=await load('section-cap-profile.js');
const {analyzeMeshTopology}=await load('mesh-topology.js');
const available=atlas.parts.filter(part=>!part.suppressed),entries=hierarchyEntries(atlas);
assert.equal(entries.flatMap(entry=>entry.parts).length,available.length,'Hierarchy must contain every renderable part once');
assert(entries.every(entry=>entry.parts.every(part=>!part.suppressed)),'Suppressed aliases cannot enter group checkboxes');
assert.equal(available.filter(part=>depthLayerFor(part)==='skin').length,702,'All derived skin and topographic skin regions must share Skin');
assert(atlas.parts.filter(part=>part.system==='regions').every(part=>depthLayerFor(part)==='skin'),'Topographic skin cannot fall under annotations');
assert(atlas.parts.filter(part=>part.sourceId==='Anterior region of thigh.l').every(part=>part.regions?.includes('lower-left')),'Thigh skin must follow its source lower-limb assignment');
const connective=available.filter(part=>part.system==='connective');
assert(connective.length>0&&connective.every(part=>!part.suppressed),'Connective group must reach fully-on state');
const structuralName=/\b(?:bursa|bursae|iliopectineal arch|retinaculum|aponeurosis|tendon|tendinous|ligament)\b/i;
for(const part of available){
 if(/\.(?:[oei]\d*[lr])$/i.test(part.sourceId))continue;
 if(/\bfascia\b/i.test(part.sourceId)&&!/tensor fasciae latae/i.test(part.sourceId))assert.equal(part.system,'fascia',`${part.sourceId} belongs under Fascia`);
 if(structuralName.test(part.sourceId)&&!/\bnode of\b/i.test(part.sourceId))assert.equal(part.system,'connective',`${part.sourceId} belongs under Connective tissue`);
}
for(const part of available.filter(item=>/\bmuscle\.[lr]$/i.test(item.sourceId)&&!/\b(?:tendon|ligament|bursa|bursae|fascia|retinaculum|aponeurosis|nerve|artery|vein|node|branch)\b/i.test(item.sourceId))){
 assert.equal(part.system,'muscular',`${part.sourceId} must be a muscle, not an attachment overlay`);
 const material=atlas.materials?.[part.material];
 if(material&&!material.map){assert(material.color[0]>=material.color[1],`${part.sourceId} retains the teal source overlay color`);assert.equal(material.opacity,undefined,`${part.sourceId} remains translucent`);}
}
assert.equal(available.filter(part=>part.system==='muscular'&&!atlas.materials?.[part.material]?.map&&atlas.materials?.[part.material]?.color?.[1]>atlas.materials?.[part.material]?.color?.[0]).length,0,'Green untextured muscle material remains');
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
// Material primitives of a single source muscle/lung have open boundaries on
// their own. Joining the intersected pieces must yield a solid cut surface.
for(const [concept,height] of [['REF:Rectus abdominis muscle.l',1.05257],['REF:Superior lobe of right lung',1.355]]){
 const fragments=atlas.parts.filter(part=>part.conceptId===concept&&part.bounds[0][1]<=height&&part.bounds[1][1]>=height).map(part=>{
  const bytes=readFileSync(`.local-models/reference/chunk-${part.chunk}.bin`);
  const geometry=new T.BufferGeometry();
  geometry.setAttribute('position',new T.BufferAttribute(new Float32Array(bytes.buffer,bytes.byteOffset+part.positions,part.vertexCount*3),3));
  geometry.setIndex(new T.BufferAttribute(new Uint32Array(bytes.buffer,bytes.byteOffset+part.indices,part.indexCount),1));
  return {part,geometry};
 });
 assert(fragments.length>1,`${concept}: expected material primitives`);
 const joined=mergeGeometries(fragments.map(item=>item.geometry));assert(joined);
 const topology=analyzeMeshTopology(joined),profile=sectionCapProfile(fragments[0].part,topology);
 assert(!profile.thinShell,`${concept}: material pieces still appear to be an open sheet`);
 const cap=sectionCapGeometry(joined,new T.Plane(new T.Vector3(0,-1,0),height),new T.Matrix4(),profile);
 assert(cap?.userData.fillIndexCount>0,`${concept}: joined cross-section has no solid cap`);
 cap.dispose();joined.dispose();for(const item of fragments)item.geometry.dispose();
}
console.log(`Reference import validated: ${available.length} hierarchy pieces, all ${spinalCordTissues.size} cord tissues inside the spinal dura.`);
