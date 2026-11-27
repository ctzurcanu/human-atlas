import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/preferred-source-representations';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const path='/local-models/geometry-workbench/biodigital-male-circulation-v1/hand-current-palm-fit-v981/ordinary-atlas.json';
const read=p=>JSON.parse(readFileSync(p));
const atlas=read('.local-models/'+path.slice('/local-models/'.length));
const base=read('.local-models/organ-models/spinal-cord/male-body.json');
const source=read('.local-models/geometry-workbench/biodigital-male-circulation-v1/atlas-cardiovascular-review-v37/atlas.json');

test('ordinary male routes include both detailed hands while public and source diagnostics remain unchanged',()=>{
 for(const query of ['', 'bowel=layers','bowel=rectal'])assert.equal(api.localTa98LayerCatalogue('local-ta98',new URLSearchParams(query),true),path);
 for(const query of ['bowel=source','delivery=0'])assert.equal(api.localTa98LayerCatalogue('local-ta98',new URLSearchParams(query),true),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',new URLSearchParams(),false),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams(),true),'/local-models/organ-models/spinal-cord/female-body.json');
});

test('the repaired left palm retains original topology and materials; other imported and body anatomy stay literal',()=>{
 const marker=atlas.handLocalInstallation,parts=new Map(source.parts.map(p=>[p.id,p]));
 assert.equal(marker.importedParts.length,52);
 assert.deepEqual(atlas.parts.slice(0,base.parts.length),base.parts);
 for(const id of marker.importedParts){
  const installed=atlas.parts.find(p=>p.id===id),original=parts.get(id);
  if(atlas.handCurrentPalmFit.arterialParts.includes(id)){
   const geometry=new Set(['chunk','positions','normals','indices','normalSigns','geometryOrigin','bounds','uvs','colors']);
   const metadata=p=>Object.fromEntries(Object.entries(p).filter(([k])=>!geometry.has(k)));
   assert.deepEqual(metadata(installed),metadata(original),id);
   const bytes=p=>readFileSync('.local-models/'+p.url.slice('/local-models/'.length));
   const a=bytes(atlas.chunks[installed.chunk]),b=bytes(source.chunks[original.chunk]);
   assert.deepEqual(a.subarray(installed.indices,installed.indices+installed.indexCount*4),b.subarray(original.indices,original.indices+original.indexCount*4),id);
  }else{
   assert.deepEqual({...installed,chunk:original.chunk},original,id);
   assert.deepEqual(atlas.chunks[installed.chunk],source.chunks[original.chunk],id);
  }
 }
 const ids=new Set(marker.importedParts);
 assert.equal(atlas.parts.filter(p=>ids.has(p.id)&&p.system==='arterial').length,34);
 assert.equal(atlas.parts.filter(p=>ids.has(p.id)&&p.system==='venous').length,18);
 assert.equal(marker.anatomicalAcceptance,false);
 assert.equal(marker.arterialInterfaces,90);
});

test('the actual packed fitted left palm has literal shared coordinates and effective normals at all 45 interfaces',()=>{
 const candidate=read('.local-models/geometry-workbench/biodigital-male-circulation-v1/hand-connected-palm-source-repair-v979/report.json');
 const cache=new Map();
 const decode=id=>{
  if(cache.has(id))return cache.get(id);
  const p=atlas.parts.find(p=>p.id===id),bank=readFileSync('.local-models/'+atlas.chunks[p.chunk].url.slice('/local-models/'.length));
  const row=i=>Array.from({length:3},(_,j)=>bank.readFloatLE(p.positions+(i*3+j)*4)+p.geometryOrigin[j]);
  const normal=i=>Array.from({length:3},(_,j)=>bank.readInt16LE(p.normals+(i*3+j)*2)*bank.readInt8(p.normalSigns+i));
  const value={row,normal};cache.set(id,value);return value;
 };
 let count=0;
 for(const pair of candidate.pairs){
  if(!atlas.handCurrentPalmFit.arterialParts.includes(pair.firstPart))continue;
  const a=decode(pair.firstPart),b=decode(pair.secondPart);
  for(let i=0;i<pair.firstRaw.length;i++){
   assert.deepEqual(a.row(pair.firstRaw[i]),b.row(pair.secondRaw[i]));
   assert.deepEqual(a.normal(pair.firstRaw[i]),b.normal(pair.secondRaw[i]));
  }
  count++;
 }
 assert.equal(count,45);assert.equal(atlas.handCurrentPalmFit.sourceNetworkSelfPairs,0);
 assert.equal(atlas.handCurrentPalmFit.anatomicalAcceptance,false);
});

test('detailed palms are visible and recoverable old alternatives cannot duplicate them during Show All',()=>{
 const display=api.withHiddenReviewAlternatives(atlas),parts=new Map(display.parts.map(p=>[p.id,p]));
 const state={selected:[],visible:['arterial','venous'],view:'front',skinOpacity:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:false,explode:0,isolate:false};
 for(const id of atlas.handLocalInstallation.importedParts)assert.equal(api.partVisible(parts.get(id),state,display.parts),true,id);
 for(const id of atlas.handLocalInstallation.preferredOver)assert.equal(api.partVisible(parts.get(id),state,display.parts),false,id);
 const binding=atlas.ta98TermBindings.filter(b=>b.representation==='source-identity'&&b.partIds.every(id=>atlas.handLocalInstallation.importedParts.includes(id)));
 assert.equal(binding.length,48);
 assert.equal(new Set(binding.map(b=>b.term)).size,24);
 for(const b of binding)assert.ok(atlas.concepts.some(c=>c.id===b.conceptId&&b.partIds.every(id=>c.elements.includes(id))));
});

test('ordinary palms use the already reviewed native nerve and muscle sources without modifying their meshes',()=>{
 const display=api.withHiddenReviewAlternatives(atlas),parts=new Map(display.parts.map(p=>[p.id,p]));
 for(const key of ['handSourceSelection','handAdductorSourceSelection']){
  assert.deepEqual(atlas[key],source[key]);
  for(const id of atlas[key].hiddenReviewAlternatives)assert.equal(parts.get(id).suppressed,true,id);
  for(const id of atlas[key].preferred){
   assert.ok(!parts.get(id).suppressed,id);
   assert.deepEqual(atlas.parts.find(p=>p.id===id),base.parts.find(p=>p.id===id),id);
  }
 }
});
