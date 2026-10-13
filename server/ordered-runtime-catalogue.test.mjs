import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
import {orderedRuntimeCatalogue} from './ordered-runtime-catalogue.mjs';
const read=p=>JSON.parse(fs.readFileSync(p)),sha=b=>crypto.createHash('sha256').update(b).digest('hex');
test('Ordinary male atlas retains source URL indices and exact runtime buffer references for every delivered source mesh',()=>{
 const raw=fs.readFileSync('.local-models/ta98-male.json'),source=JSON.parse(raw),runtime=read('.local-models/ta98-runtime.json'),atlas=read('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'),r=read('reports/ta98-male/small-intestine-local-atlas-pack.json'),ordered=orderedRuntimeCatalogue(source,runtime,sha(raw)),packed=new Map(runtime.parts.map(p=>[p.id,p]));
 assert.equal(atlas.delivery.mode,'on-demand');assert.equal(atlas.delivery.sourceAtlasSha256,sha(raw));assert.equal(r.sourceDeliverySha256,sha(fs.readFileSync('.local-models/ta98-runtime.json')));
 assert.deepEqual(atlas.parts.slice(0,source.parts.length).map(p=>p.id),source.parts.map(p=>p.id));assert.deepEqual(ordered.concepts,source.concepts);assert.deepEqual(ordered.defaultHidden,source.defaultHidden);assert.deepEqual(ordered.materials,source.materials);
 const fields=['positions','normals','indices','uvs','colors','normalSigns','gastricLumen'];let live=0,historical=0;
 for(let i=0;i<source.parts.length;i++){
  const p=atlas.parts[i],original=source.parts[i],delivered=packed.get(p.id);
  if(delivered){live++;for(const field of fields)assert.equal(p[field],delivered[field],p.id+':'+field);for(const field of ['vertexCount','indexCount','bounds','geometryOrigin','normalSigns','sourceOffset','material','tissueDrawGroups'])assert.deepEqual(p[field],delivered[field],p.id+':'+field);assert.deepEqual(atlas.chunks[p.chunk],runtime.chunks[delivered.chunk]);}
  else{historical++;assert.equal(p.suppressed,true);for(const field of fields)assert.equal(p[field],original[field]);assert.deepEqual(atlas.chunks[p.chunk],source.chunks[original.chunk]);}
 }
 assert.equal(live,6301);assert.equal(historical,177);
 const active=new Set(atlas.parts.filter(p=>!p.suppressed).map(p=>p.chunk));assert.equal(active.size,r.deliveryResources.reachableChunks);assert.equal([...active].reduce((n,i)=>n+atlas.chunks[i].bytes,0),r.deliveryResources.reachableGeometryBytes);assert.ok(r.deliveryResources.reachableGeometryBytes<r.deliveryResources.sourceGeometryBytes*.3);
 for(let i=runtime.chunks.length;i<atlas.chunks.length-1;i++)assert.ok(!active.has(i),'historical source buffers must never be requested');
 assert.equal(atlas.parts.filter(p=>p.pairedTissue&&p.id.startsWith('GASTRIC:male')).length,35);assert.equal(atlas.bowelWallPrototype.partIds.length,36);
});
test('Ordered delivery rejects stale, incomplete or unrelated runtime catalogues',()=>{
 const source={parts:[{id:'one',chunk:0},{id:'old',chunk:0,suppressed:true}],chunks:[{url:'original'}],materials:{}},runtime={parts:[{id:'one',chunk:0}],chunks:[{url:'packed'}],materials:{},delivery:{mode:'on-demand',sourceAtlasSha256:'current'}};
 assert.throws(()=>orderedRuntimeCatalogue(source,runtime,'stale'),/bound/);
 assert.throws(()=>orderedRuntimeCatalogue(source,{...runtime,parts:[]},'current'),/omitted/);
 assert.throws(()=>orderedRuntimeCatalogue(source,{...runtime,parts:[...runtime.parts,{id:'unrelated'}]},'current'),/unrelated/);
 assert.throws(()=>orderedRuntimeCatalogue(source,{...runtime,materials:{changed:true}},'current'),/materials/);
 const result=orderedRuntimeCatalogue(source,runtime,'current');assert.deepEqual(result.parts.map(p=>p.id),['one','old']);assert.equal(result.parts[1].chunk,1);assert.deepEqual(source.chunks,[{url:'original'}]);
});
