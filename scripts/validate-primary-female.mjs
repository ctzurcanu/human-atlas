import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import {gunzipSync} from 'node:zlib';

const filename=path.resolve(process.argv[2]??'public/models/atlas-hra-female.json');
const atlas=JSON.parse(fs.readFileSync(filename,'utf8'));
assert.equal(atlas.source,'Human Reference Atlas');
assert.equal(atlas.sex,'female');
assert.equal(atlas.provenance.sourceSha256,'95f0c3d2f918582608692ca1139e8bdb18c147a16470e9ee9af8b276bd77c422');
assert.ok(!/brianp|pridgen|slorksmo/i.test(JSON.stringify(atlas)));
assert.equal(atlas.parts.length,1030);
const embryo=JSON.parse(fs.readFileSync(path.join(path.dirname(filename),'atlas-embryo.json'),'utf8'));
assert.equal(embryo.scope,'embryo');
assert.equal(embryo.parts.length,8);
assert.equal(embryo.parts.some(part=>part.name==='Amnion'),true);
assert.equal(atlas.parts.some(part=>embryo.parts.some(embryoPart=>embryoPart.id===part.id)),false);
const embryoBuffer=fs.readFileSync(path.join(path.dirname(filename),path.basename(embryo.chunks[0].url)));
assert.equal(embryoBuffer.length,embryo.chunks[0].bytes);
assert.ok(gunzipSync(fs.readFileSync(path.join(path.dirname(filename),path.basename(embryo.chunks[0].gzip)))).equals(embryoBuffer));
for(const part of embryo.parts){
 assert.ok(embryo.concepts.some(concept=>concept.elements.includes(part.id)));
 for(const [offset,length] of [[part.positions,part.vertexCount*12],[part.normals,part.vertexCount*6],[part.indices,part.indexCount*4]])assert.ok(offset%4===0&&offset+length<=embryoBuffer.length,part.id);
}
assert.equal(atlas.provenance.supplementParts.length,8);
assert.equal(atlas.provenance.muscleSourceSha256,'9886eda040f6087bbb65182530f8bd262be456f366b155025c12b3dc56da774e');
assert.equal(atlas.provenance.muscleCount,74);
assert.equal(atlas.provenance.muscleArchiveCount,76);
assert.equal(atlas.sourceTriangles,7756626);
const ids=new Set(atlas.parts.map(p=>p.id));assert.equal(ids.size,atlas.parts.length);
const chunks=atlas.chunks.map(chunk=>{const raw=fs.readFileSync(path.join(path.dirname(filename),path.basename(chunk.url)));assert.equal(raw.length,chunk.bytes);const gz=fs.readFileSync(path.join(path.dirname(filename),path.basename(chunk.gzip)));assert.equal(gz.length,chunk.gzipBytes);assert.ok(gunzipSync(gz).equals(raw));return raw;});
let triangles=0;
for(const part of atlas.parts){
 assert.ok(atlas.materials[part.material]);
 assert.ok(part.provenance?.url.includes('humanatlas.io/')||part.provenance?.url.includes('digitalcommons.du.edu/'));
 assert.ok(part.indexCount>0&&part.indexCount%3===0&&part.vertexCount>0);
 const chunk=chunks[part.chunk];
 for(const [offset,length] of [[part.positions,part.vertexCount*12],[part.normals,part.vertexCount*6],[part.indices,part.indexCount*4]]){assert.equal(offset%4,0);assert.ok(offset+length<=chunk.length,part.id);}
 const positions=new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3);
 const indices=new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount);
 for(const value of positions)assert.ok(Number.isFinite(value),part.id);
 for(const value of indices)assert.ok(value<part.vertexCount,part.id);
 triangles+=part.indexCount/3;
}
assert.equal(triangles,atlas.triangles);
for(const concept of atlas.concepts){assert.ok(concept.elements.length);for(const id of concept.elements)assert.ok(ids.has(id),id);}
for(const name of ['Skin of body','Uterus','Ovary (left)','Ovary (right)','Rectus femoris (left)','Rectus femoris (right)'])assert.ok(atlas.concepts.some(c=>c.name===name),name);
assert.ok(!atlas.parts.some(p=>/prostate|testis|penis/i.test(p.name)));
console.log(`Validated ${atlas.parts.length} direct-source female meshes, ${atlas.concepts.length} concepts, ${triangles.toLocaleString()} triangles, ${chunks.length} matching raw/gzip buffers.`);
