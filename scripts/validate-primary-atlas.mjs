/** Validate an independently built atlas before installing it in public/models.
 * Usage: node scripts/validate-primary-atlas.mjs /path/to/atlas-z-anatomy.json
 */
import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import {gunzipSync} from 'node:zlib';

const filename=path.resolve(process.argv[2]??'public/models/atlas-z-anatomy.json');
const atlas=JSON.parse(fs.readFileSync(filename,'utf8'));
assert.equal(atlas.source,'Z-Anatomy + Open 3D Model');
assert.match(atlas.provenance.archiveSha256,/^[a-f0-9]{64}$/);
assert.match(atlas.provenance.sourceFileSha256,/^[a-f0-9]{64}$/);
assert.equal(atlas.provenance.additionalSources?.[0]?.archiveSha256,'5af0190a6d7bf47393447ac30021e4f3ba619721c7f3a620c39a895947078432');
assert.equal(atlas.provenance.additionalSources?.[0]?.sourceFileSha256,'e440c84c794239d1850e62b4ede0195d81bcff4f0078c7945528c388ab72fdb4');
assert.equal(atlas.parts.length,5254);
assert.equal(atlas.parts.filter(part=>part.id.startsWith('O3M:')).length,492);
for(const [name,system] of [
 ['Ureter.l','urinary'],['Ureter.r','urinary'],
 ['Common iliac artery.l','arterial'],['Common iliac vein.r','venous'],
 ['Femoral artery.l','arterial'],['Femoral vein.r','venous'],
 ['Sciatic nerve.l','nervous'],['Femoral nerve.r','nervous'],
 ['Pudendal nerve.l','nervous'],['Tibial nerve.r','nervous'],
]){
 const part=atlas.parts.find(part=>part.sourceId===name);
 assert.ok(part,`Missing original source curve ${name}`);
 assert.equal(part.system,system,name);
}
assert.ok(atlas.triangles>5_000_000);
assert.ok(!JSON.stringify(atlas).match(/brianp|pridgen|slorksmo/i));
const ids=new Set(atlas.parts.map(p=>p.id));
assert.equal(ids.size,atlas.parts.length);
const chunks=atlas.chunks.map(chunk=>{
 const raw=fs.readFileSync(path.join(path.dirname(filename),path.basename(chunk.url)));
 assert.equal(raw.length,chunk.bytes);
 const compressed=fs.readFileSync(path.join(path.dirname(filename),path.basename(chunk.gzip)));
 assert.equal(compressed.length,chunk.gzipBytes);
 assert.ok(gunzipSync(compressed).equals(raw));
 return raw;
});
let triangles=0;
for(const part of atlas.parts){
 assert.ok(part.provenance?.url.includes(part.id.startsWith('O3M:')?'anatomytool.org/open3dmodel-create':'github.com/Z-Anatomy/'));
 assert.ok(atlas.materials[part.material]);
 assert.ok(['CC BY-SA 4.0','CC BY-NC 4.0','CC BY-NC-SA 4.0'].includes(part.license));
 assert.ok(part.indexCount>0&&part.indexCount%3===0&&part.vertexCount>0);
 const chunk=chunks[part.chunk];
 for(const [offset,count,bytes] of [[part.positions,part.vertexCount*3,4],[part.normals,part.vertexCount*3,2],[part.indices,part.indexCount,4]]){
  assert.equal(offset%4,0);
  assert.ok(offset+count*bytes<=chunk.length,part.id);
 }
 const positions=new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3);
 const indices=new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount);
 for(const number of positions)assert.ok(Number.isFinite(number),part.id);
 for(const index of indices)assert.ok(index<part.vertexCount,part.id);
 for(let axis=0;axis<3;axis++)assert.ok(part.bounds[0][axis]<=part.bounds[1][axis]);
 triangles+=part.indexCount/3;
}
assert.equal(triangles,atlas.triangles);
assert.equal(new Set(atlas.concepts.map(concept=>concept.id)).size,atlas.concepts.length);
for(const concept of atlas.concepts){assert.ok(concept.elements.length);for(const id of concept.elements)assert.ok(ids.has(id),id);}
for(const label of ['Stomach','Kidney.l','Kidney.r','Sternocostal head of pectoralis major muscle.l']){
 assert.ok(atlas.parts.some(p=>p.sourceId===label),label);
}
assert.equal(atlas.parts.filter(p=>p.license==='CC BY-NC 4.0').length,2);
assert.equal(atlas.parts.filter(p=>p.license==='CC BY-NC-SA 4.0').length,4);
assert.equal(atlas.parts.filter(p=>p.id.startsWith('O3M:')&&p.provenance.mirrored).length,246);
console.log(`Validated ${atlas.parts.length} direct-source surfaces, ${atlas.concepts.length} concepts, ${triangles.toLocaleString()} triangles, ${chunks.length} matching raw/gzip buffers.`);
