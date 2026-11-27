import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {readFile} from 'node:fs/promises';
import {resolve} from 'node:path';
import {gunzipSync} from 'node:zlib';
import {catalogueAttribution, projectCatalogueAttribution} from './catalogue-attribution.mjs';

const root=resolve(import.meta.dirname,'..'), site=resolve(root,process.argv[2]??'outputs/male-ta98-pages');
const original=JSON.parse(await readFile(resolve(root,'.local-models/geometry-workbench/dmta-v1/dmtm-v8-skin-clearance-v1/circulation-atlas.json')));
const atlas=JSON.parse(await readFile(resolve(site,'models/male-ta98/atlas.json')));
const expected=projectCatalogueAttribution(original,await catalogueAttribution(resolve(root,'ATTRIBUTION.md')));
assert.equal(atlas.name,'Male TA98');
assert.equal(atlas.parts.length,expected.parts.length);
assert.deepEqual(atlas.concepts,expected.concepts);
assert.deepEqual(atlas.ta98TermBindings,expected.ta98TermBindings);
assert.deepEqual(atlas.defaultHidden,expected.defaultHidden);
assert.deepEqual(atlas.sourceGeometry.chunks,expected.chunks);
const widths={positions:12,normals:6,indices:4,uvs:8,colors:3,normalSigns:1,tissueChart:6,gastricLumen:1};
const pairs=Map.groupBy(atlas.parts.map((part,i)=>({part,before:expected.parts[i]})),pair=>pair.before.chunk);
const systemCounts=new Map();
const verified=new Set();
for (const [index,rows] of pairs) {
 const url=original.chunks[index].url.split('?')[0];
 const raw=await readFile(resolve(root,url.startsWith('/local-models/')?'.local-models/'+url.slice(14):'public'+url));
 const packed=new Map();
 for (const {part,before} of rows) {
  const {chunk,...sourceLayout}=atlas.sourceGeometry.parts[part.id];
  assert.equal(chunk,index);
  assert.equal(part.id,before.id);
  for (const key of ['name','system','conceptId','ta98Term','ta98Kind','bounds','suppressed','indexCount','vertexCount']) assert.deepEqual(part[key],before[key],part.id);
  systemCounts.set(part.system,(systemCounts.get(part.system)??0)+1);
  const bank=atlas.chunks[part.chunk];
  assert.equal(bank.system,part.system);
  assert.ok(bank.gzip.startsWith('/models/male-ta98/geometry/'+part.system+'/'));
  if (!packed.has(part.chunk)) {
   const bytes=gunzipSync(await readFile(resolve(site,'.'+bank.gzip)));
   assert.equal(bytes.length,bank.bytes);
   assert.equal(createHash('sha256').update(bytes).digest('hex'),bank.sha256);
   packed.set(part.chunk,bytes); verified.add(part.chunk);
  }
  for (const [key,width] of Object.entries(widths)) {
   assert.equal(sourceLayout[key],before[key]);
   assert.equal(key in before,key in part);
   if (!(key in part)) continue;
   const length=(key==='indices'?part.indexCount:part.vertexCount)*width;
   const target=packed.get(part.chunk).subarray(part[key],part[key]+length);
   assert.equal(target.length,length);
   assert.ok(target.equals(raw.subarray(before[key],before[key]+length)),part.id+' '+key);
  }
 }
}
assert.equal(verified.size,atlas.chunks.length);
assert.deepEqual(Object.fromEntries(systemCounts),Object.fromEntries(atlas.systems.map(s=>[s.id,s.items])));
console.log('Verified original geometry and labels for',atlas.parts.length,'parts in',systemCounts.size,'systems.');
