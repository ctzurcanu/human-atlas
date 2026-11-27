import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile,stat} from 'node:fs/promises';
import {gunzipSync} from 'node:zlib';
import {resolve,sep} from 'node:path';
import {createHash} from 'node:crypto';

const root=resolve('public/models/female-ta98');
const base='/models/female-ta98/';
function file(url){
 assert.ok(url.startsWith(base),`Unexpected asset URL: ${url}`);
 const result=resolve(root,url.slice(base.length).split('?')[0]);
 assert.ok(result.startsWith(root+sep),'Asset path escapes the model folder');
 return result;
}
test('the female catalogue references complete, size-limited system assets',async()=>{
 const atlas=JSON.parse(await readFile(resolve(root,'atlas.json'),'utf8'));
 const systems=atlas.chunks.map(()=>new Set());
 for(const part of atlas.parts){
  assert.ok(Number.isInteger(part.chunk)&&systems[part.chunk],part.id);
  systems[part.chunk].add(part.system);
  const chunk=atlas.chunks[part.chunk];
  for(const [key,width] of Object.entries({positions:12,normals:6,indices:4,uvs:8,normalSigns:1,tissueChart:6,gastricLumen:1,colors:3})){
   if(part[key]===undefined)continue;
   const n=(key==='indices'?part.indexCount:part.vertexCount)*width;
   assert.ok(Number.isInteger(part[key])&&part[key]>=0&&part[key]+n<=chunk.bytes,`${part.id}: ${key}`);
  }
  if(part.hairField){
   const raw=await readFile(file(part.hairField.url));
   assert.equal(raw.length,part.hairField.bytes);
   assert.deepEqual(gunzipSync(await readFile(file(part.hairField.url)+'.gz')),raw);
  }
 }
 assert.ok(systems.every(s=>s.size===1));
 for(const chunk of atlas.chunks){
  const path=file(chunk.url),info=await stat(path);
  assert.equal(info.size,chunk.bytes);assert.ok(info.size<50_000_000);
  const digest=createHash('sha256').update(await readFile(path)).digest('hex');
  assert.equal(digest,chunk.sha256);
 }
 const urls=new Set();
 function visit(value){if(typeof value==='string'&&value.startsWith(base))urls.add(value);else if(Array.isArray(value))value.forEach(visit);else if(value&&typeof value==='object')Object.values(value).forEach(visit);}
 visit(atlas);
 await Promise.all([...urls].map(url=>stat(file(url))));
});
