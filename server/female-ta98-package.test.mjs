import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {readFile, stat} from 'node:fs/promises';
import {test} from 'node:test';
import {gunzipSync} from 'node:zlib';

test('female TA98 system banks and textures decode with their catalogue offsets', async () => {
 const prefix='/models/female-ta98/';
 const atlas=JSON.parse(await readFile('public'+prefix+'atlas.json','utf8'));
 assert.equal(atlas.parts.length,7136);
 assert.equal(new Set(atlas.parts.map(part=>part.id)).size,atlas.parts.length);
 assert.ok(!JSON.stringify(atlas).includes('/local-models/'));
 const widths={positions:12,normals:6,indices:4,uvs:8,normalSigns:1,tissueChart:6,gastricLumen:1,colors:3};
 for(const [i,chunk] of atlas.chunks.entries()){
  const members=atlas.parts.filter(part=>part.chunk===i);
  assert.equal(new Set(members.map(part=>part.system)).size,1);
  assert.ok(chunk.url.startsWith(prefix+'geometry/'+members[0].system+'-'));
  assert.equal(chunk.gzip,chunk.url);
  const packed=await readFile('public'+chunk.gzip), bank=gunzipSync(packed);
  assert.equal(packed.length,chunk.gzipBytes);
  assert.equal(bank.length,chunk.bytes);
  assert.equal(createHash('sha256').update(bank).digest('hex'),chunk.sha256);
  for(const part of members){
   for(const [key,width] of Object.entries(widths)){
    if(part[key]===undefined)continue;
    const count=key==='indices'?part.indexCount:part.vertexCount;
    assert.ok(Number.isInteger(part[key])&&part[key]>=0&&part[key]+width*count<=bank.length,part.id+' '+key);
   }
   if(part.hairField){
    const field=part.hairField, data=gunzipSync(await readFile('public'+field.url+'.gz'));
    assert.equal(createHash('sha256').update(data).digest('hex'),field.sha256);
   }
   assert.doesNotMatch(part.description??'',/\bSources?\s*:/i);
  }
 }
 for(const part of atlas.parts)assert.ok(atlas.chunks[part.chunk],part.id);
 const textures=new Set(JSON.stringify(atlas).match(/\/models\/female-ta98\/textures\/[^"?]+/g));
 assert.ok(textures.size>0);
 for(const texture of textures)assert.ok((await stat('public'+texture)).isFile(),texture);
});
