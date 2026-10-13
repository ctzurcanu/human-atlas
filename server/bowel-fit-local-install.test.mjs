import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
const atlas=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
const proof=JSON.parse(fs.readFileSync('reports/bowel-local-installed-v86.json'));
const before=JSON.parse(fs.readFileSync(proof.snapshot));
const after=new Map(atlas.parts.map(p=>[p.id,p]));
test('Only the promoted local J/I wall parts change; original anatomy and controls survive',()=>{
 assert.equal(proof.installed,true); assert.equal(proof.partCount,24);
 const changed=new Set(proof.parts);
 for(const part of before.parts)if(!changed.has(part.id))assert.deepEqual(after.get(part.id),part);
 assert.deepEqual(atlas.concepts,before.concepts);assert.deepEqual(atlas.defaultHidden,before.defaultHidden);
 const parts=proof.parts.map(id=>after.get(id));assert.equal(new Set(parts.map(p=>p.chunk)).size,1);assert.equal(new Set(parts.map(p=>p.positions)).size,1);
 for(const part of parts){assert.ok(part.pairedTissue);assert.equal(part.tissueDrawGroups[0].count,part.bowelOuterFaceCount*3);}
});
test('The restored alimentary tour binds the installed geometry and crosses the shared J/I mouth exactly',()=>{
 const route=JSON.parse(fs.readFileSync('.local-models/ta98-review/visits/male.json'));
 assert.equal(route.bridges.find(b=>b.from==='Jejunum'&&b.to==='Ileum').distanceMm,0);
 for(const name of ['Jejunum','Ileum']){
  const segment=route.segments.find(s=>s.name===name),binding=route.bindings.find(b=>b.id===segment.parts[0]),part=after.get(binding.id);
  assert.equal(binding.chunkUrl,atlas.chunks[part.chunk].url);assert.equal(binding.positions,part.positions);assert.equal(binding.indices,part.indices);assert.equal(binding.vertexCount,part.vertexCount);assert.ok(segment.points.length>250);
 }
});
