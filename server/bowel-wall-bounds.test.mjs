import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/bowel-wall-bounds.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {bowelWallBoundingSphere}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
for(const sex of ['female','male'])test(`${sex}: stored-cell spheres enclose every actual indexed wall vertex`,async()=>{
 const root=`.local-models/ta98-review/bowel-wall-main-prototype/${sex}`,atlas=JSON.parse(await fs.readFile(root+'/atlas.json')),raw=Uint8Array.from(await fs.readFile(root+'/wall.bin')).buffer,parts=atlas.parts.filter(p=>p.bowelWallSource);assert.equal(parts.length,36);
 let bypassed=0,used=0;
 for(const part of parts){const sphere=bowelWallBoundingSphere(part);assert.ok(sphere);assert.ok(sphere.radius>0);const p=new Float32Array(raw,part.positions,part.vertexCount*3),ids=new Uint32Array(raw,part.indices,part.indexCount),referenced=new Set(ids);
  for(const id of referenced){const v=Array.from(p.slice(id*3,id*3+3));assert.ok(Math.hypot(...v.map((x,d)=>x-sphere.center.getComponent(d)))<=sphere.radius,`${part.id}: vertex ${id} outside picking sphere`);for(let d=0;d<3;d++)assert.ok(v[d]+part.geometryOrigin[d]>=part.bounds[0][d]-1e-10&&v[d]+part.geometryOrigin[d]<=part.bounds[1][d]+1e-10);}
  used+=referenced.size;bypassed+=part.vertexCount*2;
 }
 assert.ok(used<bypassed/2);assert.equal(bypassed,sex==='female'?1725192:3903552);
});
test('ordinary, incomplete, translated or malformed bounds keep the ordinary geometry path',async()=>{
 const a=JSON.parse(await fs.readFile('.local-models/ta98-review/bowel-wall-main-prototype/female/atlas.json')),part=a.parts.find(p=>p.bowelWallSource);
 for(const p of [{...part,bowelWallSource:undefined},{...part,pairedTissue:false},{...part,geometryOrigin:undefined},{...part,tissueDrawGroups:undefined},{...part,sourceOffset:[0,0,0]},{...part,bounds:[[NaN,0,0],[1,1,1]]},{...part,bounds:[[2,0,0],[1,1,1]]},{...part,bounds:[[0,0],[1,1,1]]}])assert.equal(bowelWallBoundingSphere(p),null);
 for(const p of a.parts.filter(p=>!p.bowelWallSource))assert.equal(bowelWallBoundingSphere(p),null);
});
