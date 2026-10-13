import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
const read=p=>JSON.parse(fs.readFileSync(p));
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
test('Actual female banks share nine junction contours, area-weighted normals and reciprocal removable tissue rims',()=>{
 const base='.local-models/ta98-review/bowel-wall-main-prototype/female',atlas=read(base+'/atlas.json'),report=read('reports/ta98-female/small-intestine-main-prototype-pack.json'),blob=fs.readFileSync(base+'/wall.bin'),source=read('.local-models/ta98-review/small-intestine-layer-v1/manifest.json'),raw=fs.readFileSync('.local-models/ta98-review/small-intestine-layer-v1/wall.bin');
 assert.equal(sha(blob),report.payloadSha256);assert.equal(sha(raw),report.sourceWallSha256);assert.equal(blob.length,5623404,'no added mesh payload');assert.equal(report.joinedWallPorts.length,2);
 const byId=new Map(atlas.parts.map(p=>[p.id,p]));
 const organ=id=>{const spec=source.organs.find(o=>o.id===id),packed=report.organs.find(o=>o.organ===id);return {...spec,...packed};};
 const normalSum=(o,k,i)=>{const n=[0,0,0],pos=v=>[0,1,2].map(d=>blob.readFloatLE(o.positionBank+((k*o.vertices+v)*3+d)*4));for(let f=0;f<o.triangles;f++){const ids=[0,1,2].map(d=>raw.readUInt32LE(o.faces+(f*3+d)*4));if(!ids.includes(i))continue;const [a,b,c]=ids.map(pos),u=b.map((v,d)=>v-a[d]),v=c.map((v,d)=>v-a[d]),cross=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];for(let d=0;d<3;d++)n[d]+=cross[d];}return n;};
 for(const join of report.joinedWallPorts){const [A,B]=join.organs.map(organ);assert.equal(join.pairs.length,12);
  for(let k=0;k<9;k++)for(const [a,b] of join.pairs){const sa=normalSum(A,k,a),sb=normalSum(B,k,b),sum=sa.map((v,d)=>v+sb[d]),len=Math.hypot(...sum);for(let d=0;d<3;d++){const pa=A.positionBank+((k*A.vertices+a)*3+d)*4,pb=B.positionBank+((k*B.vertices+b)*3+d)*4;assert.deepEqual(blob.subarray(pa,pa+4),blob.subarray(pb,pb+4));const na=blob.readInt16LE(A.normalBank+((k*A.vertices+a)*3+d)*2),nb=blob.readInt16LE(B.normalBank+((k*B.vertices+b)*3+d)*2);assert.equal(na,nb);assert.equal(na,Math.round(sum[d]/len*32767));}}
  for(const [own,other] of [[A,B],[B,A]]){let annularIndices=0;const controls=new Set(join.pairs.map(p=>p[own===A?0:1]));for(const id of own.newLayerParts){const p=byId.get(id);assert.equal(p.tissueDrawGroups.reduce((n,g)=>n+g.count,0),p.indexCount);for(const g of p.tissueDrawGroups){const q=byId.get(g.neighbor);if(!q||!other.newLayerParts.includes(q.id))continue;annularIndices+=g.count;assert.ok(q.tissueDrawGroups.some(h=>h.neighbor===p.id));for(let at=g.start;at<g.start+g.count;at+=6){const ids=Array.from({length:6},(_,d)=>blob.readUInt32LE(p.indices+(at+d)*4));assert.equal(new Set(ids.map(i=>i%own.vertices)).size,2);assert.ok(ids.every(i=>controls.has(i%own.vertices)));assert.equal(new Set(ids.map(i=>Math.floor(i/own.vertices))).size,2,'annulus only, never lumen disk');}}}assert.equal(annularIndices,576);}
 }
 assert.equal(atlas.bowelWallPrototype.anatomicallyAccepted,false);
});
test('Unregistered male source mouths are never welded or hidden against a different organ',()=>{
 const atlas=read('.local-models/ta98-review/bowel-wall-main-prototype/male/atlas.json'),report=read('reports/ta98-male/small-intestine-main-prototype-pack.json'),parts=new Map(atlas.parts.map(p=>[p.id,p]));assert.deepEqual(report.joinedWallPorts,[]);
 for(const p of atlas.parts.filter(p=>p.bowelWallSource))for(const g of p.tissueDrawGroups)if(g.neighbor)assert.equal(parts.get(g.neighbor).bowelWallSource,p.bowelWallSource);
});
