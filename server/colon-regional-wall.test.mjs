import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export {createRegionalTissueWall} from './app/regional-tissue-wall';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {createRegionalTissueWall}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const folder='.local-models/ta98-review/colon-layer-v1/',m=JSON.parse(await readFile(folder+'manifest.json','utf8')),bytes=await readFile(folder+'wall.bin'),raw=bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength);
const outer=new Float32Array(raw,m.outer,m.vertices*3),inner=new Float32Array(raw,m.inner,m.vertices*3),interfaces=[outer];let t=0;
for(const l of m.layers){t+=l.fraction;interfaces.push(Float32Array.from(outer,(x,i)=>x+(inner[i]-x)*Math.min(t,1)));}interfaces[8]=inner;
for(const p of m.interfaceCorrections)interfaces[p.boundary].set(p.position,p.vertex*3);
const wall=createRegionalTissueWall({interfaces,faces:new Uint32Array(raw,m.faces,m.triangles*3),rims:new Uint32Array(raw,m.rimEdges,m.rims*2),layerIds:m.layers.map(l=>l.id)},new Uint8Array(raw,m.regionOwners,m.triangles),m.regions);
test('All four source partitions share coordinates/normals without losing tissue faces',()=>{
 assert.deepEqual(wall.regionFaceCounts,[2064,3360,2688,1202]);assert.equal(wall.regionFaceCounts.reduce((a,b)=>a+b),m.triangles);
 for(const g of [...wall.full,...wall.displayed]){assert.equal(g.attributes.position,wall.position);assert.equal(g.attributes.normal,wall.normal);}
 assert.equal(wall.position.count,m.vertices*9);
 for(let r=0;r<4;r++)for(let k=0;k<8;k++){
  const g=wall.full[r*8+k],edges=new Map();
  for(let j=0;j<g.index.count;j+=3){const f=[0,1,2].map(a=>g.index.getX(j+a));for(let a=0;a<3;a++){const x=f[a],y=f[(a+1)%3],key=x<y?`${x}:${y}`:`${y}:${x}`,row=edges.get(key)??[0,0];row[0]++;row[1]+=x<y?1:-1;edges.set(key,row);}}
  assert.ok([...edges.values()].every(([count,winding])=>count===2&&winding===0),'each region coat is a closed counting shell');
  for(let j=wall.regionFaceCounts[r]*6;j<g.index.count;j+=3){const controls=[0,1,2].map(a=>g.index.getX(j+a)%m.vertices);assert.ok(new Set(controls).size<=2,'rim spans wall thickness only');}
 }
});
test('Assembled regional joins have no displayed barriers; hiding neighbours exposes annuli',()=>{
 const visible=Array(32).fill(true),indices=wall.displayed.map(g=>g.index);wall.update(visible);
 const rimCount=r=>wall.displayed[r*8].groups.filter(g=>g.materialIndex===2).reduce((n,g)=>n+g.count,0);
 assert.equal(rimCount(0),24*6);assert.equal(rimCount(1),0);assert.equal(rimCount(2),0);assert.equal(rimCount(3),24*6);
 visible.fill(false,8,16);wall.update(visible);assert.equal(rimCount(0),48*6);assert.equal(rimCount(2),24*6);
 for(let k=8;k<16;k++)assert.equal(wall.displayed[k].drawRange.count,0);
 wall.displayed.forEach((g,k)=>assert.equal(g.index,indices[k]));
 for(let r=0;r<4;r++)for(let k=0;k<8;k++)assert.equal(wall.full[r*8+k].groups.length,2+wall.ports[r].length);
});
