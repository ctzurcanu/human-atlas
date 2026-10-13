import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
import {readFile} from 'node:fs/promises';
const bundle=await build({stdin:{contents:"export {createPartitionedTissueWall} from './app/partitioned-tissue-wall';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {createPartitionedTissueWall}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
function closed(g){const edges=new Map();for(let i=0;i<g.index.count;i+=3){const tri=[0,1,2].map(j=>g.index.getX(i+j));for(let j=0;j<3;j++){const a=tri[j],b=tri[(j+1)%3],key=a<b?`${a}:${b}`:`${b}:${a}`,row=edges.get(key)??[0,0];row[0]++;row[1]+=a<b?1:-1;edges.set(key,row);}}assert.ok([...edges.values()].every(([count,winding])=>count===2&&winding===0));}
test('Removing half a middle coat exposes only its actual two neighboring patches and the remaining tissue edge',()=>{
 const faces=Uint32Array.of(0,1,2,0,2,3),rims=Uint32Array.of(0,1,1,2,2,3,3,0),interfaces=[0,-1,-2,-3].map(z=>Float32Array.of(0,0,z,1,0,z,1,1,z,0,1,z));
 const wall=createPartitionedTissueWall({interfaces,faces,rims,layerIds:['cover','muscle','inner']},[Uint16Array.of(0,0),Uint16Array.of(1,2),Uint16Array.of(3,3)],[{id:'cover',layer:0},{id:'band',layer:1},{id:'remainder',layer:1},{id:'inner',layer:2}]);
 wall.full.forEach(closed);const originalIndices=wall.displayed.map(g=>g.index);
 wall.update([true,true,true,true]);assert.equal(wall.displayed.reduce((n,g)=>n+g.groups.reduce((n,c)=>n+c.count,0),0)/3,28);
 wall.update([true,false,true,true]);assert.equal(wall.displayed[1].drawRange.count,0);
 function surface(owner,material){const g=wall.displayed[owner],result=[];for(const c of g.groups.filter(c=>c.materialIndex===material))for(let i=c.start;i<c.start+c.count;i++)result.push(g.index.getX(i)%4);return result;}
 assert.deepEqual(surface(0,1),[0,2,1],'only the removed band triangle on cover underside');
 assert.deepEqual(surface(3,0),[0,1,2],'only the removed band triangle on inner coat exterior');
 assert.deepEqual(surface(2,0),[]);assert.deepEqual(surface(2,1),[]);
 assert.equal(wall.displayed[2].groups.filter(c=>c.materialIndex===2).reduce((n,c)=>n+c.count,0)/3,6,'two outer edge quads plus one newly exposed diagonal edge quad');
 wall.full.forEach(closed);wall.displayed.forEach((g,i)=>assert.equal(g.index,originalIndices[i]));
});
for(const sex of ['female','male'])test(`${sex}: eight coats and bands share one bank; all 44 section-counting shells remain closed during partial peeling`,async()=>{
 const folder='.local-models/ta98-review/colon-integrated-wall/',m=JSON.parse(await readFile(folder+sex+'-manifest.json','utf8')),bytes=await readFile(folder+sex+'-wall.bin'),raw=bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength),interfaces=Array.from({length:9},(_,i)=>new Float32Array(raw,m.interfaces+i*m.vertices*12,m.vertices*3)),r=new Uint8Array(raw,m.regionOwners,m.triangles),b=new Uint8Array(raw,m.bandOwners,m.triangles);
 const subsets=[],ownership=m.layers.map(()=>new Uint16Array(m.triangles));
 for(let region=0;region<4;region++)for(let layer=0;layer<8;layer++)for(const band of layer===2?[0,1,2,3]:[-1]){const owner=subsets.length;subsets.push({id:`${region}:${layer}:${band}`,layer});for(let fi=0;fi<m.triangles;fi++)if(r[fi]===region&&(band<0||b[fi]===band))ownership[layer][fi]=owner;}
 const wall=createPartitionedTissueWall({interfaces,faces:new Uint32Array(raw,m.faces,m.triangles*3),rims:new Uint32Array(raw,m.rimsOffset,m.rims*2),layerIds:m.layers.map(l=>l.id)},ownership,subsets);assert.equal(wall.full.length,44);
 wall.full.forEach(closed);for(const g of [...wall.full,...wall.displayed]){assert.equal(g.attributes.position,wall.position);assert.equal(g.attributes.normal,wall.normal);}
 const mask=Array(44).fill(true);wall.update(mask);const band=subsets.findIndex(s=>s.id==='1:2:1');mask[band]=false;wall.update(mask);
 const above=subsets.findIndex(s=>s.id==='1:1:-1'),below=subsets.findIndex(s=>s.id==='1:3:-1');
 for(const [owner,material] of [[above,1],[below,0]]){const expected=wall.contacts[owner].filter(c=>c.neighbor===band&&c.material===material);assert.equal(expected.length,1);assert.ok(expected[0].count>0);assert.ok(wall.displayed[owner].groups.some(g=>g.start===expected[0].start&&g.count===expected[0].count));}
 assert.equal(wall.displayed[band].drawRange.count,0);wall.full.forEach(closed);
});
