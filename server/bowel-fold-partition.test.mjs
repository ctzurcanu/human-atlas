import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
import {partitionFoldShell} from '../scripts/bowel-feature-partition.mjs';
const compiled=await build({entryPoints:['web/bowel-fold-fields.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {buildBowelInterfaces}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const components=await build({entryPoints:['app/anatomical-components.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalComponentIds,storedAssemblyParent}=await import('data:text/javascript;base64,'+Buffer.from(components.outputFiles[0].contents).toString('base64'));
for(const sex of ['female','male'])test(`${sex}: fold footprint is a disjoint closed mucosal/submucosal partition on unchanged banks`,async()=>{
 const root=`.local-models/ta98-review/bowel-wall-main-prototype/${sex}`,atlas=JSON.parse(await readFile(root+'/atlas.json')),bank=Uint8Array.from(await readFile(root+'/wall.bin')).buffer;
 const source='.local-models/ta98-review/small-intestine-layer-v1'+(sex==='male'?'-male':''),meta=JSON.parse(await readFile(source+'/manifest.json')),raw=Uint8Array.from(await readFile(source+'/wall.bin')).buffer,hr=Uint8Array.from(await readFile(source+'/fold-height.bin')).buffer;
 for(const organ of meta.organs){
  const n=organ.vertices,faces=new Uint32Array(raw,organ.faces,organ.triangles*3),height=new Float32Array(hr,meta.foldRelief.organs.find(o=>o.id===organ.id).offset,n),p=partitionFoldShell(faces,height,n);
  const parts=atlas.parts.filter(p=>p.bowelWallSource&&p.id.includes(`:BOWELWALL:${organ.id}:`));
  const fractions=[0];for(const l of meta.layers)fractions.push(fractions.at(-1)+l.fraction);fractions[8]=1;
  const expected=buildBowelInterfaces(new Float32Array(raw,organ.outer,n*3),new Float32Array(raw,organ.inner,n*3),fractions,meta.foldRelief.gains,height),actual=new Float32Array(bank,parts[0].positions,n*9*3);
  for(let k=0;k<9;k++)assert.deepEqual(actual.subarray(k*n*3,(k+1)*n*3),expected[k]);
  for(const attr of ['positions','normals','tissueChart'])assert.equal(new Set(parts.map(p=>p[attr])).size,1);
  const visible=new Set(parts.map(p=>p.id));let assembledIndices=0;for(const part of parts)for(const group of part.tissueDrawGroups)if(!visible.has(group.neighbor))assembledIndices+=group.count;
  assert.equal(assembledIndices,organ.triangles*6+organ.rims*6*8);
  const key=t=>Array.from(t).join(','),original=new Set(Array.from({length:faces.length/3},(_,i)=>key(faces.subarray(i*3,i*3+3)))),union=new Set();
  for(let layer=4;layer<8;layer++){
  const partition=partitionFoldShell(faces,height,n,layer,layer+1),layerName=meta.layers[layer].id;const layerUnion=new Set();
  for(let which=0;which<2;which++){
   const patch=partition.patches[which],part=atlas.parts.find(part=>part.bowelWallSource&&part.id.endsWith(`:${organ.id}:`+(which?'circular-folds'+(layer===7?'':':'+layerName):layerName)));
   assert.ok(part);assert.deepEqual(new Uint32Array(bank,part.indices,part.indexCount),patch.indices);
   for(let j=0;j<patch.faces.length;j+=3){const t=key(patch.faces.subarray(j,j+3));assert.ok(original.has(t));assert.ok(!layerUnion.has(t));layerUnion.add(t);}
   const edges=new Map();for(let j=0;j<patch.indices.length;j+=3)for(let d=0;d<3;d++){const a=patch.indices[j+d],b=patch.indices[j+(d+1)%3],k=Math.min(a,b)+':'+Math.max(a,b),e=edges.get(k)??[0,0];e[0]++;e[1]+=a<b?1:-1;edges.set(k,e);}
   for(const e of edges.values())assert.deepEqual(e,[2,0]);
   // Closed edges alone can hide two solids touching at a bow-tie vertex.
   const links=new Map();for(let j=0;j<patch.indices.length;j+=3)for(let d=0;d<3;d++){
    const v=patch.indices[j+d],a=patch.indices[j+(d+1)%3],b=patch.indices[j+(d+2)%3],link=links.get(v)??new Map();
    for(const [from,to] of [[a,b],[b,a]]){const neighbors=link.get(from)??[];neighbors.push(to);link.set(from,neighbors);}links.set(v,link);
   }
   for(const link of links.values()){
    for(const neighbors of link.values())assert.equal(neighbors.length,2);
    const reached=new Set(),queue=[link.keys().next().value];while(queue.length){const v=queue.pop();if(reached.has(v))continue;reached.add(v);queue.push(...link.get(v));}
    assert.equal(reached.size,link.size);
   }
   for(let j=0;j<patch.faces.length;j+=3)if(which)assert.ok((height[patch.faces[j]]+height[patch.faces[j+1]]+height[patch.faces[j+2]])/3>partition.threshold);

   for(const group of part.tissueDrawGroups){if(!group.count)continue;assert.ok(group.start+group.count<=part.indexCount);}
   assert.equal(part.tissueDrawGroups.reduce((sum,g)=>sum+g.count,0),part.indexCount);
  }
  assert.equal(layerUnion.size,original.size);
  const coat=atlas.concepts.find(c=>c.id===`SIW:${organ.id}:${layerName}`);assert.equal(coat.elements.length,2);
  const pair=coat.elements.map(id=>atlas.parts.find(p=>p.id===id));assert.equal(pair[0].normalSigns,pair[1].normalSigns);
  }
  assert.equal(p.patches[1].portEdges,0); assert.equal(p.patches[1].portEdges,0);assert.ok(p.patches[1].seamEdges>0);
  const ep=atlas.concepts.find(c=>c.id===`SIW:${organ.id}:epithelium`),fold=atlas.concepts.find(c=>c.id===`SIW:${organ.id}:circular-folds`);
  assert.equal(ep.elements.length,2);assert.equal(fold.elements.length,4);assert.equal(fold.elements.filter(id=>ep.elements.includes(id)).length,1);
  assert.equal(fold.partOf,`SIW:${organ.id}:mucosa`);assert.match(fold.coverageLimitation,/submucosal core/);
  const partMap=new Map(atlas.parts.map(p=>[p.id,p]));assert.deepEqual(new Set(anatomicalComponentIds(fold,[],partMap,atlas.concepts)),new Set(fold.elements));
  for(const id of fold.elements){const part=partMap.get(id),tissue=atlas.concepts.find(c=>c.id===part.conceptId);assert.equal(storedAssemblyParent(tissue,atlas.concepts,partMap).id,fold.id);}
 }
});


test('threshold saddle removes the smaller sector while retaining the stronger crest',()=>{
 const faces=Uint32Array.from(Array.from({length:8},(_,i)=>[0,i+1,(i+1)%8+1]).flat()),heights=Float32Array.from([0,1,1,0,0,.8,.8,0,0]);
 const p=partitionFoldShell(faces,heights,9);
 assert.equal(p.initialCrestFaces,2);assert.equal(p.removedSaddleFaces,1);assert.equal(p.patches[1].faces.length/3,1);
 assert.deepEqual(p.patches[1].faces,Uint32Array.from([0,1,2]));
});
