import test from 'node:test';import assert from 'node:assert/strict';import {build} from 'esbuild';import fs from 'node:fs';
const b=await build({stdin:{contents:"export * from './app/tissue-interface-visibility';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});const {tissueNeighborOccludes}=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
test('Selected opaque coat keeps both walls against invisible or faint context; assembled peers suppress shared faces',()=>{
 assert.equal(tissueNeighborOccludes(1,0),false);assert.equal(tissueNeighborOccludes(1,.18),false);assert.equal(tissueNeighborOccludes(1,1),true);assert.equal(tissueNeighborOccludes(.18,1),true);assert.equal(tissueNeighborOccludes(.18,.18),true);
 for(const x of [NaN,Infinity,-1]){assert.equal(tissueNeighborOccludes(1,x),false);assert.equal(tissueNeighborOccludes(x,1),false);}
});
for(const sex of ['female','male'])test(sex+': isolated longitudinal union retains outer/inner surfaces and suppresses only same-coat seams',()=>{
 const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),ids=new Set(a.concepts.find(c=>c.id==='CW:colon:longitudinal-muscle').elements),parts=new Map(a.parts.map(p=>[p.id,p]));
 for(const id of ids){const p=parts.get(id),k=2,n=p.vertexCount/9,chunk=a.chunks[p.chunk],raw=fs.readFileSync('.'+chunk.url.replace('/local-models/','/.local-models/'));const kept=p.tissueDrawGroups.filter(g=>!g.neighbor||!tissueNeighborOccludes(1,ids.has(g.neighbor)?1:0));let outer=0,inner=0;
  for(const g of kept){const indices=new Uint32Array(raw.buffer,raw.byteOffset+p.indices+g.start*4,g.count);for(let j=0;j<indices.length;j+=3){const face=Array.from(indices.subarray(j,j+3));if(face.every(v=>v>=k*n&&v<(k+1)*n))outer++;if(face.every(v=>v>=(k+1)*n&&v<(k+2)*n))inner++;}}
  assert.ok(outer>0&&inner>0,id+' retains both envelope surfaces');
 }
});
