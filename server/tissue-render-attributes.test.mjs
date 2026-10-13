import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync,existsSync} from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/tissue-wall-material';export * from './app/safe-surface-normals';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {tissueWallMaterial,applySafeSurfaceNormals}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));

test('Procedural angle and missing-normal fallbacks compose with the carrier shader',()=>{
 const {material}=tissueWallMaterial({color:[.7,.5,.4],pattern:'longitudinal',scale:.001,relief:.000008,courseLength:.3});
 applySafeSurfaceNormals(material);
 const shader={uniforms:{},vertexShader:'#include <beginnormal_vertex>\n#include <begin_vertex>',fragmentShader:'#include <color_fragment>\n#include <normal_fragment_begin>\n#include <normal_fragment_maps>'};
 material.onBeforeCompile(shader,{});
 assert.ok(shader.fragmentShader.includes('dot(tissueCourse.yz,tissueCourse.yz)>1e-12?atan'));
 assert.ok(shader.fragmentShader.includes('float tissueAngle(vec2 p)'));
 assert.ok(shader.fragmentShader.includes('dot(vNormal,vNormal)>1e-12?normalize(vNormal)'));
 assert.ok(shader.fragmentShader.includes('cross(dFdx(vViewPosition),dFdy(vViewPosition))'));
 assert.ok(shader.fragmentShader.includes('float height=tissueDetail*'));
 material.dispose();
});

for(const sex of ['female','male'])test(`${sex}: every indexed procedural chart has a valid circumferential coordinate`,{skip:!existsSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)},()=>{
 const atlas=JSON.parse(readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),chunks=new Map();
 let checked=0;
 for(const p of atlas.parts){
  if(p.suppressed||p.tissueChart===undefined)continue;
  if(!chunks.has(p.chunk))chunks.set(p.chunk,readFileSync(atlas.chunks[p.chunk].url.replace('/local-models/','.local-models/')));
  const b=chunks.get(p.chunk),ids=new Set();
  for(let i=0;i<p.indexCount;i++)ids.add(b.readUInt32LE(p.indices+i*4));
  for(const i of ids){
   const y=b.readInt16LE(p.tissueChart+i*6+2),z=b.readInt16LE(p.tissueChart+i*6+4);
   assert.ok(y*y+z*z>1,p.id+' has invalid chart at vertex '+i);
  }
  checked++;
 }
 assert.equal(checked,sex==='female'?43:50);
});
