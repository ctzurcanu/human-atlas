import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const compiled=await build({stdin:{contents:"export * from './app/organ-render-pose';export * from './app/organ-frame';export {Vector3,Quaternion,Box3,ShaderLib} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const shader=()=>({uniforms:{geometryOrigin:{value:new api.Vector3(99,99,99)}},vertexShader:api.ShaderLib.standard.vertexShader});
for(const sex of ['male','female']){
 const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 test(`${sex}: every displayed instance uses the part pose, independently of legacy material placement`,()=>{
  let checked=0;
  for(const p of a.parts.filter(p=>p.geometryOrigin||p.geometryInstance)){
   const s=shader();api.applyOrganRenderPose(s,p);
   assert.deepEqual(s.uniforms.geometryOrigin.value.toArray(),p.geometryOrigin??[0,0,0],p.id);
   if(p.geometryInstance){
    assert.equal(s.uniforms.geometryInstanceScale.value,p.geometryInstance.scale,p.id);
    assert.deepEqual(s.uniforms.geometryInstanceRotation.value.toArray(),p.geometryInstance.rotation,p.id);
    assert.ok(s.vertexShader.includes('transformed=instanceVector(transformed*geometryInstanceScale);'));
   }else assert.equal(s.uniforms.geometryInstanceScale,undefined,p.id);
   checked++;
  }
  assert.ok(checked>40);
 });
 test(`${sex}: rendered epicardium vertices match the audited assembled bounds`,()=>{
  const p=a.parts.find(p=>p.id===`TA98HEART:${sex}:epicardium`),s=shader();api.applyOrganRenderPose(s,p);
  const bytes=fs.readFileSync(a.chunks[p.chunk].url.replace('/local-models/','.local-models/'));
  const positions=new Float32Array(bytes.buffer,bytes.byteOffset+p.positions,p.vertexCount*3),bounds=new api.Box3();
  const q=new api.Quaternion().fromArray(s.uniforms.geometryInstanceRotation.value.toArray());
  const indices=new Uint32Array(bytes.buffer,bytes.byteOffset+p.indices,p.indexCount);
  // Tissue cells share a position bank; only this cell's indexed vertices render.
  for(const vertex of new Set(indices)){const i=vertex*3;
   const raw=new api.Vector3(positions[i],positions[i+1],positions[i+2]).multiplyScalar(p.geometryScale??1).add(new api.Vector3(...(p.sourceOffset??[0,0,0])));
   const rendered=raw.clone().multiplyScalar(s.uniforms.geometryInstanceScale.value).applyQuaternion(q).add(s.uniforms.geometryOrigin.value);
   assert.ok(rendered.distanceTo(api.organAnatomicalPoint(p,raw))<1e-12);bounds.expandByPoint(rendered);
  }
  for(const [actual,wanted] of [[bounds.min,p.bounds[0]],[bounds.max,p.bounds[1]]])assert.ok(actual.distanceTo(new api.Vector3(...wanted))<1e-8);
 });
}
test('material sharing separates different part placements and keeps identical poses reusable',()=>{
 const p={geometryOrigin:[0,1,0],geometryInstance:{rotation:[0,0,0,1],scale:.66}},same=structuredClone(p),other=structuredClone(p);other.geometryInstance.scale=1.04;
 assert.equal(api.organRenderPoseKey(p),api.organRenderPoseKey(same));assert.notEqual(api.organRenderPoseKey(p),api.organRenderPoseKey(other));
 assert.notEqual(api.organRenderPoseKey(p),api.organRenderPoseKey({geometryOrigin:[0,1.02,0],geometryInstance:p.geometryInstance}));
 const s=shader();api.applyOrganRenderPose(s);assert.deepEqual(s.uniforms.geometryOrigin.value.toArray(),[0,0,0]);assert.equal(s.uniforms.geometryInstanceScale,undefined);
});
test('scene materials bind both their shader and sharing key to the authoritative part pose',()=>{
 const scene=fs.readFileSync('app/scene.tsx','utf8');
 assert.ok(scene.includes('applyOrganRenderPose(shader,pose)'));
 assert.ok(scene.includes('pose:part'));
 assert.ok(scene.includes('+organRenderPoseKey(part)'));
 assert.ok(!scene.includes('if(source?.geometryInstance)'));
});
