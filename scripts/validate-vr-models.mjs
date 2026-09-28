import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {gunzipSync} from 'node:zlib';
import {build} from 'esbuild';
import * as T from 'three';

const bundle=await build({entryPoints:['app/vr-anatomy.ts'],bundle:true,platform:'node',format:'esm',write:false,define:{'import.meta.env.BASE_URL':'"/"'}});
const {loadVrAnatomy}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);
const originalFetch=globalThis.fetch;
globalThis.fetch=async url=>new Response(await readFile(`public${url}`));
try{
 for(const name of ['atlas-male-complete','atlas','atlas-hra-female','atlas-embryo','atlas-cell']){
  const atlas=JSON.parse(await readFile(`public/models/${name}.json`,'utf8')),model=JSON.parse(await readFile(`public/models/${name}.vr.json`,'utf8'));
  assert.equal(model.parts.length,atlas.parts.length,'every anatomical part must remain available');
  assert.ok(model.triangles<=450_000,'standalone VR triangle budget');
  assert.ok(model.maxRelativeError<=.030001,'simplification must respect its geometric error limit');
  const data=gunzipSync(await readFile(`public${model.chunk.url}`));assert.equal(data.byteLength,model.chunk.bytes);
  let triangles=0;
  for(let i=0;i<model.parts.length;i++){
   const part=model.parts[i],source=atlas.parts[i];assert.equal(part.id,source.id);assert.equal(part.sourceVertexCount,source.vertexCount);assert.equal(part.sourceIndexCount,source.indexCount);
   assert.ok(part.vertexCount>=3&&part.indexCount>0&&part.indexCount<=source.indexCount);assert.equal(part.indexCount%3,0);triangles+=part.indexCount/3;
   const indices=new Uint32Array(data.buffer,data.byteOffset+part.indices,part.indexCount),positions=new Float32Array(data.buffer,data.byteOffset+part.positions,part.vertexCount*3);
   for(const index of indices)assert.ok(index<part.vertexCount,`${part.id}: valid triangle indices`);
   for(let v=0;v<positions.length;v++){assert.ok(Number.isFinite(positions[v]));assert.ok(positions[v]>=source.bounds[0][v%3]-.001&&positions[v]<=source.bounds[1][v%3]+.001,`${part.id}: preserved placement`);}
  }
  assert.equal(model.triangles,triangles);
  const uniforms={rotationState:new T.DataTexture(new Uint8Array(T.MathUtils.ceilPowerOfTwo(atlas.parts.length)*4),T.MathUtils.ceilPowerOfTwo(atlas.parts.length),1),rotation:{value:new T.Vector4(0,0,0,1)},pivot:{value:new T.Vector3()},planes:[new T.Plane(new T.Vector3(0,1,0),10000),new T.Plane(new T.Vector3(0,1,0),10000)]};
  const vr=await loadVrAnatomy(atlas,`/models/${name}.vr.json`,uniforms,new AbortController().signal);
  assert.equal(vr.pickers.length,atlas.parts.length);assert.equal(vr.mesh.material.type,'MeshLambertMaterial');assert.equal(vr.mesh.material.transparent,false);
  const offsets=new Float32Array(atlas.parts.length*4);for(let i=0;i<atlas.parts.length;i++)offsets[i*4+3]=1;
  const selected=atlas.parts.findIndex(part=>(atlas.materials?.[part.material]?.opacity??1)>.999),state={selected:[atlas.parts[selected].id]};
  const sourceMesh=new T.Mesh();sourceMesh.matrixWorld.makeTranslation(1,2,3);const sourcePickers=[];sourcePickers[selected]=sourceMesh;
  vr.sync(state,offsets,sourcePickers,0,1);
  assert.equal(vr.mesh.geometry.drawRange.count,model.parts[selected].indexCount,'zero-opacity context must be excluded from rendering');
  assert.deepEqual(vr.pickers[selected].matrixWorld.elements,sourceMesh.matrixWorld.elements,'controller picking must use the visible part transform');
  offsets[selected*4+3]=0;vr.sync(state,offsets,sourcePickers,0,1);assert.equal(vr.mesh.geometry.drawRange.count,0,'hidden geometry must not be submitted');
  offsets[selected*4+3]=1;vr.sync(state,offsets,sourcePickers,0,1);
  const picker=vr.pickers[selected],positions=picker.geometry.getAttribute('position'),indices=picker.geometry.getIndex(),a=new T.Vector3(),b=new T.Vector3(),c=new T.Vector3();let hit=false;
  for(let i=0;i<indices.count;i+=3){a.fromBufferAttribute(positions,indices.getX(i));b.fromBufferAttribute(positions,indices.getX(i+1));c.fromBufferAttribute(positions,indices.getX(i+2));const normal=b.clone().sub(a).cross(c.clone().sub(a));if(normal.lengthSq()<1e-12)continue;normal.normalize();const center=a.clone().add(b).add(c).multiplyScalar(1/3).applyMatrix4(picker.matrixWorld);hit=new T.Raycaster(center.clone().addScaledVector(normal,.01),normal.negate()).intersectObject(picker,false).length>0;break;}
  assert.ok(hit,'the reduced mesh must support controller ray selection');
  vr.dispose();uniforms.rotationState.dispose();
  console.log(`${name}: ${atlas.parts.length} selectable parts, ${triangles.toLocaleString()} triangles; visibility and controller picking passed.`);
 }
}finally{globalThis.fetch=originalFetch;}
