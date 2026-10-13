import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/organ-frame';export {Vector3,Quaternion,Matrix4} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {organPose,organBounds,Vector3,Quaternion}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
test('Organ-local thin surfaces retain separation lost in full-body float32 coordinates',()=>{
 const origin=1.2,a=Math.fround(.00002),b=Math.fround(.00002001);
 assert.notEqual(a,b);assert.equal(Math.fround(origin+a),Math.fround(origin+b));
 const part={geometryOrigin:[0,origin,0],bounds:[[0,origin+a,0],[.01,origin+b,.01]]};
 const matrix=organPose(part,new Vector3());const bounds=organBounds(part,matrix);
 assert.ok(Math.abs(bounds.min.y-(origin+a))<1e-14);assert.ok(bounds.max.y>bounds.min.y);
});
test('Bounds and points use the same local frame through offsets and anatomical rotations',()=>{
 const part={geometryOrigin:[.1,1.2,.02],bounds:[[.09,1.19,.01],[.11,1.21,.03]]};
 const offset=new Vector3(.03,-.02,.01),pivot=new Vector3(.04,1.18,.01),rotation=new Quaternion().setFromAxisAngle(new Vector3(0,0,1),Math.PI/2);
 const matrix=organPose(part,offset,rotation,pivot),bounds=organBounds(part,matrix);
 for(let x=0;x<2;x++)for(let y=0;y<2;y++)for(let z=0;z<2;z++){
  const world=new Vector3(part.bounds[x][0],part.bounds[y][1],part.bounds[z][2]);
  const expected=world.clone().add(offset).sub(pivot).applyQuaternion(rotation).add(pivot);
  const actual=world.clone().sub(new Vector3(...part.geometryOrigin)).applyMatrix4(matrix);
  assert.ok(actual.distanceTo(expected)<1e-14);assert.ok(bounds.clone().expandByScalar(1e-14).containsPoint(actual));
 }
});
test('Legacy world-packed organs preserve their existing transforms',()=>{
 const part={bounds:[[0,0,0],[1,2,3]]},offset=new Vector3(.2,.3,.4);
 const box=organBounds(part,organPose(part,offset));assert.deepEqual(box.min.toArray(),offset.toArray());assert.deepEqual(box.max.toArray(),[1.2,2.3,3.4]);
});
