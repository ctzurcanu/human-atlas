import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {build} from 'esbuild';
import {Box3,Vector3,PerspectiveCamera} from 'three';
const compiled=await build({stdin:{contents:"export * from './app/explosion-plane';export * from './app/hierarchical-explosion';export * from './app/anatomy-hierarchy';export * from './app/view-framing';export * from './app/explosion-view';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {viewPlaneLayout,createHierarchicalExplosionLayout,explosionBoundsAt,hierarchyEntries,setFrameOffset,captureExplosionView,fitExplosionView}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
for(const sex of ['female','male'])test(`${sex}: ordinary layered atlas explosion preserves its POV plane, viewport and exact return`,async()=>{
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),ids=new Set(atlas.parts.filter(p=>!p.suppressed).map(p=>p.id));
 for(const aspect of [.6,1.78]){
  const area={left:.23,right:.75,top:.08,bottom:.9},camera=new PerspectiveCamera(34,aspect,.00001,100),target=new Vector3(.06,.8,.11);
  camera.position.copy(target).add(new Vector3(.3,.2,.5));camera.up.set(.15,1,.2).normalize();camera.lookAt(target);setFrameOffset(camera,area);camera.updateMatrixWorld(true);
  const original=captureExplosionView(camera,target),layout=createHierarchicalExplosionLayout(atlas,hierarchyEntries(atlas),ids,'systems',aspect),mapped=viewPlaneLayout(layout,atlas.parts,camera.quaternion),normal=new Vector3(0,0,1).applyQuaternion(camera.quaternion);
  assert.equal(mapped.stages[0],layout.stages[0]);
  for(const stage of mapped.stages.slice(1))for(const members of stage.clusters){const box=new Box3();for(const i of members)box.union(new Box3(new Vector3(...atlas.parts[i].bounds[0]),new Vector3(...atlas.parts[i].bounds[1])));const groupCenter=box.getCenter(new Vector3());let previous;
   for(const i of members){const center=new Box3(new Vector3(...atlas.parts[i].bounds[0]),new Vector3(...atlas.parts[i].bounds[1])).getCenter(new Vector3()),placed=new Vector3().fromArray(stage.positions,i*3).sub(center.sub(groupCenter));assert.ok(Math.abs(placed.clone().sub(new Vector3(0,.85,0)).dot(normal))<2e-5);if(previous)assert.ok(placed.distanceTo(previous)<2e-5);previous=placed;}
  }
  for(const direction of [1,-1])for(let frame=0;frame<=16;frame++){const amount=direction===1?frame/16:1-frame/16,b=explosionBoundsAt(mapped,atlas.parts,ids,amount),box=new Box3(new Vector3(...b.min),new Vector3(...b.max));fitExplosionView(camera,target,original,box,area,amount);
   if(amount===0){assert.deepEqual(captureExplosionView(camera,target),original);continue;}
   for(let i=0;i<8;i++){const p=new Vector3(box[i&1?'max':'min'].x,box[i&2?'max':'min'].y,box[i&4?'max':'min'].z).project(camera),x=(p.x+1)/2,y=(1-p.y)/2;assert.ok(x>=area.left-1e-8&&x<=area.right+1e-8&&y>=area.top-1e-8&&y<=area.bottom+1e-8&&p.z>=-1&&p.z<=1);}
  }
 }
});
