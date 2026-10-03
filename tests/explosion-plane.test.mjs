import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
import {Quaternion,Euler,Vector3} from 'three';
const compiled=await build({entryPoints:['app/explosion-plane.ts','app/hierarchical-explosion.ts','app/anatomy-hierarchy.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/explosion-plane-test'});
const load=async name=>import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles.find(f=>f.path.endsWith(name+'.js')).contents).toString('base64'));
const {viewPlaneLayout}=await load('explosion-plane');
const {createHierarchicalExplosionLayout}=await load('hierarchical-explosion');
const {hierarchyEntries}=await load('anatomy-hierarchy');
const {readFile}=await import('node:fs/promises');
for(const name of ['atlas-male-complete','atlas-hra-female','atlas-embryo','atlas-cell','embryo-3month/atlas'])test(`${name}: every exploded cluster stays on one POV plane through orbit and roll`,async()=>{
 const atlas=JSON.parse(await readFile(`public/models/${name}.json`));
 const layout=createHierarchicalExplosionLayout(atlas,hierarchyEntries(atlas),new Set(atlas.parts.map(p=>p.id)),'systems',1.7);
 for(const rotation of [[0,0,0],[.6,1.4,.9],[Math.PI,0,Math.PI/2],[-1.1,-2.2,2.5]]){
  const q=new Quaternion().setFromEuler(new Euler(...rotation)),normal=new Vector3(0,0,1).applyQuaternion(q),right=new Vector3(1,0,0).applyQuaternion(q),up=new Vector3(0,1,0).applyQuaternion(q);
  const mapped=viewPlaneLayout(layout,atlas.parts,q);
  assert.equal(mapped.stages[0],layout.stages[0],'Unexploded assembly must remain unchanged');
  for(let level=1;level<mapped.stages.length;level++)for(const members of mapped.stages[level].clusters){
   const min=[Infinity,Infinity,Infinity],max=[-Infinity,-Infinity,-Infinity];
   for(const i of members)for(let a=0;a<3;a++){min[a]=Math.min(min[a],atlas.parts[i].bounds[0][a]);max[a]=Math.max(max[a],atlas.parts[i].bounds[1][a]);}
   let group;
   for(const i of members){
    const offset=new Vector3(...[0,1,2].map(a=>(atlas.parts[i].bounds[0][a]+atlas.parts[i].bounds[1][a]-min[a]-max[a])/2));
    const center=new Vector3().fromArray(mapped.stages[level].positions,i*3).sub(offset);
    assert.ok(Math.abs(center.clone().sub(new Vector3(0,.85,0)).dot(normal))<2e-5,'Common plane must stay perpendicular to POV');
    if(group)assert.ok(center.distanceTo(group)<2e-5,'Cluster anatomy must retain its relative positions');else group=center;
    const original=new Vector3().fromArray(layout.stages[level].positions,i*3).sub(offset).sub(new Vector3(0,.85,0));
    const placed=center.clone().sub(new Vector3(0,.85,0));
    assert.ok(Math.abs(placed.dot(right)-original.x)<2e-5&&Math.abs(placed.dot(up)-original.y)<2e-5,'Screen packing must remain stable');
   }
  }
 }
});
