import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';

async function load(entryPoint){
 const result=await build({entryPoints:[entryPoint],bundle:true,platform:'node',format:'esm',write:false});
 return import(`data:text/javascript;base64,${Buffer.from(result.outputFiles[0].contents).toString('base64')}`);
}
const {hierarchyEntries}=await load('app/anatomy-hierarchy.ts');
const {createHierarchicalExplosionLayout}=await load('app/hierarchical-explosion.ts');

function verifyPacking(atlas,stage){
 const rectangles=atlas.parts.map((part,i)=>{
  const x=stage.positions[i*3],y=stage.positions[i*3+1];
  return {minX:x-(part.bounds[1][0]-part.bounds[0][0])/2,maxX:x+(part.bounds[1][0]-part.bounds[0][0])/2,minY:y-(part.bounds[1][1]-part.bounds[0][1])/2,maxY:y+(part.bounds[1][1]-part.bounds[0][1])/2};
 }).sort((a,b)=>a.minX-b.minX);
 let active=[];
 for(const rect of rectangles){
  active=active.filter(other=>other.maxX>rect.minX+1e-6);
  for(const other of active)assert(other.maxY<=rect.minY+1e-6||rect.maxY<=other.minY+1e-6,'Atomic pieces overlap at full explosion');
  active.push(rect);
 }
}

for(const name of ['atlas-male-complete','atlas-hra-female','atlas-embryo','atlas-cell']){
 const atlas=JSON.parse(await readFile(new URL(`../public/models/${name}.json`,import.meta.url)));
 const entries=hierarchyEntries(atlas),visible=new Set(atlas.parts.map(part=>part.id));
 const modes=name==='atlas-cell'?['systems','regions']:['systems','regions','depth'];
 const layouts=Object.fromEntries(modes.map(mode=>[mode,createHierarchicalExplosionLayout(atlas,entries,visible,mode,1.7)]));
 for(const [mode,layout] of Object.entries(layouts)){
  assert.equal(layout.stages.length,layout.steps+1);
  assert(layout.steps>=2);
  for(const stage of layout.stages){
   assert(stage.positions.every(Number.isFinite));assert(stage.width>0&&stage.height>0);
   assert.equal(stage.clusters.flat().length,atlas.parts.length,'Every visible piece must belong to one rotation cluster');
   assert.equal(new Set(stage.clusters.flat()).size,atlas.parts.length,'Rotation clusters must not overlap');
   stage.clusters.forEach((members,cluster)=>members.forEach(i=>assert.equal(stage.clusterIds[i],cluster)));
  }
  assert(layout.stages.at(-1).clusters.every(members=>members.length===1),'The final explosion must rotate atomic pieces independently');
  verifyPacking(atlas,layout.stages.at(-1));
  console.log(`${name} ${mode}: ${layout.steps} hierarchy levels; all ${atlas.parts.length} final pieces packed without overlap`);
 }
 if(name==='atlas-male-complete'){
  const bone=atlas.parts.findIndex(part=>part.name.startsWith('Femur (left) · Bone'));
  const otherBone=atlas.parts.findIndex(part=>part.name.startsWith('Humerus (left) · Bone'));
  const muscle=atlas.parts.findIndex(part=>part.name.startsWith('Trapezius muscle (left)'));
  const position=(mode,i)=>layouts[mode].stages[1].positions.slice(i*3,i*3+3);
  const displacement=(mode,i)=>[...position(mode,i)].map((value,axis)=>value-layouts[mode].stages[0].positions[i*3+axis]);
  assert(displacement('systems',bone).every((value,axis)=>Math.abs(value-displacement('systems',otherBone)[axis])<1e-6),'Skeleton must move as one group at the first level');
  assert.equal(layouts.systems.stages[1].clusterIds[bone],layouts.systems.stages[1].clusterIds[otherBone],'Selecting either bone must rotate the Skeleton cluster');
  assert.notEqual(layouts.systems.stages[1].clusterIds[bone],layouts.systems.stages[1].clusterIds[muscle],'Separate systems must stay fixed when another cluster rotates');
  const muscleClusterSize=stage=>stage.clusters[stage.clusterIds[muscle]].length;
  assert(muscleClusterSize(layouts.systems.stages[1])>muscleClusterSize(layouts.systems.stages[2]),'A deeper tick must rotate a smaller group around the selected muscle');
  assert.equal(muscleClusterSize(layouts.systems.stages.at(-1)),1,'The selected muscle must rotate alone at the last tick');
  assert(layouts.systems.stages[1].groups.filter(group=>group.depth===1).every(group=>group.z===0),'Exploded groups must share a front-facing plane');
  assert.notDeepEqual([...position('systems',bone)],[...position('regions',bone)]);
  assert.notDeepEqual([...position('systems',muscle)],[...position('depth',muscle)]);
 }
}
