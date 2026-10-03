import {Quaternion,Vector3} from 'three';
import type {Part} from './anatomy';
import type {ExplosionStage,HierarchicalExplosionLayout} from './hierarchical-explosion';

// All exploded cluster centers share this plane. Camera roll/orbit changes
// its basis, never the relative layout. Tissue inside each cluster stays intact.
const origin=new Vector3(0,.85,0);
const cache=new WeakMap<ExplosionStage,{offsets:Float64Array;centers:Float64Array;positions:Float32Array;orientation:Quaternion|null}>();
export function viewPlanePositions(stage:ExplosionStage,parts:readonly Part[],orientation:Quaternion){
 let entry=cache.get(stage);
 if(!entry){
  const offsets=new Float64Array(stage.positions.length),centers=new Float64Array(stage.positions.length);
  for(const members of stage.clusters){
   const min=[Infinity,Infinity,Infinity],max=[-Infinity,-Infinity,-Infinity];
   for(const i of members)for(let a=0;a<3;a++){min[a]=Math.min(min[a],parts[i].bounds[0][a]);max[a]=Math.max(max[a],parts[i].bounds[1][a]);}
   for(const i of members)for(let a=0;a<3;a++){
    offsets[i*3+a]=(parts[i].bounds[0][a]+parts[i].bounds[1][a]-min[a]-max[a])/2;
    centers[i*3+a]=stage.positions[i*3+a]-offsets[i*3+a];
   }
  }
  entry={offsets,centers,positions:stage.positions.slice(),orientation:null};cache.set(stage,entry);
 }
 if(!entry.orientation?.equals(orientation)){
  const point=new Vector3();
  for(const members of stage.clusters)for(const i of members){const j=i*3;
   point.set(entry.centers[j],entry.centers[j+1]-.85,0).applyQuaternion(orientation).add(origin);
   entry.positions[j]=point.x+entry.offsets[j];entry.positions[j+1]=point.y+entry.offsets[j+1];entry.positions[j+2]=point.z+entry.offsets[j+2];
  }
  entry.orientation=orientation.clone();
 }
 return entry.positions;
}
export function viewPlaneLayout(layout:HierarchicalExplosionLayout,parts:readonly Part[],orientation:Quaternion):HierarchicalExplosionLayout{
 return {...layout,stages:layout.stages.map((stage,i)=>i===0?stage:{...stage,positions:viewPlanePositions(stage,parts,orientation)})};
}
export function viewPlanePoint(x:number,y:number,orientation:Quaternion){return new Vector3(x,y-.85,0).applyQuaternion(orientation).add(origin);}
