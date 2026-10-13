import type {Atlas,SceneState,Part} from './anatomy';
import {sectionPartVisible} from './viewer-state';
import type {Bounds3,Point3} from './section-plane';

export function sectionSetBounds(atlas:Atlas,state:SceneState,displayBounds?:(part:Part,index:number)=>Bounds3):Bounds3{
 const min:Point3=[Infinity,Infinity,Infinity],max:Point3=[-Infinity,-Infinity,-Infinity];
 let found=false;
 for(const [index,part] of atlas.parts.entries()){
  if(!sectionPartVisible(part,state))continue;
  const bounds=displayBounds?.(part,index)??{min:part.bounds[0],max:part.bounds[1]};
  found=true;
  for(let axis=0;axis<3;axis++){
   min[axis]=Math.min(min[axis],bounds.min[axis]);
   max[axis]=Math.max(max[axis],bounds.max[axis]);
  }
 }
 if(found&&!displayBounds&&atlas.version==='Anatomy Atlas adapted GLB local import'&&state.region==='torso'){
  // The publisher's torso section extent is anatomical. A few long vessels
  // and nerves continue below it and must not pull a 38% chest cut into the
  // abdomen when the checked structure set changes.
  min[1]=.709716765779;max[1]=1.852567015646;
 }
 return found?{min,max}:{min:[-.5,0,-.5],max:[.5,1.8,.5]};
}
