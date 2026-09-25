import type {Atlas,SceneState} from './anatomy';
import {sectionPartVisible} from './viewer-state';
import type {Bounds3,Point3} from './section-plane';

export function sectionSetBounds(atlas:Atlas,state:SceneState):Bounds3{
 const min:Point3=[Infinity,Infinity,Infinity],max:Point3=[-Infinity,-Infinity,-Infinity];
 let found=false;
 for(const part of atlas.parts){
  if(!sectionPartVisible(part,state))continue;
  found=true;
  for(let axis=0;axis<3;axis++){
   min[axis]=Math.min(min[axis],part.bounds[0][axis]);
   max[axis]=Math.max(max[axis],part.bounds[1][axis]);
  }
 }
 return found?{min,max}:{min:[-.5,0,-.5],max:[.5,1.8,.5]};
}
