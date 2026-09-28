import assert from 'node:assert/strict';
import {Box3,MathUtils,Matrix4,PerspectiveCamera,Quaternion,Vector3} from 'three';
import {applyFrame,captureFrame,frameDistance,immersiveArea,interpolateFrame,preferredAnchor,screenAnchor,setFrameOffset,usableViewArea} from '../app/view-framing.ts';
import {validCamera} from '../shared/camera-frame.mjs';

const bounds=new Box3(new Vector3(-.35,0,-.2),new Vector3(.35,1.8,.2)),target=bounds.getCenter(new Vector3());
const desktop=usableViewArea([{kind:'panel',left:0,right:.2,top:.08,bottom:.86},{kind:'panel',left:.76,right:1,top:.08,bottom:1},{kind:'top',left:0,right:1,top:0,bottom:.07},{kind:'bottom',left:.21,right:.75,top:.85,bottom:1}]);
const camera=new PerspectiveCamera(34,16/9,.0001,100),direction=new Vector3(.35,.06,1).normalize();
camera.position.copy(target).addScaledVector(direction,frameDistance(bounds,direction,camera.up,camera.aspect,desktop,camera.fov)*1.08);camera.lookAt(target);setFrameOffset(camera,desktop);camera.updateMatrixWorld(true);
const saved=captureFrame(camera,target,bounds,desktop);assert.ok(validCamera(saved));assert.equal(saved.length,12);
const near=(a,b)=>assert.ok(Math.abs(a-b)<1e-7,`${a} != ${b}`);
const same=(a,b)=>a.forEach((value,i)=>near(value,b[i]));
for(const [aspect,area] of [[9/19.5,usableViewArea([{kind:'top',left:0,right:1,top:0,bottom:.11},{kind:'bottom',left:0,right:1,top:.7,bottom:1}])],[19.5/9,usableViewArea([{kind:'panel',left:.7,right:1,top:.05,bottom:.95},{kind:'bottom',left:0,right:.7,top:.7,bottom:1}])],[1,immersiveArea]]){
 const local=new PerspectiveCamera(34,aspect,.0001,100),localTarget=new Vector3();applyFrame(local,localTarget,saved,bounds,area);
 same(captureFrame(local,localTarget,bounds,area),saved);
 const projected=target.clone().project(local),anchor=screenAnchor(area,saved.slice(6,8));near((projected.x+1)/2,anchor[0]);near((1-projected.y)/2,anchor[1]);assert.ok(anchor[1]<.5);
 for(let corner=0;corner<8;corner++){
  const point=new Vector3(bounds[(corner&1)?'max':'min'].x,bounds[(corner&2)?'max':'min'].y,bounds[(corner&4)?'max':'min'].z).project(local),x=(point.x+1)/2,y=(1-point.y)/2;
  assert.ok(x>=area.left&&x<=area.right&&y>=area.top&&y<=area.bottom,'Whole-model framing stays inside the local usable area.');
 }
}
const scaled=new Box3(bounds.min.clone().multiplyScalar(10).addScalar(25),bounds.max.clone().multiplyScalar(10).addScalar(25)),scaledCamera=new PerspectiveCamera(34,16/9,.0001,1000),scaledTarget=new Vector3();applyFrame(scaledCamera,scaledTarget,saved,scaled,desktop);same(captureFrame(scaledCamera,scaledTarget,scaled,desktop),saved);
const moreTools=usableViewArea([{kind:'top',left:0,right:1,top:0,bottom:.07},{kind:'bottom',left:0,right:1,top:.64,bottom:1}]);assert.ok(screenAnchor(moreTools,preferredAnchor(moreTools))[1]<screenAnchor(desktop,preferredAnchor(desktop))[1]);
const axial=new PerspectiveCamera(34,1,.0001,100);axial.up.set(0,0,-1);axial.position.copy(target).add(new Vector3(0,3,0));axial.lookAt(target);setFrameOffset(axial,desktop);const axialFrame=captureFrame(axial,target,bounds,desktop),restored=new PerspectiveCamera(34,1,.0001,100),restoredTarget=new Vector3();applyFrame(restored,restoredTarget,axialFrame,bounds,desktop);same(restored.up.toArray(),axial.up.toArray());
interpolateFrame(restored,restoredTarget,saved,axialFrame,bounds,desktop,1);same(captureFrame(restored,restoredTarget,bounds,desktop),axialFrame);
assert.ok(validCamera([1,1,3,0,.8,0,0,0]));for(const invalid of [[...saved.slice(0,8),1,...saved.slice(9)],[.5,.5,.5,...saved.slice(3)],[...saved.slice(0,6),-1,...saved.slice(7)]])assert.equal(validCamera(invalid),false);
// The XR presentation transform places the starting head at the requested view,
// then preserves an independently tracked head movement relative to that view.
const initialHead=new Matrix4().compose(new Vector3(.1,1.6,.2),new Quaternion().setFromAxisAngle(new Vector3(0,1,0),.3),new Vector3(1,1,1));
const presentation=new Matrix4().compose(new Vector3(0,.9,2),new Quaternion().setFromAxisAngle(new Vector3(0,1,0),-.2),new Vector3(1,1,1)),rig=presentation.clone().multiply(initialHead.clone().invert());
same(rig.clone().multiply(initialHead).elements,presentation.elements);
const movement=new Matrix4().makeTranslation(.1,0,0);same(rig.clone().multiply(initialHead.clone().multiply(movement)).elements,presentation.clone().multiply(movement).elements);
console.log('Normalized framing round-trips across portrait, landscape, XR and model scale; tools move the center upward; axial views and tracked head movement are preserved.');
