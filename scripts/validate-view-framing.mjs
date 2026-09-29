import assert from 'node:assert/strict';
import {Box3,MathUtils,Matrix4,PerspectiveCamera,Quaternion,Vector3} from 'three';
import {applyFrame,captureCameraPose,cameraPoseForFrame,captureFrame,frameDistance,immersiveArea,interpolateCameraPose,interpolateFrame,preferredAnchor,screenAnchor,setFrameOffset,usableViewArea,zoomInToFill,zoomOutToFit} from '../app/view-framing.ts';
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

// The next exploded mark must fit even with a panned target, depth, a tilted
// view, camera zoom, and narrow space between panels. Never zoom back in.
const exploded=new Box3(new Vector3(-2,-1,-.7),new Vector3(3,3,.9));
for(const [aspect,area] of [[16/9,desktop],[9/19.5,moreTools],[1,immersiveArea]]){
 const local=new PerspectiveCamera(34,aspect,.01,100),pivot=new Vector3(.6,.8,.1),back=new Vector3(.35,.12,1).normalize();
 local.zoom=1.3;local.position.copy(pivot).addScaledVector(back,.5);local.lookAt(pivot);setFrameOffset(local,area);local.updateMatrixWorld(true);
 const orientation=local.quaternion.toArray(),offset={...local.view},start=local.position.clone(),distance=local.position.distanceTo(pivot);
 assert.equal(zoomOutToFit(local,exploded,area),true);
 assert.ok(local.position.distanceTo(pivot)>distance);
 same(local.quaternion.toArray(),orientation);assert.deepEqual(local.view,offset);
 const delta=local.position.clone().sub(start).normalize();same(delta.toArray(),back.toArray());
 for(let corner=0;corner<8;corner++){
  const point=new Vector3(exploded[(corner&1)?'max':'min'].x,exploded[(corner&2)?'max':'min'].y,exploded[(corner&4)?'max':'min'].z).project(local),x=(point.x+1)/2,y=(1-point.y)/2;
  assert.ok(x>=area.left&&x<=area.right&&y>=area.top&&y<=area.bottom&&point.z>=-1&&point.z<=1,'Every corner of the upcoming explosion fits.');
 }
 const fitted=local.position.toArray();
 assert.equal(zoomOutToFit(local,exploded,area),false);same(local.position.toArray(),fitted);
 assert.equal(zoomOutToFit(local,bounds,area),false);same(local.position.toArray(),fitted);
}
const coveredCenter=new PerspectiveCamera(34,16/9,.01,100);coveredCenter.position.set(0,.9,1);coveredCenter.lookAt(target);
assert.equal(zoomOutToFit(coveredCenter,exploded,{left:.65,right:.98,top:.1,bottom:.8}),true);
const uncovered=target.clone().project(coveredCenter);assert.ok((uncovered.x+1)/2>.65,'An explicit fit can uncover a projection center hidden by a newly opened panel.');

const projectedArea=(camera,box)=>{
 let left=Infinity,right=-Infinity,top=Infinity,bottom=-Infinity;
 for(let corner=0;corner<8;corner++){
  const point=new Vector3(box[(corner&1)?'max':'min'].x,box[(corner&2)?'max':'min'].y,box[(corner&4)?'max':'min'].z).project(camera),x=(point.x+1)/2,y=(1-point.y)/2;
  left=Math.min(left,x);right=Math.max(right,x);top=Math.min(top,y);bottom=Math.max(bottom,y);
 }
 return {left,right,top,bottom};
};
// During implosion the content grows on screen as its world bounds shrink.
// Matching the content to each viewport also covers portrait and XR ratios.
for(const [aspect,area] of [[16/9,desktop],[9/19.5,moreTools],[1,immersiveArea]]){
 const local=new PerspectiveCamera(34,aspect,.01,100),pivot=new Vector3(.1,.9,0),halfWidth=aspect*(area.right-area.left)/(area.bottom-area.top);
 const box=new Box3(pivot.clone().sub(new Vector3(halfWidth,1,.08)),pivot.clone().add(new Vector3(halfWidth,1,.08)));
 local.position.copy(pivot).add(new Vector3(0,0,25));local.lookAt(pivot);setFrameOffset(local,area);local.updateMatrixWorld(true);
 const orientation=local.quaternion.toArray(),offset={...local.view};let previousDistance=local.position.distanceTo(pivot);
 for(const scale of [1,.65,.3]){
  const shrunk=new Box3(box.min.clone().sub(pivot).multiplyScalar(scale).add(pivot),box.max.clone().sub(pivot).multiplyScalar(scale).add(pivot));
  assert.equal(zoomInToFill(local,shrunk,area),true);
  const visible=projectedArea(local,shrunk);
  assert.ok((visible.right-visible.left)/(area.right-area.left)>.5&&(visible.bottom-visible.top)/(area.bottom-area.top)>.5,'Imploding content fills more than half of both usable dimensions');
  assert.ok(visible.left>=area.left&&visible.right<=area.right&&visible.top>=area.top&&visible.bottom<=area.bottom,'Zoom-in retains the entire shrinking view');
  const distance=local.position.distanceTo(pivot);assert.ok(distance<previousDistance);previousDistance=distance;
  same(local.quaternion.toArray(),orientation);assert.deepEqual(local.view,offset);
  assert.equal(zoomInToFill(local,shrunk,area),false,'Stop zooming when the minimum size is reached');
 }
}
const slender=new Box3(new Vector3(-.15,-1,-.1),new Vector3(.15,1,.1)),slenderCamera=new PerspectiveCamera(34,16/9,.01,100),slenderArea=usableViewArea([]);
slenderCamera.position.set(0,0,25);slenderCamera.lookAt(new Vector3());setFrameOffset(slenderCamera,slenderArea);slenderCamera.updateMatrixWorld(true);
assert.equal(zoomInToFill(slenderCamera,slender,slenderArea),true);
const slenderScreen=projectedArea(slenderCamera,slender);
assert.ok(slenderScreen.top>=slenderArea.top&&slenderScreen.bottom<=slenderArea.bottom,'A tall body stays fully visible even when both half-size goals cannot coexist');
near(slenderScreen.top,slenderArea.top+(slenderArea.bottom-slenderArea.top)*.02);
assert.equal(zoomInToFill(slenderCamera,slender,slenderArea),false,'Do not push a constrained dimension beyond the viewport');
const limited=new PerspectiveCamera(34,1,.01,100);limited.position.set(0,0,20);limited.lookAt(new Vector3());limited.updateMatrixWorld(true);
zoomInToFill(limited,slender,slenderArea,2);near(limited.position.z,18);
// A slide can change from the whole body to a distant, small selected organ.
// Keep the actual source camera and decode the destination only once.
const sourcePose=captureCameraPose(camera,target),organ=new Box3(new Vector3(2,4,1),new Vector3(2.1,4.2,1.1));
const destination=cameraPoseForFrame(saved,organ,camera,desktop),animated=new PerspectiveCamera(34,16/9,.0001,100),animatedTarget=new Vector3();
interpolateCameraPose(animated,animatedTarget,sourcePose,destination,0);
same(animated.position.toArray(),sourcePose.position.toArray());same(animatedTarget.toArray(),sourcePose.target.toArray());
interpolateCameraPose(animated,animatedTarget,sourcePose,destination,.5);const midpoint=animated.position.toArray();
organ.translate(new Vector3(100,100,100));interpolateCameraPose(animated,animatedTarget,sourcePose,destination,.5);same(animated.position.toArray(),midpoint);
interpolateCameraPose(animated,animatedTarget,sourcePose,destination,1);same(animated.position.toArray(),destination.position.toArray());same(animatedTarget.toArray(),destination.target.toArray());
const inverted={...destination,up:destination.up.clone().negate()};
for(let i=0;i<=100;i++){interpolateCameraPose(animated,animatedTarget,destination,inverted,i/100);assert.ok(animated.position.toArray().every(Number.isFinite));assert.ok(animated.up.length()>.99,'Opposing camera up vectors never collapse during a transition');}
console.log('Normalized framing preserves pan and orientation; explosion fits ahead, and implosion fills half of both dimensions or the largest fully visible size.');
