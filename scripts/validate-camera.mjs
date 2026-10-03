import assert from 'node:assert/strict';
import {PerspectiveCamera,Vector3,Vector2,Raycaster,Plane} from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';
import {selectionCenter} from '../app/camera-pivot.ts';
const parts=[{id:'left',bounds:[[1,2,3],[3,4,5]]},{id:'right',bounds:[[-3,0,0],[-1,2,2]]}];
const offsets=new Float32Array([10,0,0,1,0,6,0,1]);
assert.deepEqual(selectionCenter(parts,['left'],offsets).toArray(),[12,3,4]);
assert.deepEqual(selectionCenter(parts,['left','right'],offsets).toArray(),[5,5,2.5]);
assert.equal(selectionCenter(parts,[],offsets),null);
// Exercise the installed camera controller against a synthetic event target.
const listeners=new Map();const ownerDocument={addEventListener:()=>{},removeEventListener:()=>{}};const element={style:{},ownerDocument,getRootNode:()=>ownerDocument,clientWidth:1000,clientHeight:800,addEventListener:(n,f)=>listeners.set(n,f),removeEventListener:()=>{},getBoundingClientRect:()=>({left:0,top:0,width:1000,height:800})};
const camera=new PerspectiveCamera(34,1.25,.005,100);camera.position.set(0,0,4);const controls=new OrbitControls(camera,element);controls.zoomToCursor=true;controls.update();
const pointer=new Vector2(.4,.2),ray=new Raycaster();ray.setFromCamera(pointer,camera);const point=new Vector3();ray.ray.intersectPlane(new Plane(new Vector3(0,0,1),0),point);
const before=point.clone().project(camera),oldZ=camera.position.z;
listeners.get('wheel')({clientX:700,clientY:320,deltaY:-100,preventDefault:()=>{}});
const after=point.clone().project(camera);assert.ok(camera.position.z<oldZ);assert.ok(Math.abs(before.x-after.x)<1e-6&&Math.abs(before.y-after.y)<1e-6,'Zoom must retain the cursor anchor');
const pivot=selectionCenter(parts,['left'],offsets);controls.target.copy(pivot);camera.position.copy(pivot).add(new Vector3(0,0,4));controls.update();const radius=camera.position.distanceTo(pivot);controls.autoRotate=true;controls.update();assert.ok(camera.position.x!==pivot.x);assert.ok(Math.abs(camera.position.distanceTo(pivot)-radius)<1e-6);assert.ok(pivot.clone().project(camera).length()<1.01);controls.dispose();
console.log('Cursor-anchored zoom, selected/group/exploded pivots, and orbit radius passed.');

// The small embryo must cross both world poles and turn upside down.
const {freeOrbit}=await import('../app/pointer-orbit.ts');
const freeCamera=new PerspectiveCamera(34,1.25,.005,100);
const freeTarget=new Vector3(0,.97,.17);freeCamera.position.copy(freeTarget).add(new Vector3(0,0,.5));
const freeControls=new OrbitControls(freeCamera,element);freeControls.enableDamping=false;freeControls.target.copy(freeTarget);freeControls.update();
const startOffset=freeCamera.position.clone().sub(freeTarget),startUp=freeCamera.up.clone();
for(let i=0;i<120;i++){
 freeOrbit(freeCamera,freeTarget,0,5,600);freeControls.update();
 assert.ok(Math.abs(freeCamera.position.distanceTo(freeTarget)-.5)<1e-5,'Free rotation preserves radius through the poles');
 if(i===59){assert.ok(freeCamera.up.y<-.999,'Embryo can turn upside down');assert.ok(freeCamera.position.z<freeTarget.z,'Vertical orbit reaches the opposite side');}
}
assert.ok(freeCamera.position.clone().sub(freeTarget).distanceTo(startOffset)<.001,'Full vertical turn returns to its initial position');
assert.ok(freeCamera.up.distanceTo(startUp)<.002,'Full turn preserves orientation');
freeOrbit(freeCamera,freeTarget,40,25,600);freeControls.update();
assert.ok(Math.abs(freeCamera.position.distanceTo(freeTarget)-.5)<1e-5);
freeControls.dispose();
console.log('Embryo free rotation: poles, upside-down view, complete vertical turn and diagonal drag passed.');

const {panPinch}=await import('../app/pointer-orbit.ts');
const pinchCamera=new PerspectiveCamera(34,1.25,.005,100),pinchTarget=new Vector3(0,.97,.17);
pinchCamera.position.copy(pinchTarget).add(new Vector3(0,0,.5));pinchCamera.lookAt(pinchTarget);freeOrbit(pinchCamera,pinchTarget,40,120,600);
const pinchUp=pinchCamera.up.clone(),pinchDirection=pinchCamera.getWorldDirection(new Vector3());
const beforePinch=[new Vector2(100,200),new Vector2(200,200)],afterPinch=[new Vector2(60,220),new Vector2(260,220)];
panPinch(pinchCamera,pinchTarget,beforePinch,afterPinch,600,.003,Infinity);
assert.ok(Math.abs(pinchCamera.position.distanceTo(pinchTarget)-.25)<1e-10,'Doubling finger separation halves camera distance');
assert.ok(pinchCamera.up.distanceTo(pinchUp)<1e-10,'Pinch retains freely rotated up vector');
assert.ok(pinchCamera.getWorldDirection(new Vector3()).distanceTo(pinchDirection)<1e-10,'Two-finger pan keeps orientation');
assert.ok(pinchTarget.distanceTo(new Vector3(0,.97,.17))>.001,'Moving the finger midpoint pans the view');
panPinch(pinchCamera,pinchTarget,afterPinch,beforePinch,600,.003,Infinity);
assert.ok(Math.abs(pinchCamera.position.distanceTo(pinchTarget)-.5)<1e-10,'Reverse pinch zooms back out');
console.log('Embryo two-finger pan/pinch zoom and free orientation preservation passed.');

// Off-axis surface/spine pivots keep their screen position through free orbit.
const {freeOrbitAroundPivot}=await import('../app/pointer-orbit.ts');
for(const pivot of [new Vector3(.16,.9,.12),new Vector3(0,.95,.1)]){
 const c=new PerspectiveCamera(34,1.25,.005,100),t=new Vector3(0,.97,.17);c.position.copy(t).add(new Vector3(0,0,.5));c.lookAt(t);c.updateMatrixWorld(true);
 const projected=pivot.clone().project(c),radius=c.position.distanceTo(pivot),initial=c.position.clone();
 freeOrbitAroundPivot(c,t,pivot,0,0,600);assert.ok(c.position.equals(initial),'Choosing a pivot does not jump');
 for(let i=0;i<120;i++){freeOrbitAroundPivot(c,t,pivot,0,5,600);assert.ok(Math.abs(c.position.distanceTo(pivot)-radius)<1e-9);const actual=pivot.clone().project(c);assert.ok(Math.abs(actual.x-projected.x)<1e-9&&Math.abs(actual.y-projected.y)<1e-9,'Off-axis pivot stays at its screen point');}
 assert.ok(c.position.distanceTo(initial)<1e-9,'Full pivot orbit returns to its initial pose');
}
console.log('Embryo surface/spine pivots: no initial jump, fixed projected point, constant radius and complete free turn passed.');
