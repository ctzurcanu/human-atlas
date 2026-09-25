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
