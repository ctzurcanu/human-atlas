import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {existsSync} from 'node:fs';
import ts from 'typescript';
import * as T from 'three';
const source=await readFile(new URL('../app/model-visit.ts',import.meta.url),'utf8');
const js=ts.transpileModule(source,{compilerOptions:{module:ts.ModuleKind.ESNext}}).outputText.replace("from 'three'",`from '${new URL('../node_modules/three/build/three.module.js',import.meta.url).href}'`);
const {ModelVisit,VisitJourney,initialVisit,visitCurve,dorsalFrames}=await import('data:text/javascript;base64,'+Buffer.from(js).toString('base64'));
const bounds=new T.Box3(new T.Vector3(-.3,0,-.2),new T.Vector3(.3,1.8,.2));
test('Paused visit leaves camera untouched, including opening and path changes',()=>{
 const camera=new T.PerspectiveCamera(),target=new T.Vector3(.1,.2,.3);camera.position.set(1,2,3);const before=camera.toJSON();
 new ModelVisit().update(camera,target,bounds,{...initialVisit,path:'back',restart:2},1);
 assert.deepEqual(camera.toJSON(),before);assert.deepEqual(target.toArray(),[.1,.2,.3]);
});
test('Forward and reverse trace the same path; look changes orientation without moving path position',()=>{
 const a=new ModelVisit(),b=new ModelVisit(),ca=new T.PerspectiveCamera(),cb=new T.PerspectiveCamera(),ta=new T.Vector3(),tb=new T.Vector3();
 a.update(ca,ta,bounds,{...initialVisit,path:'front',playing:true},6);
 b.update(cb,tb,bounds,{...initialVisit,path:'front',playing:true,direction:-1},24);
 assert.ok(ca.position.distanceTo(cb.position)<1e-8);const position=ca.position.clone(),q=ca.quaternion.clone();a.look(60,-20);a.update(ca,ta,bounds,{...initialVisit,path:'front',playing:true},0);
 assert.ok(position.distanceTo(ca.position)<1e-12);assert.ok(q.angleTo(ca.quaternion)>.1);
});
test('Speed scales progress; open paths finish and closed paths wrap smoothly',()=>{
 const v=new ModelVisit(),camera=new T.PerspectiveCamera(),target=new T.Vector3();
 assert.equal(v.update(camera,target,bounds,{...initialVisit,path:'front',playing:true,speed:2},16),true);assert.equal(v.progress,1);
 const orbit=visitCurve(bounds,'orbit');assert.ok(orbit.getPointAt(0).distanceTo(orbit.getPointAt(1))<1e-10);
 assert.equal(v.update(camera,target,bounds,{...initialVisit,playing:true,loop:true,restart:1},31),false);assert.ok(v.progress<.04);
});
test('Smoothed alimentary traversal retains endpoints and reverses exactly',()=>{
 const points=[[0,1,0],[0,1,.2],[.1,1,.2],[.1,.8,.2]],curve=visitCurve(bounds,'alimentary',points);
 assert.ok(curve.getPointAt(0).distanceTo(new T.Vector3(...points[0]))<1e-10);assert.ok(curve.getPointAt(1).distanceTo(new T.Vector3(...points.at(-1)))<1e-10);
 assert.ok(curve.getPointAt(0).distanceTo(new T.Vector3(...points[0]))<1e-8);
 const a=new ModelVisit(),b=new ModelVisit(),ca=new T.PerspectiveCamera(),cb=new T.PerspectiveCamera(),ta=new T.Vector3(),tb=new T.Vector3();
 a.update(ca,ta,bounds,{...initialVisit,path:'alimentary',playing:true,lookAlong:true},60,points);
 b.update(cb,tb,bounds,{...initialVisit,path:'alimentary',playing:true,direction:-1,lookAlong:true},60,points);
 assert.ok(ca.position.distanceTo(cb.position)<1e-8);assert.ok(ca.getWorldDirection(new T.Vector3()).dot(cb.getWorldDirection(new T.Vector3()))<-.999);
});
test('Organ seeks move a paused camera once and progress remains frozen afterward',()=>{
 const points=[[0,1,0],[0,1,.2],[0,1,.2],[.1,1,.2],[.1,.8,.2]],v=new ModelVisit(),camera=new T.PerspectiveCamera(),target=new T.Vector3();
 const curve=visitCurve(bounds,'alimentary',points);assert.ok(curve.waypointProgress[1]>0);assert.equal(curve.waypointProgress[1],curve.waypointProgress[2]);
 const state={...initialVisit,path:'alimentary',lookAlong:true,seek:.5,seekId:1};v.update(camera,target,bounds,state,.016,points);assert.equal(v.progress,.5);const pose=camera.toJSON();v.update(camera,target,bounds,state,.016,points);assert.deepEqual(camera.toJSON(),pose);
});
test('Curved path has continuous tangents at each landmark',()=>{
 const curve=visitCurve(bounds,'alimentary',[[0,1,0],[0,1,.2],[.1,1,.2],[.1,.8,.2]]);
 for(const t of [1/3,2/3])assert.ok(curve.getTangent(t-1e-6).dot(curve.getTangent(t+1e-6))>.99999);
});
test('Dorsal frame stays anterior-up through horizontal bends without roll flips',()=>{
 const curve=visitCurve(bounds,'alimentary',[[0,1,0],[.1,.9,0],[0,.8,0],[-.1,.7,0],[0,.6,0]]),frames=dorsalFrames(curve);
 for(let i=100;i<frames.length;i++){assert.ok(new T.Vector3(0,1,0).applyQuaternion(frames[i]).z>.99);assert.ok(frames[i-1].angleTo(frames[i])<.08);}
});
test('Loop repeats an open path and terminal alimentary view faces the anal outlet',()=>{
 const camera=new T.PerspectiveCamera(),target=new T.Vector3(),v=new ModelVisit(),points=[[0,1,.6],[0,1,0],[0,.2,0],[0,.2,.6]];
 assert.equal(v.update(camera,target,bounds,{...initialVisit,path:'front',playing:true,loop:true},31),false);assert.ok(v.progress<.04);
 v.update(camera,target,bounds,{...initialVisit,path:'alimentary',restart:1,loop:true,lookAlong:true,seek:1,seekId:1},0,points);
 assert.ok(camera.getWorldDirection(new T.Vector3()).dot(new T.Vector3(0,0,-1))>.999);
});
test('Alimentary guides approach 60 cm anterior to the mouth and end between anus and toes',{skip:!['male','female'].every(sex=>existsSync(new URL(`../.local-models/ta98-review/visits/${sex}.json`,import.meta.url)))},async()=>{
 for(const sex of ['male','female']){
  const route=JSON.parse(await readFile(new URL(`../.local-models/ta98-review/visits/${sex}.json`,import.meta.url),'utf8'));
  const points=route.points,entry=route.segments[1].points[0],outlet=route.segments.at(-2).points.at(-1);
  assert.ok(Math.abs(points[0][2]-entry[2]-.6)<1e-7);assert.ok(Math.abs(points[0][1]-entry[1])<1e-7);assert.ok(points.at(-1)[1]<outlet[1]);assert.ok(points.at(-1)[2]>outlet[2]);
  const curve=visitCurve(bounds,'alimentary',points);assert.ok(curve.getPointAt(0).distanceTo(new T.Vector3(...points[0]))<1e-8);assert.ok(curve.getPointAt(1).distanceTo(new T.Vector3(...points.at(-1)))<1e-8);
  const boundaries=route.segments.map(s=>curve.waypointProgress[s.startIndex]);assert.ok(boundaries.every((n,i)=>Number.isFinite(n)&&(!i||n>=boundaries[i-1])));
 }
});

test('Playback never hands the camera back to OrbitControls or runs a covering scan',async()=>{
 const scene=await readFile(new URL('../app/scene.tsx',import.meta.url),'utf8');
 assert.ok(scene.includes("if(!presenting&&!visitRef.current?.playing&&(travel.phase!=='done'||controls.autoRotate))controls.update()"));
 assert.ok(scene.includes('if(!presenting&&ready&&!visitRef.current?.playing){const key='));
});
test('Actual male and female guides have bounded frame-to-frame camera turning', {skip:!['male','female'].every(sex=>existsSync(new URL(`../.local-models/ta98-review/visits/${sex}.json`,import.meta.url)))},async()=>{
 for(const sex of ['male','female']){const route=JSON.parse(await readFile(new URL(`../.local-models/ta98-review/visits/${sex}.json`,import.meta.url),'utf8')),v=new ModelVisit(),c=new T.PerspectiveCamera(),target=new T.Vector3(),state={...initialVisit,path:'alimentary',playing:true,lookAlong:true};
 v.update(c,target,bounds,state,0,route.points);
 for(let i=0;i<2400;i++){const q=c.quaternion.clone();v.update(c,target,bounds,state,.05,route.points);assert.ok(q.angleTo(c.quaternion)<=.062501);assert.ok(c.position.toArray().every(Number.isFinite));}
 }
});

test('Unchecked Loop ends every route, including Around the model',()=>{
 for(const path of ['orbit','front','back','alimentary'])for(const direction of [1,-1]){
 const camera=new T.PerspectiveCamera(),target=new T.Vector3(),v=new ModelVisit(),state={...initialVisit,path,direction,playing:true};
 assert.equal(v.update(camera,target,bounds,state,121,[[0,1,0],[0,.2,.1]]),true);assert.equal(v.progress,direction===1?1:0);
 }
});
test('Visit glides preserve their first frame, freeze on pause and return the exact requesting pose',()=>{
 for(const direction of [1,-1]){
 const camera=new T.PerspectiveCamera(47,1.6,.001,100),target=new T.Vector3(.23,.64,-.1),j=new VisitJourney();camera.position.set(1.2,.83,2.3);camera.quaternion.setFromEuler(new T.Euler(.2,-.4,.3));camera.up.set(.1,.9,.2).normalize();camera.zoom=1.7;camera.setViewOffset(1600,1000,123,-85,1400,1000);camera.updateProjectionMatrix();
 const pose=()=>({position:camera.position.toArray(),quaternion:camera.quaternion.toArray(),up:camera.up.toArray(),target:target.toArray(),zoom:camera.zoom,view:{...camera.view},projection:camera.projectionMatrix.toArray()});const original=pose(),state={...initialVisit,path:'front',direction,playing:true};
 assert.equal(j.update(camera,target,bounds,state,.016,undefined),false);assert.deepEqual(pose(),original);assert.equal(j.phase,'approach');
 j.update(camera,target,bounds,state,.5,undefined);const moving=pose();assert.notDeepEqual(moving,original);j.update(camera,target,bounds,{...state,playing:false},1,undefined);assert.deepEqual(pose(),moving);
 for(let i=0;i<70&&j.phase!=='done';i++)j.update(camera,target,bounds,state,.5,undefined);
 assert.equal(j.phase,'done');assert.deepEqual(pose(),original);
 }
});
