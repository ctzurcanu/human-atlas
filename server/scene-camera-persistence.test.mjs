import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
import * as T from 'three';

const built=await build({stdin:{contents:"export * from './app/scene-camera-persistence';export {captureFrame,applyFrame} from './app/view-framing';export {explicitCameraAction} from './app/camera-intent';export {captureExplosionView,fitExplosionView} from './app/explosion-view';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const state={view:'front',reset:2,focus:1,camera:[.5,.5,1,.5,.5,.5,.5,.45,.35,.5,1,.5]};
function fixture(){
 const camera=new T.PerspectiveCamera(41,1.6,.003,7),target=new T.Vector3(.001,1.25,.03);camera.position.set(-.012,1.265,.02);camera.up.set(.2,.91,.15).normalize();camera.lookAt(target);camera.zoom=1.23;camera.filmGauge=37;camera.filmOffset=.3;camera.setViewOffset(1600,1000,-136,-180,1600,1000);camera.updateProjectionMatrix();camera.updateMatrixWorld(true);
 const pose=api.captureRawSceneCamera(camera,target),intent=api.sceneCameraIntent('local-ta98',3,state,7),record={intent,pose,minDistance:.003,maxDistance:.12,explosion:{amount:.47,velocity:.3,requested:.8,start:api.captureExplosionView(camera,target),fitting:true}};
 return {camera,target,record,owner:{}};
}

test('Same-model restart preserves the full physical pose exactly, independent of selection bounds or panel framing',()=>{
 const {camera,target,record,owner}=fixture();api.rememberSceneCamera(owner,record);const saved=api.resumeSceneCamera(owner,record.intent);assert.ok(saved);
 const rebuilt=new T.PerspectiveCamera(34,1,.01,100),rebuiltTarget=new T.Vector3();api.restoreRawSceneCamera(rebuilt,rebuiltTarget,saved.pose);
 assert.deepEqual(api.captureRawSceneCamera(rebuilt,rebuiltTarget),record.pose);assert.deepEqual(rebuilt.projectionMatrix.elements,camera.projectionMatrix.elements);
 const area={left:.02,right:.67,top:.08,bottom:.92},oldBounds=new T.Box3(new T.Vector3(-.01,1.22,.02),new T.Vector3(.01,1.28,.04)),newBounds=new T.Box3(new T.Vector3(-.03,1.15,-.02),new T.Vector3(.03,1.35,.06));
 const frame=api.captureFrame(camera,target,oldBounds,area),relative=new T.PerspectiveCamera(41,1.6,.003,7),relativeTarget=new T.Vector3();api.applyFrame(relative,relativeTarget,frame,newBounds,{left:.26,right:.98,top:.08,bottom:.92});assert.notDeepEqual(relative.position.toArray(),camera.position.toArray(),'normalized frames cannot substitute for the raw checkpoint');
 for(let reload=0;reload<8;reload++){api.rememberSceneCamera(owner,{...record,pose:api.captureRawSceneCamera(rebuilt,rebuiltTarget)});api.restoreRawSceneCamera(rebuilt,rebuiltTarget,api.resumeSceneCamera(owner,record.intent).pose);assert.deepEqual(api.captureRawSceneCamera(rebuilt,rebuiltTarget),record.pose);}
 assert.notEqual(saved.pose.position,camera.position);assert.notEqual(saved.pose.view,camera.view,'checkpoint owns its projection values');
});

test('Only implicit restarts resume: explicit navigation, model, frame, Reset, view and Center requests invalidate the checkpoint',()=>{
 const {record,owner}=fixture();
 for(const change of [{model:'local-female-ta98'},{navigation:4},{view:'back'},{reset:3},{focus:2},{camera:[...state.camera]},{transition:8}]){
  api.rememberSceneCamera(owner,record);assert.equal(api.resumeSceneCamera(owner,{...record.intent,...change}),undefined,JSON.stringify(change));assert.equal(api.resumeSceneCamera(owner,record.intent),undefined,'an invalidated model/frame cannot return through a later retry');
 }
 api.rememberSceneCamera(owner,record);assert.ok(api.resumeSceneCamera(owner,{...record.intent,camera:undefined,focus:0,transition:0}),'clearing a saved frame or focus is not camera movement');
 assert.equal(api.resumeSceneCamera({},record.intent),undefined,'fresh studio/document has no checkpoint');
 const other={};api.rememberSceneCamera(other,{...record,intent:{...record.intent,model:'embryo'}});assert.equal(api.resumeSceneCamera(owner,record.intent).intent.model,'local-ta98','instances cannot exchange camera state');
});

test('Restart preserves the existing explosion motion and exact original Implode checkpoint',()=>{
 const {record,owner}=fixture();api.rememberSceneCamera(owner,record);const resumed=api.resumeSceneCamera(owner,record.intent);assert.deepEqual(resumed.explosion,record.explosion);
 const camera=new T.PerspectiveCamera(),target=new T.Vector3();api.restoreRawSceneCamera(camera,target,resumed.pose);
 const box=new T.Box3(new T.Vector3(-5,-5,-5),new T.Vector3(5,5,5)),area={left:.2,right:.8,top:.1,bottom:.9};
 for(const amount of [.47,.2,.7,.3,0])api.fitExplosionView(camera,target,resumed.explosion.start,box,area,amount);
 assert.deepEqual(api.captureExplosionView(camera,target),record.explosion.start,'full Implode still returns to the exact pre-explosion pose');
});

test('Actual scene restart initialization consumes old camera intent instead of requesting startup fitting',async()=>{
 const scene=await readFile('app/scene.tsx','utf8'),page=await readFile('app/page.tsx','utf8');
 assert.ok(scene.indexOf('if(restart){restoreRawSceneCamera(')>scene.indexOf('camera.position.set(1.4'),'restoration follows default camera initialization');
 for(const expected of ["lastView=restart?latest.current.view:''","lastReset=restart?latest.current.reset:-1","areaInitialized=!!restart","lastCamera:number[]|undefined=restart?latest.current.camera:undefined","lastFocus=restart?latest.current.focus??0:0","lastExplode=restart?.explosion.requested??0","amount=restart?.explosion.amount??0"]){assert.ok(scene.includes(expected),expected);}
 assert.ok(page.includes('const pop=()=>{cameraNavigationSerial.current++;'),'deliberate URL/saved-view navigation invalidates the checkpoint');
 assert.ok(page.includes('cameraNavigation={cameraNavigationSerial.current}'));
 assert.equal(api.explicitCameraAction(state,{view:state.view,reset:state.reset,focus:state.focus},false),false);
 assert.equal(api.explicitCameraAction(state,{view:'',reset:-1,focus:0},false),true,'a fresh deliberate navigation retains normal initial fitting');
});
