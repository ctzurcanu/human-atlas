import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {PerspectiveCamera,Vector3} from 'three';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/scene-world-pose';export {compactViewUrl,viewParameters} from './app/view-url';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {captureSceneWorldPose,restoreSceneWorldPose,validSceneWorldPose,readWorldView,writeWorldView,compactViewUrl,viewParameters}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const frame=[.5,.5,.9,.5,.5,.5,.5,.5,.5,.5,1,.5];
function source(){const camera=new PerspectiveCamera(34,1.6,.001,100),target=new Vector3(-1.99e-9,1.69511857,.05273111);camera.position.set(-1.99e-9,1.701201146,.356859901);camera.lookAt(target);camera.zoom=1.24;camera.filmOffset=.03;camera.setViewOffset(1.6,1,-.24,.03,1.6,1);camera.updateProjectionMatrix();return{camera,target,pose:{model:'local-ta98',camera:frame,up:camera.up.toArray(),rotation:[0,0,0,1],world:captureSceneWorldPose(camera,target)}};}
function close(a,b,label){assert.equal(a.length,b.length);a.forEach((value,i)=>assert.ok(Math.abs(value-b[i])<1e-13,`${label}[${i}]: ${value} vs ${b[i]}`));}

test('A compressed saved URL restores the actual world camera independently of changed model bounds',()=>{
 const {camera,target,pose}=source(),url=new URL('https://atlas.example/?model=local-ta98&atlasSource=%2Flocal-models%2Fcurrent.json&frame='+frame.join(','));writeWorldView(url,'local-ta98',pose);const saved=compactViewUrl(url.href),restored=readWorldView(viewParameters(saved),'local-ta98');assert.ok(restored);assert.deepEqual(restored,pose);assert.ok(saved.includes('#atlas-view='));
 const next=new PerspectiveCamera(34,1.6,.01,500),nextTarget=new Vector3(0,.85821,0);next.position.set(0,3.90196,152.18755);next.lookAt(nextTarget);restoreSceneWorldPose(next,nextTarget,restored.world);close(next.position.toArray(),camera.position.toArray(),'position');close(nextTarget.toArray(),target.toArray(),'target');close(next.quaternion.toArray(),camera.quaternion.toArray(),'quaternion');assert.equal(next.zoom,camera.zoom);assert.equal(next.fov,camera.fov);assert.equal(next.filmOffset,camera.filmOffset);assert.deepEqual(next.view,camera.view);close(next.projectionMatrix.elements,camera.projectionMatrix.elements,'projection');assert.ok(Math.abs(next.position.distanceTo(nextTarget)-camera.position.distanceTo(target))<1e-14);
});

test('A world snapshot is model-bound and never carries an old male source into a chosen female model',()=>{
 const {pose}=source(),url=new URL('https://atlas.example/?model=female');writeWorldView(url,'female',pose);assert.equal(url.searchParams.get('world'),null);const params=new URLSearchParams({world:JSON.stringify(pose)});assert.equal(readWorldView(params,'female'),null);assert.deepEqual(readWorldView(params,'local-ta98'),pose);assert.equal(readWorldView(new URLSearchParams(),'local-ta98'),null,'Old normalized-only links remain valid fallbacks');
});

test('Malformed world cameras and projection offsets are rejected without changing source data',()=>{
 const {pose}=source(),before=structuredClone(pose);for(const world of [{...pose.world,position:[0,NaN,1]},{...pose.world,target:pose.world.position},{...pose.world,up:[0,0,0]},{...pose.world,quaternion:[0,0,0,0]},{...pose.world,near:-1},{...pose.world,zoom:0},{...pose.world,fov:180},{...pose.world,aspect:0},{...pose.world,filmOffset:Infinity},{...pose.world,view:{...pose.world.view,width:0}}])assert.equal(validSceneWorldPose(world),false);assert.equal(readWorldView(new URLSearchParams({world:'not-json'}),'local-ta98'),null);assert.equal(readWorldView(new URLSearchParams({world:'x'.repeat(6001)}),'local-ta98'),null);assert.deepEqual(pose,before);
});

test('Scene wiring captures standalone world poses and restores once after startup, without pinning later gestures',async()=>{
 const scene=await readFile('app/scene.tsx','utf8'),page=await readFile('app/page.tsx','utf8');assert.match(scene,/world:captureSceneWorldPose\(camera,controls\.target\)/);assert.match(scene,/saved\?\.world&&saved\.model===sceneModel&&ready&&!following/);assert.match(scene,/restoreSceneWorldPose\(camera,controls\.target,saved\.world\)/);assert.match(scene,/controls\.addEventListener\('start',discardSavedPose\)/);assert.match(scene,/syncSelectedRotation\(\);discardSavedPose\(\);dirty=true/);assert.match(page,/writeWorldView\(url,model,savedWorldPoseRef\.current\?\?scenePoseRef\.current\)/);assert.match(page,/savedPoseRef=\{savedWorldPoseRef\}/);
});

test('Shared world motion preserves each pane projection instead of recomputing it when panels restore',async()=>{
 const scene=await readFile('app/scene.tsx','utf8'),start=scene.indexOf('const ownProjection='),block=scene.slice(start,scene.indexOf('lastComparisonPose=compared',start));assert.match(block,/savedUrlPose\.current\?\.current\?\.world\?\.view/);assert.match(block,/if\(ownProjection!==undefined\)camera\.view=ownProjection\?\{\.\.\.ownProjection\}:null/);assert.doesNotMatch(block,/setFrameOffset/);
});
