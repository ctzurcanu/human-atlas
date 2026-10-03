import assert from 'node:assert/strict';import test from 'node:test';import {readFile} from 'node:fs/promises';import {existsSync} from 'node:fs';import {build} from 'esbuild';import {Box3,Vector3,PerspectiveCamera} from 'three';
const bundle=await build({entryPoints:['app/explosion-plane.ts','app/hierarchical-explosion.ts','app/anatomy-hierarchy.ts','app/view-framing.ts','app/explosion-view.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/explosion-viewport-test'});
const load=async name=>import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles.find(f=>f.path.endsWith(name+'.js')).contents).toString('base64'));
const {viewPlaneLayout}=await load('explosion-plane'),{createHierarchicalExplosionLayout,explosionBoundsAt}=await load('hierarchical-explosion'),{hierarchyEntries}=await load('anatomy-hierarchy'),{setFrameOffset}=await load('view-framing'),{captureExplosionView,fitExplosionView}=await load('explosion-view');
const catalogs=['public/models/atlas-male-complete.json','public/models/atlas.json','public/models/atlas-hra-female.json','public/models/atlas-embryo.json','public/models/atlas-cell.json','public/models/embryo-3month/atlas.json',...['ta98-male','ta98-female-runtime','reference','male','female'].map(name=>`.local-models/${name}.json`).filter(existsSync)];
for(const model of catalogs)test(`${model}: explode and implode remain inside the usable viewport`,async()=>{
 const atlas=JSON.parse(await readFile(model)),ids=new Set(atlas.parts.filter(p=>!p.suppressed).map(p=>p.id));
 for(const aspect of [1.78,.6]){
  const area={left:.23,right:.75,top:.08,bottom:.9},camera=new PerspectiveCamera(34,aspect,.00001,100),target=new Vector3(.06,.8,.11);
  camera.position.copy(target).add(new Vector3(.3,.2,.5));camera.up.set(.15,1,.2).normalize();camera.lookAt(target);setFrameOffset(camera,area);camera.updateMatrixWorld(true);
  const orientation=camera.quaternion.clone(),originalTarget=target.clone(),checkpoint=captureExplosionView(camera,target),layout=createHierarchicalExplosionLayout(atlas,hierarchyEntries(atlas),ids,'systems',aspect);
  for(const direction of [1,-1])for(let frame=0;frame<=24;frame++){
   const amount=direction===1?frame/24:1-frame/24,mapped=viewPlaneLayout(layout,atlas.parts,camera.quaternion),current=explosionBoundsAt(mapped,atlas.parts,ids,amount);
   const box=new Box3(new Vector3(...current.min),new Vector3(...current.max));
   fitExplosionView(camera,target,checkpoint,box,area,amount);
   assert.ok(1-Math.abs(camera.quaternion.dot(orientation))<1e-10,'Fitting preserves orientation');assert.ok(target.equals(originalTarget),'Fitting preserves orbit target');
   if(amount===0){assert.deepEqual(captureExplosionView(camera,target),checkpoint,'The original view returns exactly');continue;}
   for(let corner=0;corner<8;corner++){const point=new Vector3(box[(corner&1)?'max':'min'].x,box[(corner&2)?'max':'min'].y,box[(corner&4)?'max':'min'].z).project(camera),x=(point.x+1)/2,y=(1-point.y)/2;assert.ok(x>=area.left-1e-8&&x<=area.right+1e-8&&y>=area.top-1e-8&&y<=area.bottom+1e-8&&point.z>=-1&&point.z<=1,`Frame ${frame}, direction ${direction}: every moving extent stays visible`);}
  }
 }
});

test('Lossless return after partial reversals with cropped pan, roll, zoom and low projection center',()=>{
 const camera=new PerspectiveCamera(41,1.6,.001,100),target=new Vector3(2,3,4);
 camera.position.set(2.1,3.2,4.3);camera.up.set(.3,.7,.2).normalize();camera.lookAt(target);camera.zoom=1.7;camera.setViewOffset(1600,1000,130,-180,1600,1000);camera.updateProjectionMatrix();
 const initial=captureExplosionView(camera,target),box=new Box3(new Vector3(-5,-5,-5),new Vector3(5,5,5)),area={left:.2,right:.8,top:.1,bottom:.9};
 for(const amount of [.001,.3,.8,.4,.7,.2,0])fitExplosionView(camera,target,initial,box,area,amount);
 assert.deepEqual(captureExplosionView(camera,target),initial);
});
