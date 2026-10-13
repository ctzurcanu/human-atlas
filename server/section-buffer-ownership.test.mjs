import test from 'node:test';import assert from 'node:assert/strict';import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/stencil-section';export {createStencilCaps} from './app/stencil-caps';export {BoxGeometry,Scene,Plane,Vector3,Matrix4,Texture} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
test('disposing a section counter never releases source GPU attributes, but still releases owned indices',()=>{
 const source=new api.BoxGeometry(.1,.1,.1),position=source.attributes.position,index=source.index;
 for(const closed of [true,false]){
  const shell=api.stencilShell(source,0,closed);assert.equal(shell.attributes.position,position);
  const ownIndex=closed?null:shell.index;let fired=false;
  shell.addEventListener('dispose',()=>{fired=true;assert.equal(shell.getAttribute('position'),undefined);assert.equal(shell.index,ownIndex);});
  api.disposeStencilShell(shell);assert.ok(fired);
  assert.equal(source.attributes.position,position);assert.equal(source.index,index);
 }
 source.dispose();
});
test('removing a finite-wall section preserves the display position/index arrays during the actual cap-cache release',()=>{
 const source=new api.BoxGeometry(.1,.1,.1),scene=new api.Scene(),planes=[new api.Plane(new api.Vector3(0,0,-1),0),new api.Plane(new api.Vector3(0,1,0),100)];
 const part={id:'wall',conceptId:'wall',name:'Gallbladder',system:'digestive',pairedTissue:true,sectionAssembly:'same-source',bounds:[[-.05,-.05,-.05],[.05,.05,.05]],tissueAppearance:{color:[.75,.55,.37],pattern:'solid',scale:1,relief:0}};
 const caps=api.createStencilCaps(scene,planes,()=>new api.Texture());
 for(let cycle=0;cycle<3;cycle++){
  caps.update([{part,index:0,cutIndex:0,geometry:source,matrix:new api.Matrix4(),mode:'solid',closedWall:true}]);
  const counter=scene.children.find(m=>m.material?.colorWrite===false&&m.visible);let fired=false;
  counter.geometry.addEventListener('dispose',()=>{fired=true;assert.notEqual(counter.geometry.attributes.position,source.attributes.position);assert.notEqual(counter.geometry.index,source.index);});
  caps.update([]);assert.ok(fired);assert.equal(source.index.count,36);assert.equal(source.attributes.position.count,24);
 }
 caps.dispose();source.dispose();assert.equal(scene.children.length,0);
});
