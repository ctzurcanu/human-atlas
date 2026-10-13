import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const result=await build({stdin:{contents:"export * from './app/hepatic-appearance';export {createStencilCaps} from './app/stencil-caps';export {MeshStandardMaterial,MeshBasicMaterial,ShaderLib,Texture,BoxGeometry,Scene,Plane,Vector3,Matrix4} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(result.outputFiles[0].contents).toString('base64'));
test('gross liver styling is confined to the qualified whole native liver in each actual TA98 model',()=>{
 for(const sex of ['male','female']){
  const atlas=JSON.parse(readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
  const chosen=atlas.parts.filter(p=>api.hepaticSurfaceFinish(p));assert.equal(chosen.length,1);
  assert.equal(chosen[0].id,sex==='male'?'ZA:Liver:Organ':'LOCAL:female:liver:1032');
  for(const part of atlas.parts.filter(p=>p.id!==chosen[0].id&&p.material===chosen[0].material))
   assert.equal(api.hepaticSurfaceFinish(part),undefined,'shared gallbladder, duct and other organ materials are not restyled');
  assert.equal(api.hepaticSurfaceFinish({...chosen[0],suppressed:true}),undefined);
  assert.equal(api.hepaticSurfaceFinish({...chosen[0],sourceId:'unqualified-donor'}),undefined);
 }
});
test('both real Three.js shader paths apply bounded gross tissue grain before lighting and preserve cap stencil state',()=>{
 for(const type of ['basic','standard']){
  const m=type==='basic'?new api.MeshBasicMaterial():new api.MeshStandardMaterial();
  m.stencilWrite=true;m.stencilRef=9;m.clippingPlanes=[new api.Plane(new api.Vector3(0,0,1),-.2)];
  const before={ref:m.stencilRef,write:m.stencilWrite,planes:m.clippingPlanes};
  api.applyHepaticCut(m);
  const shader={vertexShader:api.ShaderLib[type].vertexShader,fragmentShader:api.ShaderLib[type].fragmentShader,uniforms:{}};
  m.onBeforeCompile(shader,{});
  assert.ok(shader.uniforms.hepaticCutColor);assert.ok(shader.vertexShader.includes('hepaticCutPoint=(modelMatrix*vec4(transformed,1.0)).xyz'));
  const field=shader.fragmentShader.indexOf('diffuseColor.rgb=hepaticCutColor*hepaticGrain');
  const lighting=shader.fragmentShader.indexOf(type==='basic'?'reflectedLight.indirectDiffuse *= diffuseColor.rgb':'#include <lights_physical_fragment>');
  assert.ok(field>0&&field<lighting);assert.ok(shader.fragmentShader.includes('fwidth(q)'));
  assert.equal(shader.fragmentShader.includes('sampler2D'),api.ShaderLib[type].fragmentShader.includes('sampler2D'));
  assert.deepEqual({ref:m.stencilRef,write:m.stencilWrite,planes:m.clippingPlanes},before);
  m.dispose();
 }
});
test('a closed liver cap uses only its actual boundary: two coverage counts, one fill, no invented capsule offset or image',()=>{
 const geometry=new api.BoxGeometry(.1,.1,.1),original=Array.from(geometry.attributes.position.array);
 const scene=new api.Scene(),planes=[new api.Plane(new api.Vector3(0,0,-1),0),new api.Plane(new api.Vector3(0,1,0),100)];
 const requested=[];const caps=api.createStencilCaps(scene,planes,url=>{requested.push(url);return new api.Texture();});
 const p={id:'LOCAL:female:liver:1032',conceptId:'LOCAL:female:liver',sourceId:'liver',name:'Liver',system:'digestive',bounds:[[-.05,-.05,-.05],[.05,.05,.05]]};
 for(const position of [-.03,0,.03]){
  planes[0].constant=position;caps.update([{part:p,index:0,cutIndex:0,geometry,matrix:new api.Matrix4(),mode:'solid',closedWall:true}]);
  const shown=scene.children.filter(m=>m.visible);
  assert.equal(shown.filter(m=>m.material.colorWrite===false).length,2);
  assert.equal(shown.filter(m=>m.material.colorWrite!==false).length,1);
  assert.ok(caps.planeError()<1e-12);
 }
 assert.deepEqual(requested,[]);assert.deepEqual(Array.from(geometry.attributes.position.array),original);
 caps.dispose();geometry.dispose();assert.equal(scene.children.length,0);
});
