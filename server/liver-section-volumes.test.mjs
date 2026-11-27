import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/liver-section-volumes';export * from './app/stencil-caps';export * from './app/local-ta98-layer-catalogue';export {Scene,Plane,Matrix4,BoxGeometry,Vector3} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const path='.local-models/geometry-workbench/liver-completion/liver-section-volumes-v141/candidate.json';
const source=fs.existsSync(path)?JSON.parse(fs.readFileSync(path)):null;

test('closed liver classifiers require explicit bank, part and coordinate identities',{skip:!source},()=>{
 const registry=api.liverSectionRegistry(source);assert.equal(registry.byId.size,4);assert.equal(registry.exclusions.length,5);
 for(const mutate of [a=>a.liverSectionIntegration.bankSha256='0'.repeat(64),a=>a.liverSectionIntegration.volumes.pop(),a=>a.liverSectionIntegration.exclusions[0].indices=-1,a=>a.liverSectionIntegration.volumes[0].partId='ZA:Liver']){
  const wrong=structuredClone(source);mutate(wrong);assert.equal(api.liverSectionRegistry(wrong),null);
 }
 for(const id of registry.byId.keys())assert.equal(source.parts.find(p=>p.id===id).sectionSolid,false,'open display surfaces are not falsely certified closed');
});

test('subtractive stencils resolve body first and erase overhanging or overlapping channels',()=>{
 const scene=new api.Scene(),planes=[new api.Plane(new api.Vector3(0,0,1),0),new api.Plane(new api.Vector3(0,1,0),100)];
 const caps=api.createStencilCaps(scene,planes,()=>{throw Error('No texture needed');});
 const part={id:'liver-test',name:'Liver',system:'digestive',pairedTissue:true,tissue:'liver',bounds:[[-1,-1,-1],[1,1,1]]};
 const source=new api.BoxGeometry(2,2,2),hole=new api.BoxGeometry(.5,.5,4);
 caps.update([{part,index:0,cutIndex:0,geometry:source,exclusions:[hole,hole],matrix:new api.Matrix4(),mode:'solid',closedWall:true}]);
 const items=scene.children.filter(o=>o.visible).sort((a,b)=>a.renderOrder-b.renderOrder);
 const constants=items.map(o=>({func:o.material.stencilFunc,ref:o.material.stencilRef,read:o.material.stencilFuncMask,write:o.material.stencilWriteMask,op:o.material.stencilZPass}));
 assert.equal(items.length,13,'body count/resolve/clear, two independent channel batches, final fill');
 assert.equal(items.at(-1).material.stencilFuncMask,128);
 assert.equal(items.at(-1).material.stencilRef,128);
 // Emulate the GL stencil read/modify/write rules at representative pixels.
 // Each closed count produces its winding; all quadrilateral operations run.
 const apply=(value,m)=>{const ref=m.stencilRef&m.stencilFuncMask,read=value&m.stencilFuncMask;const pass=m.stencilFunc===519||m.stencilFunc===517&&ref!==read||m.stencilFunc===514&&ref===read;if(!pass)return value;return (value&~m.stencilWriteMask)|(m.stencilRef&m.stencilWriteMask);};
 for(const [body,first,second,expected] of [[1,0,0,true],[1,1,0,false],[1,1,1,false],[0,1,0,false],[0,0,0,false]]){
  let value=body;value=apply(value,items[2].material);value=apply(value,items[3].material);
  value=(value&128)|first;value=apply(value,items[6].material);value=apply(value,items[7].material);
  value=(value&128)|second;value=apply(value,items[10].material);value=apply(value,items[11].material);
  assert.equal((value&128)===128,expected,JSON.stringify({body,first,second,constants}));
 }
 caps.dispose();source.dispose();hole.dispose();
});

test('liver section candidates remain local review routes',()=>{
 const q=new URLSearchParams('hepatic=section-volumes-v141');
 assert.equal(api.localTa98LayerCatalogue('local-ta98',q,false),undefined);
 assert.match(api.localTa98LayerCatalogue('local-ta98',q,true),/liver-section-volumes-v141\/candidate.json$/);
 assert.doesNotMatch(api.localTa98LayerCatalogue('local-female-ta98',q,true),/liver-section-volumes/);
});
