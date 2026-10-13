import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
import * as T from 'three';

const result=await build({stdin:{contents:"export {createSharedTissueWall} from './app/shared-tissue-wall';export {createStencilCaps} from './app/stencil-caps';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {createSharedTissueWall,createStencilCaps}=await import('data:text/javascript;base64,'+Buffer.from(result.outputFiles[0].contents).toString('base64'));
const radii=[1,.9,.8,.7],segments=24;
function fixture(){
 const faces=[],rims=[];
 for(let j=0;j<segments;j++){const k=(j+1)%segments;faces.push(j,segments+j,k,k,segments+j,segments+k);rims.push(j,k,segments+k,segments+j);}
 const wall=createSharedTissueWall({interfaces:radii.map(r=>new Float32Array(Array.from({length:2*segments},(_,i)=>[r*Math.cos(i%segments*2*Math.PI/segments),i<segments?0:2,r*Math.sin(i%segments*2*Math.PI/segments)]).flat())),faces:new Uint32Array(faces),rims:new Uint32Array(rims),layerIds:['externa','media','intima']});
 const candidates=wall.full.map((geometry,index)=>({part:{id:String(index),conceptId:String(index),name:'Finite vascular coat',system:'venous',pairedTissue:true,sectionSolid:true,bounds:[[-1,0,-1],[1,2,1]],tissueAppearance:{color:[[.52,.43,.63],[.62,.33,.43],[.77,.62,.67]][index],pattern:'solid',scale:.01,relief:0}},index,geometry,matrix:new T.Matrix4(),mode:'solid',closedWall:true}));
 return {wall,candidates};
}
function counts(scene,fill){return scene.children.filter(m=>m.visible&&m.material.colorWrite===false&&m.renderOrder===fill.renderOrder-1);}
// Model the actual fragment counters with exact ray/triangle intersections,
// face orientation and each material's live clipping planes. No cap geometry
// or sampled color is used to decide whether tissue is filled.
function winding(meshes,point){
 const ray=new T.Raycaster(new T.Vector3(point.x,3,point.z),new T.Vector3(0,-1,0));let value=0;
 for(const mesh of meshes){for(const hit of ray.intersectObject(mesh,false))if(mesh.material.clippingPlanes.every(p=>p.distanceToPoint(hit.point)>=0))value+=mesh.material.stencilZPass===T.IncrementWrapStencilOp?1:-1;}
 return value;
}
const visibleFills=scene=>scene.children.filter(m=>m.visible&&m.userData.tissueFrame);
function cutCandidates(candidates){return candidates.flatMap(c=>[0,1].map(cutIndex=>({...c,cutIndex})));}

test('Opposing cuts fill each complete finite annulus while leaving lumen and exterior empty',()=>{
 const {wall,candidates}=fixture(),scene=new T.Scene(),planes=[new T.Plane(new T.Vector3(0,1,0),-.7),new T.Plane(new T.Vector3(0,-1,0),.76)],caps=createStencilCaps(scene,planes,()=>new T.Texture());
 caps.update(cutCandidates(candidates));scene.updateMatrixWorld(true);assert.equal(caps.activeCount(),6);
 const fills=visibleFills(scene);assert.equal(fills.length,6);
 for(const fill of fills){const cut=fill.userData.cutIndex,coat=candidates.findIndex(c=>counts(scene,fill)[0].geometry.getIndex()===c.geometry.getIndex()),counter=counts(scene,fill);assert.equal(counter.length,2);assert.ok(coat>=0);
  for(const mesh of counter)assert.deepEqual(mesh.material.clippingPlanes,[planes[cut]],'other-plane clipping must not open the counting source');
  assert.deepEqual(fill.material.clippingPlanes,[planes[1-cut]],'visible cap still obeys the other cut');
  for(let j=0;j<segments;j++){const angle=(j+.37)*2*Math.PI/segments,r=(radii[coat]+radii[coat+1])/2,point=new T.Vector3(r*Math.cos(angle),cut===0?.7:.76,r*Math.sin(angle));assert.equal(Math.abs(winding(counter,point)),1,'all circumference samples of real tissue are covered');}
  assert.equal(winding(counter,new T.Vector3(.2,fill.position.y,.03)),0,'blood lumen remains open');
  assert.equal(winding(counter,new T.Vector3(1.2,fill.position.y,.03)),0,'cap cannot exceed outer wall');
 }
 // The old two-plane predicate counts only the open slab: its two ends cancel.
 const fill=fills[0],counter=counts(scene,fill),point=new T.Vector3(.945,.7,.02);
 assert.equal(Math.abs(winding(counter,point)),1);
 for(const mesh of counter)mesh.material.clippingPlanes=[planes[0],planes[1]];
 assert.equal(winding(counter,point),0,'fixture reproduces the previous unfilled-annulus defect');
 caps.dispose();wall.dispose();
});

test('Two moving cuts keep live source membership, opposing masks and exact cap planes without rebuilding shells',()=>{
 const {wall,candidates}=fixture(),scene=new T.Scene(),planes=[new T.Plane(new T.Vector3(0,1,0),-.7),new T.Plane(new T.Vector3(0,-1,0),.76)],caps=createStencilCaps(scene,planes,()=>new T.Texture());
 caps.update(cutCandidates(candidates));scene.updateMatrixWorld(true);const initialCounters=scene.children.filter(m=>m.visible&&m.material.colorWrite===false),initialFills=visibleFills(scene);
 for(let frame=0;frame<20;frame++){
  const y=.3+frame*.06;planes[0].constant=-y;planes[1].constant=y+.06;caps.update(cutCandidates(candidates));scene.updateMatrixWorld(true);
  assert.ok(caps.planeError()<1e-12);assert.deepEqual(visibleFills(scene),initialFills);assert.deepEqual(scene.children.filter(m=>m.visible&&m.material.colorWrite===false),initialCounters);
  for(const fill of initialFills){const counter=counts(scene,fill),coat=candidates.findIndex(c=>counter[0].geometry.index===c.geometry.index),r=(radii[coat]+radii[coat+1])/2;assert.equal(Math.abs(winding(counter,new T.Vector3(r*.96,fill.position.y,r*.28))),1);}
 }
 caps.dispose();wall.dispose();assert.equal(scene.children.length,0);
});

test('Final fills clear discarded stencil coverage before another cut or coat; other-cut masks remain strict',()=>{
 const {wall,candidates}=fixture(),scene=new T.Scene(),planes=[new T.Plane(new T.Vector3(0,1,0),-.7),new T.Plane(new T.Vector3(1,0,0),0)],caps=createStencilCaps(scene,planes,()=>new T.Texture());
 caps.update(cutCandidates(candidates));scene.updateMatrixWorld(true);
 for(const fill of visibleFills(scene)){
  let clears=0;fill.onAfterRender({clearStencil:()=>clears++});assert.equal(clears,1,'discarded cap pixels cannot leave winding for the next batch');
  assert.equal(fill.material.stencilFunc,T.NotEqualStencilFunc);assert.equal(fill.material.stencilRef,0);assert.equal(fill.material.stencilZPass,T.ReplaceStencilOp);
  assert.equal(fill.material.transparent,false,'counts, fill and clear remain in one opaque sequence');
 }
 const fill=visibleFills(scene).find(m=>m.userData.cutIndex===0),counter=counts(scene,fill);
 for(const x of [-.94,.94]){const point=new T.Vector3(x,.7,.03);assert.equal(Math.abs(winding(counter,point)),1);assert.equal(fill.material.clippingPlanes.every(p=>p.distanceToPoint(point)>=0),x>0,'the other cut admits only its retained side');}
 caps.dispose();wall.dispose();
});

test('Ordinary solid rim counts finish their inner fill before clearing the shared stencil',()=>{
 const scene=new T.Scene(),geometry=new T.BoxGeometry(2,2,2),planes=[new T.Plane(new T.Vector3(0,1,0),0),new T.Plane(new T.Vector3(1,0,0),0)],caps=createStencilCaps(scene,planes,()=>new T.Texture());
 caps.update([{part:{id:'muscle',name:'Muscle',system:'muscular',bounds:[[-1,-1,-1],[1,1,1]]},index:0,cutIndex:0,geometry,matrix:new T.Matrix4(),mode:'solid'}]);
 const quads=scene.children.filter(m=>m.visible&&m.material.colorWrite!==false).sort((a,b)=>a.renderOrder-b.renderOrder);assert.equal(quads.length,2);
 let clears=0;const renderer={clearStencil:()=>clears++};quads[0].onAfterRender(renderer);assert.equal(clears,0,'rim does not clear before inset coverage');quads[1].onAfterRender(renderer);assert.equal(clears,1);
 caps.dispose();geometry.dispose();
});
