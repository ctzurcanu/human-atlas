import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
import {build} from 'esbuild';

const root=process.cwd(),stage=process.env.PULMONARY_INTERFACE_STAGE??root+'/.local-models/geometry-workbench/biodigital-male-circulation-v1/biodigital-heart-v1/pulmonary-root-interface-v1/reciprocal-cap-ownership-v1';
const manifest=JSON.parse(fs.readFileSync(stage+'/manifest.json')),layout=JSON.parse(fs.readFileSync(stage+'/host-geometry.json')),raw=fs.readFileSync(stage+'/host-geometry.bin');
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
assert.equal(sha(raw),layout.bankSha256);
const scene=fs.readFileSync('app/scene.tsx','utf8'),lines=scene.split('\n');
const opacityLine=lines.find(l=>l.includes('const displayedTissueOpacity=(i:number)=>'));
const groupLine=lines.find(l=>l.includes('for(const {index,geometry} of tissueDisplays)'));
assert.ok(opacityLine&&groupLine,'Run the real current scene group updater, not a copied policy implementation');
const imports=`import * as T from 'three';export {T};export {anatomyPartGeometry,anatomyGeometryView} from './app/shared-tissue-geometry';export {sectionCapGeometry} from './app/section-cap';export {sectionCapProfile} from './app/section-cap-profile';export {createStencilCaps} from './app/stencil-caps';export {sourceCardiacIntegration} from './app/biodigital-heart';import {tissueNeighborOccludes,tissuePosesCoincide} from './app/tissue-interface-visibility';\n`;
const harness=`export function updateActualSceneGroups(atlas,pickers,tissueDisplays,data,selectedData,contextUniform){const partIndices=new Map(atlas.parts.map((p,i)=>[p.id,i]));const hepaticRegionReplaces=()=>{throw Error('No hepatic replacement is permitted in this source-bound test');};${opacityLine}\n${groupLine}\nreturn atlas.parts.map((_,i)=>displayedTissueOpacity(i));}`;
const built=await build({stdin:{contents:imports+harness,resolveDir:root,loader:'ts'},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64')),{T}=api;
const full=JSON.parse(fs.readFileSync(manifest.baseProfile));assert.equal(sha(fs.readFileSync(manifest.baseProfile)),manifest.baseProfileSha256);
const atlas={...full,parts:layout.parts.map(p=>p.record)},pickers=[],views=[],display=[],material=new T.MeshBasicMaterial({side:T.DoubleSide});
const slices=(r,Ctor)=>new Ctor(raw.buffer,raw.byteOffset+r.byteOffset,r.scalarCount);
for(let i=0;i<layout.parts.length;i++){
 const fixture=layout.parts[i],p=fixture.record,g=api.anatomyPartGeometry(p,i,layout.parts.length);
 g.setAttribute('position',new T.BufferAttribute(slices(fixture.positions,Float32Array),3));g.setAttribute('normal',new T.BufferAttribute(slices(fixture.normals,Float32Array),3));g.setIndex(new T.BufferAttribute(slices(fixture.indices,Uint32Array),1));g.computeBoundingBox();g.computeBoundingSphere();
 const picker=new T.Mesh(g,material);picker.position.fromArray(p.geometryOrigin);picker.updateMatrixWorld(true);pickers.push(picker);const view=api.anatomyGeometryView(g);views.push(view);const mesh=new T.Mesh(view,[material]);mesh.position.copy(picker.position);mesh.updateMatrixWorld(true);display.push(mesh);
}
const tissueDisplays=views.map((geometry,index)=>({geometry,index}));const results=[],render={states:[],sections:[],centre:layout.centre_m,basis:layout.basis,normal:layout.normal};
const centre=new T.Vector3(...layout.centre_m),normal=new T.Vector3(...layout.normal),U=new T.Vector3(...layout.basis.map(r=>r[0])),V=new T.Vector3(...layout.basis.map(r=>r[1]));
const capCounts=()=>views.map((view,i)=>{const cap=new Set(layout.parts[i].sourceCapFaceIds);let n=0;for(const g of view.groups)for(let f=g.start/3;f<(g.start+g.count)/3;f++)if(cap.has(f))n++;return n;});
const capTriangles=()=>views.map((view,i)=>{const cap=new Set(layout.parts[i].sourceCapFaceIds),tri=[];for(const group of view.groups)for(let f=group.start/3;f<(group.start+group.count)/3;f++)if(cap.has(f)){tri.push([0,1,2].map(k=>{const p=new T.Vector3().fromBufferAttribute(view.getAttribute('position'),view.index.getX(f*3+k)).applyMatrix4(pickers[i].matrixWorld).sub(centre);return [p.dot(U),p.dot(V)];}));}return tri;});
function state({visible=[1,1],selected=[1,1],context=1,detached=false,rotated=false}){
 const data=new Float32Array(8),selectedData=new Uint8Array(8);
 for(let i=0;i<2;i++){
  data[i*4+3]=visible[i];selectedData[i*4]=selected[i]?255:0;selectedData[i*4+1]=255;
  const q=rotated?new T.Quaternion().setFromEuler(new T.Euler(.13,-.22,.09)):new T.Quaternion(),origin=new T.Vector3(...atlas.parts[i].geometryOrigin).applyQuaternion(q);if(rotated)origin.add(new T.Vector3(.01,-.02,.03));if(detached&&i===1)origin.x+=.01;
  for(const mesh of [pickers[i],display[i]]){mesh.position.copy(origin);mesh.quaternion.copy(q);mesh.updateMatrixWorld(true);}
 }
 return api.updateActualSceneGroups(atlas,pickers,tissueDisplays,data,selectedData,{value:context});
}
async function check(name,run){await test(name,()=>{const began=performance.now();const detail=run();results.push({name,passed:true,elapsedMs:performance.now()-began,detail});});}
await check('Exact source face ranges partition unchanged buffers and add no cap disk',()=>{
 for(const [i,p]of atlas.parts.entries()){let at=0;for(const g of p.tissueDrawGroups){assert.equal(g.start,at);assert.equal(g.start%3,0);assert.equal(g.count%3,0);at+=g.count;}assert.equal(at,p.indexCount);assert.equal(p.tissueDrawGroups.filter(g=>g.neighbor).reduce((n,g)=>n+g.count,0),198);assert.equal(views[i].index,pickers[i].geometry.index);assert.equal(views[i].getAttribute('position'),pickers[i].geometry.getAttribute('position'));assert.equal(views[i].getAttribute('normal'),pickers[i].geometry.getAttribute('normal'));}
 return {capFacesEach:66,computationalFaces:pickers.map(p=>p.geometry.index.count/3),sameIndexAndAttributeObjects:true};
});
const cases=[['both opaque',{},[0,0]],['RV isolated',{visible:[1,0]},[66,0]],['PT isolated',{visible:[0,1]},[0,66]],['equal faint context',{selected:[0,0],context:.18},[0,0]],['RV stronger than PT context',{selected:[1,0],context:.18},[66,0]],['PT stronger than RV context',{selected:[0,1],context:.18},[0,66]],['detached exploded peers',{detached:true},[66,66]],['common rotated pose',{rotated:true},[0,0]]];
for(const [name,settings,expected]of cases)await check('Actual scene drawgroups: '+name,()=>{const opacity=state(settings),counts=capCounts();assert.deepEqual(counts,expected);render.states.push({name,opacity,capFaces:counts,capTriangles:capTriangles()});return {opacity,visibleCapFaces:counts,totalDisplayFaces:views.map(g=>g.groups.reduce((n,r)=>n+r.count/3,0))};});
await check('Existing source cardiac registry and typed label guards survive drawgroup metadata',()=>{const changed=structuredClone(full);for(const fixture of layout.parts)Object.assign(changed.parts.find(p=>p.id===fixture.id),{tissueDrawGroups:fixture.record.tissueDrawGroups});assert.ok(api.sourceCardiacIntegration(full));const registry=api.sourceCardiacIntegration(changed);assert.ok(registry);if(manifest.expectedSourceCardiacTerms)assert.equal(registry.termGroups.length,manifest.expectedSourceCardiacTerms);return {registryRetained:true,typedReferences:registry.referenceSites.length,termGroups:registry.termGroups.length};});
await check('Full computational meshes retained; actual CPU and GPU-counter section centres stay unfilled',()=>{
 state({});const detail=[];
 for(const offset of [-.001,0,.001]){
  const point=centre.clone().addScaledVector(normal,offset),plane=new T.Plane().setFromNormalAndCoplanarPoint(normal,point),row={offset_m:offset,cpu:[]};
  for(let i=0;i<2;i++){
   const cap=api.sectionCapGeometry(pickers[i].geometry,plane,pickers[i].matrixWorld,api.sectionCapProfile(atlas.parts[i]));let hits=0,triangles=0,tri=[];
   if(cap){const mesh=new T.Mesh(cap,material);mesh.updateMatrixWorld(true);const ray=new T.Raycaster(point.clone().addScaledVector(normal,.05),normal.clone().negate());hits=ray.intersectObject(mesh,false).length;const pos=cap.getAttribute('position'),ind=cap.index;triangles=(ind?ind.count:pos.count)/3;for(let f=0;f<triangles;f++)tri.push([0,1,2].map(k=>{const p=new T.Vector3().fromBufferAttribute(pos,ind?ind.getX(f*3+k):f*3+k).sub(centre);return [p.dot(U),p.dot(V)];}));cap.dispose();}
   assert.equal(hits,0,'No CPU cap triangle may span the actual annular centre');row.cpu.push({part:atlas.parts[i].id,triangles,centralRayHits:hits});render.sections.push({offset_m:offset,owner:i,triangles:tri});
  }
  const scene=new T.Scene(),planes=[plane,new T.Plane(new T.Vector3(0,1,0),10000)],stencil=api.createStencilCaps(scene,planes,()=>new T.Texture());const candidates=pickers.map((picker,i)=>{const profile=api.sectionCapProfile(atlas.parts[i]);return {part:atlas.parts[i],index:i,cutIndex:0,geometry:picker.geometry,matrix:picker.matrixWorld,mode:profile.hollowWall?'wall':'solid',wallWidth:profile.width,closedWall:false,opacity:1};});stencil.update(candidates);const hit=stencil.pick(new T.Ray(point.clone().addScaledVector(normal,.05),normal.clone().negate()),()=>true);assert.equal(hit,undefined,'Actual stencil counter parity/picker must leave the lumen centre empty');row.gpuCounterCentrePick=false;row.activeGpuCounters=stencil.activeCount();stencil.dispose();detail.push(row);
 }
 return {sourceFullIndexCountsRetained:pickers.map(p=>p.geometry.index.count),cuts:detail,gpuQualification:'Executes actual GPU counter geometry and parity-picking path without a WebGL raster context; CPU cap triangles are separately rendered.'};
});
await check('Actual full picker versus assembled display: targeted hidden annulus rays',()=>{
 state({});const targets=[];
 for(let i=0;i<2;i++){const g=pickers[i].geometry;for(const f of layout.parts[i].sourceCapFaceIds){const target=new T.Vector3();for(let k=0;k<3;k++)target.add(new T.Vector3().fromBufferAttribute(g.getAttribute('position'),g.index.getX(f*3+k)));target.multiplyScalar(1/3).applyMatrix4(pickers[i].matrixWorld);targets.push(target);}}
 const dirs=[normal,normal.clone().negate(),U,U.clone().negate(),V,V.clone().negate()],wrong=[],hiddenAhead=[];let count=0;
 for(const target of targets)for(const direction of dirs){const ray=new T.Raycaster(target.clone().addScaledVector(direction,.15),direction.clone().negate());const first=meshes=>meshes.flatMap((mesh,i)=>ray.intersectObject(mesh,false).map(hit=>({...hit,owner:i}))).sort((a,b)=>a.distance-b.distance)[0];const full=first(pickers),shown=first(display);count++;if(!full||!shown)continue;const hidden=layout.parts[full.owner].sourceCapFaceIds.includes(full.faceIndex)&&!views[full.owner].groups.some(g=>full.faceIndex*3>=g.start&&full.faceIndex*3<g.start+g.count);if(hidden&&full.distance+8.2e-9<shown.distance){const witness={fullOwner:atlas.parts[full.owner].id,displayOwner:atlas.parts[shown.owner].id,fullFace:full.faceIndex,displayFace:shown.faceIndex,distanceGap_m:shown.distance-full.distance,origin:ray.ray.origin.toArray(),direction:ray.ray.direction.toArray()};hiddenAhead.push(witness);if(full.owner!==shown.owner)wrong.push(witness);}}
 // The full computational picker is intentionally retained; record any actual
 // hidden-face selection conflict without changing its global policy.
 assert.equal(wrong.length,0,'Hidden annular display faces must not select a different owner ahead of the displayed wall');
 return {targetedRayCount:count,invisibleCapAheadOfDisplayedSurface:hiddenAhead.length,wrongOwnerSelections:wrong.length,witnesses:wrong.slice(0,8),geometryPrerequisitesUnresolved:!manifest.sourceGeometryPrerequisitesResolved,pickingTieNumericalAllowance_m:8.2e-9};
});
fs.writeFileSync(stage+'/host-results.json',JSON.stringify({schema:'actual-host-pt-rv-annular-interface-tests/v2',checks:results,profileSha256:manifest.baseProfileSha256,sourceCandidateSha256:manifest.sourcePulmonaryCandidateSha256,actualPackedBankSha256:manifest.actualPackedBankSha256,sceneUpdaterSha256:sha(Buffer.from(opacityLine+'\n'+groupLine)),hostSourceBindings:manifest.hostSourceBindings,allRecordedChecksPassed:results.length===12&&results.every(x=>x.passed),noAppOrRouteChanges:true,geometryCandidateInstallable:manifest.geometryCandidateInstallable,anatomicalAcceptance:false},null,2)+'\n');
fs.writeFileSync(stage+'/host-render-data.json',JSON.stringify(render)+'\n');
