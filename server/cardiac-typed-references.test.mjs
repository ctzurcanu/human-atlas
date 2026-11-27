import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
import {Group,Mesh,BoxGeometry,MeshBasicMaterial,FrontSide} from 'three';
const compiled=await build({stdin:{contents:"export * from './app/cardiac-label-sites';export * from './app/cardiac-reference-types';export * from './app/cardiac-reference-overlay';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const directory=process.env.ATLAS_CARDIAC_LABEL_BUNDLE_STAGE??'.local-models/ta98-review/bowel-wall-local-atlas';
const expected={'finite-mesh':38,'surface-region':6,'boundary-curve':7,'opening-rim':9,'point-set':2,'space-region':3,collection:6};
const models={};
for(const sex of ['male','female']){
 const raw=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),bundle=JSON.parse(fs.readFileSync(`${directory}/${sex}/cardiac-label-sites.json`));
 const atlas=api.withCardiacLabelBundle(raw,bundle),registry=api.cardiacLabelRegistry(atlas);models[sex]={raw,bundle,atlas,registry};
 test(`${sex}: all 71 references have a validated anatomical dimension`,()=>{
  assert.ok(registry);assert.equal(registry.typedGeometryCount,71);assert.deepEqual(registry.typedGeometryCounts,expected);
  for(const g of registry.termGroups){assert.equal(g.representation,g.geometry.type,g.term);assert.ok(api.validCardiacReferenceGeometry(g.geometry,new Map(atlas.parts.map(p=>[p.id,p])),registry.bankSha256),g.term);}
  assert.strictEqual(atlas.parts,raw.parts);assert.strictEqual(atlas.chunks,raw.chunks);
 });
 test(`${sex}: commissures are points, surface regions use native facets, and mouths stay unfilled`,()=>{
  for(const g of registry.termGroups){const geometry=g.geometry;
   if(geometry.type==='point-set'){assert.equal(geometry.points.length,3);assert.equal(new Set(geometry.points.map(p=>p.join(','))).size,3);}
   if(geometry.type==='surface-region'){assert.ok(geometry.surfaces.every(s=>s.faces.length>0&&s.areaMetres2>0));assert.ok(!('vertices'in geometry));}
   if(geometry.type==='opening-rim'){assert.equal(geometry.filled,false);assert.ok(geometry.paths.every(p=>p.closed&&p.points.length>=3));assert.ok(!('faces'in geometry));}
   if(geometry.type==='space-region'){assert.equal(geometry.filled,false);assert.ok(geometry.vertices.length>=4);}
  }
  const foramen=registry.termGroups.find(g=>g.term==='A12.1.01.007');assert.equal(foramen.geometry.type,'surface-region');assert.equal(foramen.geometry.state,'closed-adult-septal-site');
  for(const term of ['A12.1.00.030','A12.1.01.005','A12.1.01.007']){
   const g=registry.termGroups.find(g=>g.term===term);assert.ok(g.geometry.surfaces.reduce((area,s)=>area+s.areaMetres2,0)>1e-5,`${term}: a macroscopic surface must not collapse to a point-sized facet`);
  }
 });
 test(`${sex}: collective references resolve to real children, with no invented collective solid`,()=>{
  const terms=new Set(api.cardiacLabeledTermParts(atlas).keys());
  for(const g of registry.termGroups.filter(g=>g.geometry.type==='collection'))for(const term of g.geometry.terms)assert.ok(terms.has(term),`${g.term}: ${term}`);
 });
 test(`${sex}: incomplete or out-of-bank typed geometry fails closed`,()=>{
  const missing=structuredClone(bundle);delete missing.integration.termGroups[0].geometry;assert.strictEqual(api.withCardiacLabelBundle(raw,missing),raw);
  const bad=structuredClone(bundle),surface=bad.integration.termGroups.find(g=>g.geometry.type==='surface-region');surface.geometry.surfaces[0].sourceBankSha256='0'.repeat(64);assert.strictEqual(api.withCardiacLabelBundle(raw,bad),raw);
 });
 test(`${sex}: finite references preserve their own depth instead of painting rear faces over the front`,()=>{
  const parent=new Group(),host=new Mesh(new BoxGeometry(),new MeshBasicMaterial()),original=host.matrixWorld.clone(),overlay=api.createCardiacReferenceOverlay(parent,[],atlas);
  overlay.update(registry.termGroups.find(g=>g.term==='A12.1.06.002'),()=>host);
  const volumes=parent.children.filter(o=>o.isMesh);assert.equal(volumes.length,6);
  for(const mesh of volumes){assert.equal(mesh.material.side,FrontSide);assert.equal(mesh.material.depthTest,true);assert.equal(mesh.material.depthWrite,true);}
  let clears=0;const renderer={info:{render:{frame:1}},clearDepth:()=>clears++};
  for(const mesh of volumes)mesh.onBeforeRender(renderer);assert.equal(clears,1);
  for(const mesh of volumes)mesh.onBeforeRender(renderer);assert.equal(clears,1);
  renderer.info.render.frame++;volumes[0].onBeforeRender(renderer);assert.equal(clears,2);
  overlay.update(undefined,()=>host);assert.equal(parent.children.length,0);assert.deepEqual(host.matrixWorld.elements,original.elements);
  overlay.dispose();host.geometry.dispose();host.material.dispose();
 });
}
test('female retains the same canonical reference shapes as the male heart',()=>{
 const male=models.male.registry,female=models.female.registry;
 for(const g of male.termGroups){const other=female.termGroups.find(f=>f.term===g.term);assert.equal(other.geometry.type,g.geometry.type);
  for(const key of ['vertices','faces','paths','points','terms'])if(g.geometry[key])assert.deepEqual(other.geometry[key],g.geometry[key],`${g.term}: ${key}`);
 }
});
