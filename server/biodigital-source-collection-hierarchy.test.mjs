import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/anatomy-hierarchy'; export * from './app/system-tree'; export * from './app/hierarchy-choice'; export * from './app/hierarchy-visibility'; export * from './app/organ-models'; export * from './app/ta98-term-bindings'; export * from './app/hierarchical-explosion';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const term='A12.2.12.038',anchor='source:left-hepatic',collection='source:lateral-segmental';
function fixture(){
 const parts=['II','III'].map((segment,i)=>({id:`source:segment:${segment}`,conceptId:`source:segment:${segment}`,sourceId:`original:segment:${segment}`,sourceGroups:[anchor],name:`Branch of left hepatic artery to segment ${segment}`,system:'arterial',chunk:0,positions:i*120,normals:i*120+36,indices:i*120+72,vertexCount:3,indexCount:3,bounds:[[0,1,0],[.01,1.01,.01]]}));
 const family={id:collection,name:'Lateral segmental artery',ta98Term:term,ta98Kind:'exact',elements:parts.map(p=>p.id),partOf:anchor,sourceCollectionId:collection,coverageLimitation:'Curated named original source members; complete anatomy and patency unverified.'};
 return {sex:'male',parts,concepts:[{id:anchor,name:'Left hepatic artery',elements:parts.map(p=>p.id)},...parts.map(p=>({id:p.conceptId,name:p.name,elements:[p.id],partOf:anchor})),family],chunks:[{url:'source.bin',bytes:400}],materials:{},defaultHidden:[],sourceCollectionIdentities:[{id:collection,term,sourceId:'original:left-hepatic',sourceAnchorConceptId:anchor,identityOnly:true,members:parts.map(p=>({partId:p.id,sourceId:p.sourceId,conceptId:p.conceptId}))}],ta98TermBindings:[{term,conceptId:collection,partIds:parts.map(p=>p.id),kind:'exact',representation:'source-collection-identity',sourceCollectionId:collection}]};
}
function entries(nodes){return nodes.flatMap(n=>n.kind==='entry'?[n.entry]:n.kind==='bilateral'?n.entries:entries(n.nodes));}
test('guarded collection has one Systems choice while retaining every original leaf and mesh identity',()=>{
 const a=fixture(),before=structuredClone(a),all=api.hierarchyEntries(a),family=all.filter(e=>e.id===collection);
 assert.equal(family.length,1);assert.equal(family[0].ta98Term,term);assert.deepEqual(family[0].parts,a.parts);
 assert.deepEqual(all.filter(e=>e.id!==collection).map(e=>e.id),a.parts.map(p=>p.conceptId));
 const tree=api.buildSystemTree(all,'systems',a),branch=tree.find(b=>b.id==='circulatory'),actual=entries(branch.nodes).filter(e=>e.id===collection);
 assert.equal(actual.length,1);assert.equal(api.meshCount(branch.parts),2);assert.deepEqual(actual[0].parts.map(p=>p.id),a.parts.map(p=>p.id));
 const state={visible:['arterial'],hidden:[],selected:[],isolate:false,camera:{position:[1,2,3]},cameraIntent:9},hidden=api.toggleSystemParts(state,actual[0].parts,a.parts);
 assert.deepEqual(hidden.hidden.sort(),a.parts.map(p=>p.id).sort());assert.equal(hidden.camera,state.camera);assert.equal(hidden.cameraIntent,9);assert.deepEqual(a,before);
});
test('Systems applies the unchanged per-binding guard to forged, hidden, partial and cross-system collections',()=>{
 const mutations=[a=>a.sourceCollectionIdentities[0].members[0].sourceId='forged',a=>a.parts[0].suppressed=true,a=>a.defaultHidden=[a.parts[0].id],a=>a.ta98TermBindings[0].partIds.pop(),a=>a.parts[0].sourceGroups=[],a=>a.parts[0].system='cardiac'];
 for(const mutate of mutations){const a=fixture();mutate(a);assert.equal(api.catalogueTa98TermParts(a).has(term),false);assert.equal(api.hierarchyEntries(a).some(e=>e.id===collection),false);}
 const a=fixture();a.ta98TermBindings.push({...a.ta98TermBindings[0]});assert.equal(api.hierarchyEntries(a).filter(e=>e.id===collection).length,1);
});
test('catalogues without collection bindings keep only their original entries',()=>{
 const a=fixture();a.ta98TermBindings=[];assert.deepEqual(api.hierarchyEntries(a).map(e=>e.id),a.parts.map(p=>p.conceptId));
});
test('logical source selectors never duplicate physical explosion ownership in any hierarchy mode',()=>{
 const a=fixture(),all=api.hierarchyEntries(a),physical=all.filter(e=>!e.selectionOnly),visible=new Set(a.parts.map(p=>p.id));
 assert.equal(all.filter(e=>e.selectionOnly).length,1);
 for(const mode of ['systems','regions','depth']){
  const actual=api.createHierarchicalExplosionLayout(a,all,visible,mode,1.7),expected=api.createHierarchicalExplosionLayout(a,physical,visible,mode,1.7);
  assert.deepEqual(actual,expected,mode);
  for(const stage of actual.stages){const members=stage.clusters.flat();assert.equal(new Set(members).size,members.length,mode);}
 }
});
const packed='.local-models/geometry-workbench/biodigital-male-circulation-v1/atlas-cardiac-source-identities-v1/atlas.json';
test('actual named source batch exposes all15 guarded family controls without new material meshes',{skip:!fs.existsSync(packed)},()=>{
 const a=JSON.parse(fs.readFileSync(packed,'utf8')),before=JSON.stringify(a),raw=api.hierarchyEntries(a),tree=api.buildSystemTree(raw,'systems',a),all=tree.flatMap(b=>entries(b.nodes));
 const bindings=a.ta98TermBindings.filter(b=>b.representation==='source-collection-identity');assert.equal(bindings.length,15);
 for(const b of bindings){const rows=all.filter(e=>e.id===b.conceptId);assert.equal(rows.length,1,b.term);assert.deepEqual(rows[0].parts.map(p=>p.id).sort(),[...b.partIds].sort(),b.term);assert.ok(rows[0].parts.every(p=>!p.suppressed&&['arterial','venous'].includes(p.system)));}
 const physical=a.parts.filter(p=>!p.suppressed);assert.equal(api.meshCount(tree.flatMap(b=>b.parts)),api.meshCount(physical));assert.equal(JSON.stringify(a),before);
});
