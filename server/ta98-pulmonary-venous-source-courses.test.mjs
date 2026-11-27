import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const root='.local-models/geometry-workbench/biodigital-male-circulation-v1';
const base=JSON.parse(fs.readFileSync(root+'/atlas-cardiovascular-review-v15/atlas.json','utf8'));
const candidate=JSON.parse(fs.readFileSync(root+'/source-course-identity-pv-v1/candidate-v1/atlas.json','utf8'));
const delta=JSON.parse(fs.readFileSync(root+'/source-course-identity-pv-v1/delta.json','utf8'));
const markers=delta.appendSourceCourseIdentities,parent=markers[0].taParentConceptId;
const built=await build({stdin:{contents:"export * from './app/ta98-term-bindings';export * from './app/anatomy-hierarchy';export * from './app/system-tree';export * from './app/guest-hierarchy';export * from './app/guest-choice';export * from './app/hierarchy-choice';export * from './app/anatomy-label-plan';export * from './app/organ-models';export * from './app/hierarchical-explosion';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const flatten=nodes=>nodes.flatMap(n=>n.kind==='entry'?[n]:n.kind==='bilateral'?n.entries.map(e=>({kind:'entry',entry:e,parts:e.parts,name:e.name})):flatten(n.nodes));
const authority=JSON.parse(fs.readFileSync('public/assets/ta98.json','utf8'));
const graph=api.parseGuestHierarchy({schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98',roots:['A12.3.02.111'],nodes:authority.nodes.filter(n=>['A12.3.02.111',...markers.map(m=>m.term)].includes(n.id))});
const fixture=()=>{
 const ids=new Set([parent,...markers.map(m=>m.sourcePartId)]),cids=new Set([...ids,...markers.map(m=>m.id),markers[0].sourceClinicalParentConceptId]);
 return structuredClone({...candidate,parts:candidate.parts.filter(p=>ids.has(p.id)),concepts:candidate.concepts.filter(c=>cids.has(c.id)),defaultHidden:[],ta98TermBindings:candidate.ta98TermBindings.filter(b=>['A12.3.02.111',...markers.map(m=>m.term)].includes(b.term)),sourceCourseIdentities:markers});
};

test('two venous course selectors preserve all v15 material meshes, banks, poses and common middle-lobe stem',()=>{
 for(const field of ['parts','chunks','materials','defaultHidden'])assert.deepEqual(candidate[field],base[field],field);
 const map=api.catalogueTa98TermParts(candidate),old=api.catalogueTa98TermParts(base);
 assert.deepEqual(map.get('A12.3.02.111').map(p=>p.id),[parent]);assert.deepEqual(map.get('A12.3.02.111'),old.get('A12.3.02.111'));
 for(const m of markers){assert.equal(old.has(m.term),false);assert.deepEqual(map.get(m.term).map(p=>p.id),[m.sourcePartId]);assert.equal(map.get(m.term)[0].ta98Term,undefined);assert.equal(map.get(m.term)[0].conceptId,m.sourcePartId);}
 for(const missing of ['A12.3.02.103','A12.3.02.104','A12.3.02.106','A12.3.02.107','A12.3.02.109','A12.3.02.110'])assert.equal(map.has(missing),false,missing);
 const all=api.hierarchyEntries(candidate),tree=api.buildSystemTree(all,'systems',candidate),rows=tree.flatMap(b=>flatten(b.nodes));
 for(const m of markers){const [entry]=rows.filter(n=>n.entry.id===m.id);assert.ok(entry);assert.equal(rows.filter(n=>n.entry.id===m.id).length,1);assert.equal(entry.entry.selectionOnly,true);assert.equal(entry.entry.taParentConcept,parent);assert.equal(entry.entry.ta98Term,m.term);assert.deepEqual(entry.parts.map(p=>p.id),[m.sourcePartId]);}
 assert.equal(api.meshCount(tree.flatMap(b=>b.parts)),api.meshCount(candidate.parts.filter(p=>!p.suppressed)));
});

test('actual TA98 parent selects three unique physical sources; normal and TA98 labels select each requested course',()=>{
 const a=fixture(),[resolved]=api.resolveGuestHierarchy(a,graph),expected=[parent,...markers.map(m=>m.sourcePartId)];
 assert.deepEqual(new Set(resolved.parts.map(p=>p.id)),new Set(expected));assert.equal(api.meshCount(resolved.parts),3);assert.deepEqual(resolved.directParts.map(p=>p.id),[parent]);
 assert.equal(a.parts.length,3);assert.equal(new Set(a.parts.map(p=>p.id)).size,3);
 for(const m of markers){
  const node=resolved.children.find(n=>n.id===m.term),logical=a.concepts.find(c=>c.id===m.id);assert.deepEqual(node.directParts.map(p=>p.id),[m.sourcePartId]);assert.match(node.description,/common middle-lobe stem is separate/);assert.doesNotMatch(node.description,/Stored finite tissue|No direct reviewed mesh correspondence/);
  const direct=api.selectedAnatomyLabels(a,[m.sourcePartId],logical,false);assert.equal(direct.length,1);assert.equal(direct[0].name,logical.name);assert.deepEqual(direct[0].members,[m.sourcePartId]);assert.deepEqual(direct[0].components,[]);
  const choice=api.guestNodeChoice('ta98',node),labels=api.selectedAnatomyLabels(a,[m.sourcePartId],choice,false);assert.equal(labels.length,1);assert.equal(labels[0].name,node.name);assert.deepEqual(labels[0].members,[m.sourcePartId]);
 }
 assert.equal(a.concepts.find(c=>c.id===parent).ta98Term,'A12.3.02.111');
});

test('source-course selectors add no physical explosion ownership or duplicate source meshes',()=>{
 const a=fixture(),entries=api.hierarchyEntries(a),physical=entries.filter(e=>!e.selectionOnly),visible=new Set(a.parts.map(p=>p.id));
 assert.equal(entries.filter(e=>e.selectionOnly).length,2);
 for(const mode of ['systems','regions','depth'])assert.deepEqual(api.createHierarchicalExplosionLayout(a,entries,visible,mode,1.7),api.createHierarchicalExplosionLayout(a,physical,visible,mode,1.7));
});

test('both course guards reject forged source, pose, parent, bank, dimension, tissue, patency and crosswalk evidence',()=>{
 for(const wanted of markers){
  const edits=[
   (a,p,m)=>a.sex='female',(a,p,m)=>delete a.biodigitalReview,(a,p,m)=>a.biodigitalReview.optInOnly=false,(a,p,m)=>a.biodigitalReview.installedAsDefault=true,(a,p,m)=>a.biodigitalReview.anatomicalAcceptance=true,
   (a,p,m)=>p.sourceId='wrong source',(a,p,m)=>p.conceptId=parent,(a,p,m)=>p.sourceGroups=[],(a,p,m)=>p.positions++,(a,p,m)=>p.normalSigns++,(a,p,m)=>p.indices++,(a,p,m)=>p.vertexCount++,(a,p,m)=>p.indexCount-=3,(a,p,m)=>p.ta98Term=m.term,
   (a,p,m)=>p.geometryOrigin[0]+=.001,(a,p,m)=>p.geometryScale=1.01,(a,p,m)=>p.sourceOffset=[.001,0,0],(a,p,m)=>p.geometryInstance={scale:1},(a,p,m)=>p.bounds[0][0]-=.001,
   (a,p,m)=>p.sectionSolid=true,(a,p,m)=>p.tissueWallLayer='venous-wall',(a,p,m)=>p.description='complete tissue wall',(a,p,m)=>p.suppressed=true,(a,p,m)=>a.defaultHidden=[p.id],
   (a,p,m)=>a.chunks[p.chunk].sha256='wrong bank',(a,p,m)=>a.chunks[p.chunk].url='other.bin',(a,p,m)=>m.logicalCarrierPose.geometryOrigin[0]+=.001,(a,p,m)=>m.dimension='point',(a,p,m)=>m.finiteWallEstablished=true,(a,p,m)=>m.lumenPatencyEstablished=true,(a,p,m)=>m.anatomicalAcceptance=true,
   (a,p,m)=>m.taParentTerm='A12.3.02.101',(a,p,m)=>m.taParentConceptId=m.sourcePartId,(a,p,m)=>m.sourceClinicalParentConceptId=parent,(a,p,m)=>m.evidence.bronchialSourceId='wrong bronchus',(a,p,m)=>m.evidence.crosswalkSha256='wrong crosswalk',(a,p,m)=>m.evidence.sourceFaces--,(a,p,m)=>m.evidence.proximalSourceMouthGapSourceUnits=0,
   (a,p,m)=>a.concepts.find(c=>c.id===m.id).partOf=m.sourcePartId,(a,p,m)=>a.concepts.find(c=>c.id===m.id).description='accepted tissue',(a,p,m)=>a.concepts.find(c=>c.id===parent).ta98Term='A12.3.02.101',(a,p,m)=>a.concepts.find(c=>c.id===parent).partOf='wrong source parent',(a,p,m)=>a.concepts.find(c=>c.id===p.conceptId).partOf=parent,
   (a,p,m)=>a.sourceCourseIdentities.push(structuredClone(m)),(a,p,m)=>a.ta98TermBindings.find(b=>b.term===m.term).representation='source-collection-identity',
  ];
  for(const edit of edits){const a=fixture(),p=a.parts.find(p=>p.id===wanted.sourcePartId),m=a.sourceCourseIdentities.find(m=>m.term===wanted.term);edit(a,p,m);assert.equal(api.catalogueTa98TermParts(a).has(wanted.term),false,wanted.term+' '+edit);assert.equal(api.hierarchyEntries(a).some(e=>e.id===wanted.id),false,wanted.term+' '+edit);assert.equal(api.resolveGuestHierarchy(a,graph)[0].children.find(n=>n.id===wanted.term).directParts.length,0,wanted.term+' '+edit);}
 }
});
