import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/source-course-identity-v1';
const candidate=JSON.parse(fs.readFileSync(base+'/candidate-v1/atlas.json','utf8'));
const term='A12.2.01.108',cid='TA98COURSE:male:pulmonary-right-posterior-ascending';
const marker=candidate.sourceCourseIdentities[0],pid=marker.sourcePartId;
const built=await build({stdin:{contents:"export * from './app/ta98-term-bindings';export * from './app/anatomy-hierarchy';export * from './app/system-tree';export * from './app/guest-hierarchy';export * from './app/guest-choice';export * from './app/hierarchy-choice';export * from './app/hierarchy-visibility';export * from './app/anatomy-label-plan';export * from './app/organ-models';export * from './app/hierarchical-explosion';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const fixture=()=>structuredClone({sex:candidate.sex,source:candidate.source,scope:candidate.scope,parts:[candidate.parts.find(p=>p.id===pid)],concepts:candidate.concepts.filter(c=>[cid,pid,marker.sourceClinicalParentConceptId].includes(c.id)),chunks:candidate.chunks,materials:candidate.materials,defaultHidden:[],biodigitalReview:{optInOnly:true,installedAsDefault:false,anatomicalAcceptance:false,sourceBaseCatalogueSha256:candidate.biodigitalReview.sourceBaseCatalogueSha256,sourceMasterSha256:candidate.biodigitalReview.sourceMasterSha256},sourceCourseIdentities:[marker],ta98TermBindings:candidate.ta98TermBindings.filter(b=>[term,'A12.2.01.107'].includes(b.term))});
const flatten=nodes=>nodes.flatMap(n=>n.kind==='entry'?[n]:n.kind==='bilateral'?n.entries.map(e=>({kind:'entry',entry:e,parts:e.parts,name:e.name})):flatten(n.nodes));
const guest=()=>api.parseGuestHierarchy({schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98',roots:['A12.2.01.107'],nodes:[{id:'A12.2.01.107',name:'posterior segmental artery',childIds:[term],matches:[]},{id:term,name:'ascending branch',matches:[],description:'TA98: '+term+'\nNo direct reviewed mesh correspondence.'}]});

test('the actual v14 candidate binds one singular whole course without duplicating or replacing its parent source',()=>{
 const map=api.catalogueTa98TermParts(candidate);
 assert.deepEqual(map.get(term).map(p=>p.id),[pid]);assert.deepEqual(map.get('A12.2.01.107').map(p=>p.id),[pid]);
 for(const missing of ['A12.2.01.105','A12.2.01.106','A12.2.01.109'])assert.equal(map.has(missing),false,missing);
 const old=JSON.parse(fs.readFileSync('.local-models/geometry-workbench/biodigital-male-circulation-v1/atlas-cardiovascular-review-v14/atlas.json','utf8'));
 assert.deepEqual(candidate.parts,old.parts);assert.deepEqual(candidate.chunks,old.chunks);assert.deepEqual(candidate.materials,old.materials);
 assert.deepEqual(map.get('A12.2.01.104').map(p=>p.id),api.catalogueTa98TermParts(old).get('A12.2.01.104').map(p=>p.id));
 assert.equal(candidate.sourceCourseIdentities.length,1);assert.equal(candidate.sourceCourseIdentities[0].dimension,'source-surface-envelope');
 const all=api.hierarchyEntries(candidate),tree=api.buildSystemTree(all,'systems',candidate),rows=tree.flatMap(b=>flatten(b.nodes));
 assert.equal(rows.filter(n=>n.entry.id===cid).length,1);assert.equal(rows.find(n=>n.entry.id===cid).entry.selectionOnly,true);
 assert.equal(api.meshCount(tree.flatMap(b=>b.parts)),api.meshCount(candidate.parts.filter(p=>!p.suppressed)));
 const authority=JSON.parse(fs.readFileSync('public/assets/ta98.json','utf8')),parent='A12.2.01.107';
 const graph=api.parseGuestHierarchy({schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98',roots:[parent],nodes:authority.nodes.filter(n=>[parent,term,'A12.2.01.109'].includes(n.id))});
 const [actualParent]=api.resolveGuestHierarchy(candidate,graph),actualCourse=actualParent.children.find(n=>n.id===term);
 assert.deepEqual(actualParent.parts.map(p=>p.id),[pid]);assert.deepEqual(actualCourse.directParts.map(p=>p.id),[pid]);assert.equal(api.meshCount(actualParent.parts),1);
 const choice=api.guestNodeChoice('ta98',actualCourse),labels=api.selectedAnatomyLabels(candidate,[pid],choice,false);
 assert.equal(labels.length,1);assert.equal(labels[0].name,actualCourse.name);assert.deepEqual(labels[0].members,[pid]);
});

test('Systems and TA98 select the existing carrier with its true parent and a qualified course label',()=>{
 const a=fixture(),before=structuredClone(a),entries=api.hierarchyEntries(a),course=entries.find(e=>e.id===cid),parent=entries.find(e=>e.id===pid);
 assert.equal(course.selectionOnly,true);assert.equal(course.taParentConcept,pid);assert.equal(course.ta98Term,term);assert.equal(parent.ta98Term,'A12.2.01.107');assert.deepEqual(course.parts,parent.parts);
 const tree=api.buildSystemTree(entries,'systems',a),branch=tree.find(t=>t.id==='circulatory'),node=flatten(branch.nodes).find(n=>n.entry.id===cid);
 assert.ok(node);assert.match(node.entry.name,/Ascending branch of posterior segmental artery \(right\)/);assert.equal(api.meshCount(branch.parts),1);
 const choice=api.anatomyNodeChoice(node,'course');assert.deepEqual(choice.elements,[pid]);assert.equal(choice.ta98Term,term);
 const [resolved]=api.resolveGuestHierarchy(a,guest());assert.deepEqual(resolved.parts.map(p=>p.id),[pid]);assert.deepEqual(resolved.children[0].directParts.map(p=>p.id),[pid]);assert.equal(api.meshCount(resolved.parts),1);
 assert.match(resolved.children[0].description,/existing whole source envelope/);assert.doesNotMatch(resolved.children[0].description,/Stored finite tissue|No direct reviewed mesh correspondence/);
 const selected=api.guestNodeChoice('ta98',resolved.children[0]);assert.deepEqual(selected.elements,[pid]);
 const state={visible:['arterial'],hidden:[],selected:[],isolate:false,camera:{position:[1,2,3]},cameraIntent:13};
 const changed=api.toggleSystemParts(state,node.parts,a.parts);assert.equal(changed.camera,state.camera);assert.equal(changed.cameraIntent,13);assert.deepEqual(changed.hidden,[pid]);assert.deepEqual(a,before);
});

test('normal selected-anatomy callouts retain the requested singular course identity',()=>{
 const a=fixture(),c=a.concepts.find(c=>c.id===cid),labels=api.selectedAnatomyLabels(a,[pid],c,false);
 assert.equal(labels.length,1);assert.equal(labels[0].id,cid);assert.equal(labels[0].name,c.name);assert.deepEqual(labels[0].members,[pid]);assert.deepEqual(labels[0].components,[]);
 const [parent]=api.resolveGuestHierarchy(a,guest()),choice=api.guestNodeChoice('ta98',parent.children[0]);
 assert.equal(api.selectedAnatomyLabels(a,[pid],choice,false)[0].name,'ascending branch');
});

test('singular logical course does not add physical explosion ownership',()=>{
 const a=fixture(),all=api.hierarchyEntries(a),physical=all.filter(e=>!e.selectionOnly),visible=new Set([pid]);assert.equal(all.filter(e=>e.selectionOnly).length,1);
 for(const mode of ['systems','regions','depth'])assert.deepEqual(api.createHierarchicalExplosionLayout(a,all,visible,mode,1.7),api.createHierarchicalExplosionLayout(a,physical,visible,mode,1.7));
});

test('frozen source, bank, logical pose, anatomical parent, dimension and tissue claims reject stale or forged bindings',()=>{
 const mutate=[
  a=>a.sex='female',a=>delete a.biodigitalReview,a=>a.biodigitalReview.optInOnly=false,a=>a.biodigitalReview.installedAsDefault=true,a=>a.biodigitalReview.anatomicalAcceptance=true,a=>a.biodigitalReview.sourceMasterSha256='forged',a=>a.biodigitalReview.sourceBaseCatalogueSha256='forged',
  a=>a.parts[0].sourceId='forged',a=>a.parts[0].conceptId='forged',a=>a.parts[0].ta98Term='A12.2.01.104',a=>a.parts[0].system='venous',a=>a.parts[0].positions++,a=>a.parts[0].normalSigns++,a=>a.parts[0].indexCount-=3,a=>a.parts[0].vertexCount++,a=>a.parts[0].sourceGroups=[],
  a=>a.parts[0].geometryOrigin[0]+=.001,a=>a.parts[0].geometryScale=1.1,a=>a.parts[0].sourceOffset=[.001,0,0],a=>a.parts[0].geometryInstance={rotation:[0,0,0,1],scale:1},a=>a.parts[0].bounds[0][0]-=.001,
  a=>a.parts[0].sectionSolid=true,a=>a.parts[0].pairedTissue=true,a=>a.parts[0].tissueWallAssembly='invented',a=>a.parts[0].tissueWallLayer='arterial-wall',a=>a.parts[0].sectionAssembly='invented-wall',a=>a.parts[0].tissue='arterial-media',a=>a.parts[0].description='Tissue wall complete',a=>a.parts[0].suppressed=true,a=>a.defaultHidden=[pid],
  a=>{a.parts[0].conceptId=cid;a.parts[0].ta98Term=term;},
  a=>a.chunks[a.parts[0].chunk].sha256='forged',a=>a.chunks[a.parts[0].chunk].url='another.bin',a=>a.chunks[a.parts[0].chunk].bytes++,
  a=>a.sourceCourseIdentities[0].dimension='point',a=>a.sourceCourseIdentities[0].finiteWallEstablished=true,a=>a.sourceCourseIdentities[0].lumenPatencyEstablished=true,a=>a.sourceCourseIdentities[0].anatomicalAcceptance=true,a=>a.sourceCourseIdentities[0].evidence.bronchialSourceId='wrong bronchus',a=>a.sourceCourseIdentities[0].evidence.crosswalkSha256='forged',a=>a.sourceCourseIdentities[0].evidence.sourceFaces--,
  a=>a.sourceCourseIdentities[0].logicalCarrierPose.geometryOrigin[0]+=.001,a=>a.sourceCourseIdentities[0].taParentTerm='A12.2.01.104',a=>a.sourceCourseIdentities[0].sourceClinicalParentConceptId='wrong parent',
  a=>a.concepts.find(c=>c.id===cid).partOf='wrong parent',a=>a.concepts.find(c=>c.id===cid).coverageLimitation='Finite wall accepted',a=>a.concepts.find(c=>c.id===cid).description='Tissue complete',a=>a.concepts.find(c=>c.id===pid).partOf='wrong parent',a=>a.concepts.find(c=>c.id===pid).ta98Term='A12.2.01.104',
  a=>a.ta98TermBindings.find(b=>b.term===term).representation='source-collection-identity',a=>a.sourceCourseIdentities.push(structuredClone(a.sourceCourseIdentities[0])),
 ];
 for(const edit of mutate){const a=fixture();edit(a);assert.equal(api.catalogueTa98TermParts(a).has(term),false,String(edit));assert.equal(api.hierarchyEntries(a).some(e=>e.id===cid),false,String(edit));assert.equal(api.resolveGuestHierarchy(a,guest())[0].children[0].directParts.length,0,String(edit));assert.equal(api.selectedAnatomyLabels(a,[pid],a.concepts.find(c=>c.id===cid),false).some(l=>l.id===cid),false,String(edit));}
});

test('unbound/default catalogues and unsupported feeder terms acquire no course control',()=>{
 const a=fixture();delete a.sourceCourseIdentities;assert.equal(api.catalogueTa98TermParts(a).has(term),false);assert.equal(api.hierarchyEntries(a).some(e=>e.id===cid),false);
 for(const code of ['A12.2.01.105','A12.2.01.106','A12.2.01.109']){const a=fixture();a.ta98TermBindings.find(b=>b.term===term).term=code;assert.equal(api.catalogueTa98TermParts(a).has(code),false);}
});
