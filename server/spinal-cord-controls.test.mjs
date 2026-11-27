import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/spinal-cord-controls';export * from './app/spinal-cord-system';export * from './app/guest-hierarchy';export * from './app/guest-choice';export * from './app/anatomical-components';export * from './app/anatomy-hierarchy';export * from './app/anatomy-label-plan';export * from './app/viewer-state';export * from './app/section-tissue';export * from './app/tissue-wall-material';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const path=process.env.ATLAS_SPINAL_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';
const read=()=>JSON.parse(fs.readFileSync(path));
const staged=!!process.env.ATLAS_SPINAL_CANDIDATE;

test('native cord tissues have finite selectable pieces and functional peel groups',t=>{
 const a=read();if(!a.finiteSpinalIntegration&&!staged){t.skip('Cord revision is not yet installed');return;}
 const r=api.spinalCordRegistry(a);assert.ok(r);const finer=r.marker.schema.endsWith('/v4'),segmented=finer||r.marker.schema.endsWith('/v3'),regional=segmented||r.marker.schema.endsWith('/v2');assert.equal(r.parts.length,finer?r.marker.partIds.length+r.marker.coveringPartIds.length:segmented?403:regional?97:12);
 const c=a.concepts.find(c=>c.id===r.marker.wholeConceptId),byId=new Map(a.parts.map(p=>[p.id,p]));
 assert.deepEqual(new Set(api.anatomicalComponentIds(c,[],byId,a.concepts)),new Set(r.marker.partIds));
 const steps=api.spinalCordPeelSteps(a,c);assert.deepEqual(steps.map(s=>s.name),['White matter','Grey matter']);assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(r.marker.partIds));assert.equal(new Set(steps.flatMap(s=>s.ids)).size,steps.flatMap(s=>s.ids).length);
 assert.equal(api.spinalCordPeelSteps(a,{...c,elements:c.elements.slice(1)}),null);
 assert.equal(api.spinalCordPeelSteps(a,{...c,id:'selection-set'}),null);
 assert.ok(byId.get(r.marker.sourcePartId).suppressed);assert.ok(a.defaultHidden.includes(r.marker.sourcePartId));
 for(const p of r.parts){
  assert.equal(api.sectionTissue(p),'cns');if(p.tissueWallLayer!=='pia')assert.ok(p.name.includes('left')||p.name.includes('right'));
  assert.ok(a.materials[p.material].color.every(v=>v>=80&&v<=255),'Base atlas colors use 0–255 RGB');
  assert.ok(p.tissueAppearance.color.every(v=>v>0&&v<=1),'Procedural colors use normalized RGB');
  const {material}=api.tissueWallMaterial(p.tissueAppearance);assert.ok(material.color.r>0);material.dispose();
 }
});

test('one label selects only its own tissue and whole controls expose its components',t=>{
 const a=read();if(!a.finiteSpinalIntegration&&!staged){t.skip();return;}
 const r=api.spinalCordRegistry(a),c=a.concepts.find(c=>c.id===r.marker.wholeConceptId);
 const labels=api.selectedAnatomyLabels(a,c.elements,c,true);
 assert.ok(labels.length>=12&&labels.length<=new Set(c.elements.map(id=>a.parts.find(p=>p.id===id).conceptId)).size);
 const finer=r.marker.schema.endsWith('/v4'),segmented=finer||r.marker.schema.endsWith('/v3'),regional=segmented||r.marker.schema.endsWith('/v2');
 for(const label of labels){const choice=api.anatomyLabelChoice(label);assert.ok(finer?choice.elements.length>0:choice.elements.length===(segmented?31:regional?8:1));assert.ok(choice.elements.every(id=>c.elements.includes(id)));}
 const terms=api.spinalCordTermParts(a);assert.equal(terms.size,finer?r.marker.termGroups.length:regional?28:16);assert.equal(terms.get('A14.1.02.020').length,finer?r.parts.filter(p=>p.tissueWallLayer==='grey').length:segmented?186:regional?48:6);assert.equal(terms.get('A14.1.02.239').length,finer?r.marker.cellRegions.filter(r=>r.nativeRegion.startsWith('white-posterior-')).length:segmented?62:regional?16:2);
 if(regional){
  const regions=['014','015','016','017','018'].map(code=>terms.get('A14.1.02.'+code));
  assert.equal(new Set(regions.flat().map(p=>p.id)).size,finer?r.marker.partIds.length:segmented?372:96);assert.equal(regions.flat().length,finer?r.marker.partIds.length:segmented?372:96);
  assert.equal(terms.get('A14.1.01.309').length,segmented?31:1);
 }
});

test('a wrong frame, broken reciprocal boundary or incomplete concept rejects the controls',t=>{
 if(!read().finiteSpinalIntegration&&!staged){t.skip();return;}
 for(const change of [
  a=>a.finiteSpinalIntegration.sex=a.sex==='male'?'female':'male',
  a=>a.parts.find(p=>p.id===a.finiteSpinalIntegration.partIds[0]).tissueDrawGroups[0].count+=3,
  a=>a.concepts.find(c=>c.id===a.finiteSpinalIntegration.wholeConceptId).elements.pop(),
 ]){const a=read();change(a);assert.equal(api.spinalCordRegistry(a),null);}
});

test('regional controls resolve in TA98 and Systems with reciprocal pial covering',t=>{
 const a=read(),r=api.spinalCordRegistry(a);if(!r||r.marker.schema.endsWith('/v1')){t.skip();return;}
 const finer=r.marker.schema.endsWith('/v4'),segmented=finer||r.marker.schema.endsWith('/v3');
 const entries=api.hierarchyEntries(a);
 for(const part of r.parts){const entry=entries.find(e=>e.parts.some(p=>p.id===part.id));assert.ok(entry,part.id);assert.ok(api.systemPathFor(entry).includes('Spinal cord'),part.id);}
 const roots=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))));
 const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};roots.forEach(walk);
 for(const group of r.marker.termGroups){const node=nodes.get(group.term);assert.ok(node,group.term);assert.deepEqual(new Set(node.directParts.map(p=>p.id)),new Set(group.partIds),group.term);}
 const covers=api.spinalCordCoveringRelations(a.parts),pia=r.marker.coveringPartIds[0];
 assert.equal(covers.get(pia).length,finer?r.marker.segments.find(s=>s.coveringPartIds.includes(pia)).partIds.length:segmented?12:96);
 for(const id of r.marker.partIds)assert.deepEqual(covers.get(id).map(r=>[r.kind,r.target.id]),[['coveredBy',segmented?r.marker.segments.find(s=>s.name===id.split('@')[1]).coveringPartIds[0]:pia]]);
 const before=api.buildAnatomyNodes(entries,e=>api.systemPathFor(e)),system=api.withSpinalCordSystem(a,before);
 assert.deepEqual(new Set(system.flatMap(n=>n.parts.map(p=>p.id))),new Set(before.flatMap(n=>n.parts.map(p=>p.id))),'Systems retains all existing native structures');
 const all=[];const descend=n=>{all.push(n);if(n.kind==='group')n.nodes.forEach(descend);};system.forEach(descend);
 const whole=all.find(n=>n.id==='spinal:A14.1.02.001');assert.equal(whole?.itemCount,119);
 assert.ok(finer?all.find(n=>n.id==='spinal:A14.1.02.121')?.parts.length>0:all.find(n=>n.id==='spinal:A14.1.02.121')?.parts.length===0,'Only constructed nuclei receive controls');
 for(const group of r.marker.termGroups){
  const node=all.find(n=>n.id==='spinal:'+group.term);assert.ok(node,group.term);
  assert.ok(group.partIds.every(id=>node.parts.some(p=>p.id===id)),group.term);
 }
 assert.deepEqual(api.spinalCordPeelSteps(a,{id:'systems:cord',name:'Spinal cord',elements:whole.parts.map(p=>p.id)}).map(s=>s.name),['Pia mater','White matter','Grey matter']);
});

test('all 31 segment controls select disjoint tissue in Systems and TA98 without inflating term coverage',t=>{
 const a=read(),r=api.spinalCordRegistry(a);if(!r||!['spinal-cord-native-core/v3','spinal-cord-native-core/v4'].includes(r.marker.schema)){t.skip();return;}const finer=r.marker.schema.endsWith('/v4');
 const segments=api.spinalCordSegments(a);assert.equal(segments.length,31);
 assert.equal(a.finiteSpinalIntegration.representedCordTerms,r.marker.termGroups.filter(g=>g.term.startsWith('A14.1.02.')).length);assert.ok(finer?a.finiteSpinalIntegration.representedCordTerms>=61:a.finiteSpinalIntegration.representedCordTerms===24);
 assert.equal(new Set(segments.flatMap(s=>s.partIds)).size,r.marker.partIds.length);
 const guest=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))));
 const guestNodes=new Map();const collect=n=>{guestNodes.set(n.id,n);n.children.forEach(collect);};guest.forEach(collect);
 const entries=api.hierarchyEntries(a),system=api.withSpinalCordSystem(a,api.buildAnatomyNodes(entries,e=>api.systemPathFor(e)));
 const systemNodes=new Map();const visit=n=>{systemNodes.set(n.id,n);if(n.kind==='group')n.nodes.forEach(visit);};system.forEach(visit);
 for(const s of segments){
  assert.ok(finer?s.partIds.length>=12:s.partIds.length===12);assert.ok(s.partIds.every(id=>id.endsWith('@'+s.name)));
  const c=a.concepts.find(c=>c.id===s.conceptId);assert.equal(c.ta98Term,undefined);assert.equal(c.partOf,r.marker.termGroups.find(g=>g.term===s.parentTerm).conceptId);
  assert.deepEqual(new Set(guestNodes.get(s.conceptId).directParts.map(p=>p.id)),new Set(s.partIds));
  assert.deepEqual(new Set(systemNodes.get(s.conceptId).parts.map(p=>p.id)),new Set(s.partIds));
  const labels=api.selectedAnatomyLabels(a,s.partIds,c,true);
  assert.deepEqual(labels.map(label=>label.name),[`White substance · ${s.name}`,`Grey substance · ${s.name}`]);
  assert.deepEqual(new Set(labels.flatMap(label=>label.members)),new Set(s.partIds));
 }
 for(const change of [
  m=>m.segments[0].partIds[0]=m.segments[1].partIds[0],
  m=>m.segments[0].coveringPartIds[0]=m.segments[1].coveringPartIds[0],
  m=>m.segments[1].yRange[1]+=.001,
 ]){const bad=read();change(bad.finiteSpinalIntegration);assert.equal(api.spinalCordRegistry(bad),null);}
});

test('the canal is a selectable non-owning boundary, never a solid tissue or a duplicate whole-organ draw',t=>{
 const a=read(),r=api.spinalCordRegistry(a);if(!r?.referenceParts.length){t.skip('No canal boundary revision in this fixture');return;}
 const ids=r.marker.referencePartIds,choice=a.concepts.find(c=>c.ta98Term==='A14.1.02.019'&&c.elements.every(id=>ids.includes(id)));
 assert.ok(choice);assert.deepEqual(new Set(choice.elements),new Set(ids));
 const physical=new Set([...r.marker.partIds,...r.marker.coveringPartIds]);
 assert.ok(ids.every(id=>!physical.has(id)));
 const state={visible:['nervous'],selected:ids,hidden:[],depthHidden:[],region:'all',isolate:true,contextOpacity:0,skinOpacity:0};
 for(const p of r.referenceParts){
  assert.equal(p.sectionSolid,false);assert.equal(p.pairedTissue,false);
  assert.equal(api.partVisible(p,state,a.parts),true);assert.equal(api.sectionPartVisible(p,state),true);
  const together={...state,selected:[...ids,...physical]};
  assert.equal(api.partVisible(p,together,a.parts),false);assert.equal(api.sectionPartVisible(p,together),false);
  assert.equal(api.partVisible(p,{...state,selected:[]},a.parts),false);
 }
 assert.equal(api.spinalCordPeelSteps(a,choice),null);
 for(const corrupt of [
  p=>p.spinalBoundaryReference.sourceIndexStart+=1,
  p=>p.spinalBoundaryReference.sourceBankSha256='0'.repeat(64),
  p=>p.indices+=4,
  p=>p.sectionSolid=true,
 ]){const bad=read(),p=bad.parts.find(p=>p.id===ids[0]);corrupt(p);assert.equal(api.spinalCordRegistry(bad),null);}
});

test('regional schema requires exhaustive native tissue ownership at every segment',t=>{
 const fixture=()=>{
  const a=read(),m=a.finiteSpinalIntegration;
  if(!m?.segments)return null;
  if(!m.schema.endsWith('/v4')){
   m.schema='spinal-cord-native-core/v4';
   m.cellRegions=[...m.partIds,...m.coveringPartIds].map(partId=>({partId,segment:partId.split('@')[1],nativeRegion:partId.replace(`TA98CORD:${a.sex}:`,'').split('@')[0]}));
  }
  return a;
 };
 const good=fixture();if(!good){t.skip();return;}
 assert.ok(api.spinalCordRegistry(good));
 for(const corrupt of [
  m=>m.cellRegions.pop(),
  m=>m.cellRegions[0].segment='Co1',
  m=>m.cellRegions[0].nativeRegion='grey-unknown-left',
  m=>m.cellRegions[1].partId=m.cellRegions[0].partId,
 ]){const bad=fixture();corrupt(bad.finiteSpinalIntegration);assert.equal(api.spinalCordRegistry(bad),null);}
});

test('the expanded spinal hierarchy contains named terms and segments without repeated internal tissue rows',t=>{
 const a=read(),r=api.spinalCordRegistry(a);if(!r){t.skip();return;}
 const raw=JSON.stringify(a),roots=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))));
 const byId=new Map();const index=n=>{byId.set(n.id,n);n.children.forEach(index);};roots.forEach(index);
 const root=byId.get('A14.1.02.001');assert.ok(root);assert.equal(root.itemCount,119);
 let named=0,sourceRows=0;const count=n=>{named++;sourceRows+=new Set(api.guestDirectAnatomyParts('ta98',n).map(p=>p.conceptId)).size;n.children.forEach(count);};count(root);
 assert.equal(named,119+r.marker.segments.length);assert.equal(sourceRows,0,'Storage pieces must not become additional anatomical rows');
 let choices=0;const visit=c=>{choices++;c.children?.forEach(visit);};visit(api.guestNodeChoice('ta98',root));assert.equal(choices,named);
 for(const g of r.marker.termGroups){
  const node=byId.get(g.term),choice=api.guestNodeChoice('ta98',node);
  assert.deepEqual(new Set(node.directParts.map(p=>p.id)),new Set(g.partIds));
  assert.deepEqual(new Set(choice.elements),new Set(node.parts.map(p=>p.id)),'Named selection retains every physical mesh piece');
 }
 for(const segment of r.marker.segments){const node=byId.get(segment.conceptId);assert.equal(node.extension,undefined);assert.equal(node.itemCount,0);assert.deepEqual(new Set(api.guestNodeChoice('ta98',node).elements),new Set(segment.partIds));}
 assert.equal(JSON.stringify(a),raw,'Presentation changes neither the installed catalogue nor its geometry');
});
