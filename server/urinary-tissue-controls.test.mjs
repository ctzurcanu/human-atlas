import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/urinary-wall';export * from './app/renal-wall';export * from './app/bowel-wall-peeling';export * from './app/guest-hierarchy';export * from './app/anatomical-components';export * from './app/urinary-relations';export * from './app/anatomical-relations';export * from './app/laterality';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const path=process.env.ATLAS_URINARY_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
const read=()=>JSON.parse(fs.readFileSync(path));
const hasFeatures=a=>['urinary-local-coats-v2','urinary-local-coats-v3','urinary-local-coats-v4'].includes(a.finiteUrinaryIntegration.schema);
const hasOpenings=a=>['urinary-local-coats-v3','urinary-local-coats-v4'].includes(a.finiteUrinaryIntegration.schema);
const bladderCount=a=>a.parts.filter(p=>!p.suppressed&&p.urinaryOrgan==='bladder').length;

test('Male urinary coats retain both complete kidneys and explicit cross-chart joins',()=>{
 const a=read(),r=api.urinaryWallRegistry(a);assert.ok(r);assert.equal(r.parts.length,(r.marker.schema==='urinary-local-coats-v4'?88:76)+bladderCount(a));
 assert.equal(api.renalWallRegistry(a)?.marker.partIds.length,r.marker.schema==='urinary-local-coats-v4'?60:66);
 if(r.marker.schema==='urinary-local-coats-v4')for(const side of ['left','right'])for(const region of ['abdominal','pelvic','intramural'])for(const coat of ['adventitia','muscularis','mucosa'])assert.ok(r.byId.has(`URETER:male:${side}:${region}:${coat}`));
 assert.ok(r.marker.externalInterfaces.length>=13);
 for(const seam of r.marker.externalInterfaces){
  const p=r.byId.get(seam.first);assert.equal(api.urinaryExternalPeer(a,p,seam.second)?.id,seam.second);
 }
});

test('Whole bladder peels every stored tissue exactly once, while partial coats do not inherit that control',()=>{
 const a=read(),c=a.concepts.find(c=>c.id==='ZA:Urinary bladder'),r=api.urinaryWallRegistry(a);
 assert.equal(c.elements.length,bladderCount(a));const steps=api.bowelWallPeelSteps(a,c);assert.ok(steps);assert.equal(steps.length,a.finiteUrinaryIntegration.apicalAttachment?8:7);
 assert.deepEqual(new Set(c.elements),new Set(a.parts.filter(p=>!p.suppressed&&p.urinaryOrgan==='bladder').map(p=>p.id)));
 assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(c.elements));assert.equal(steps.flatMap(s=>s.ids).length,bladderCount(a));
 assert.equal(api.bowelWallPeelSteps(a,{...c,elements:c.elements.slice(1)}),null);
 assert.equal(api.bowelWallPeelSteps(a,{...c,id:'selection-set'}),null);
 assert.deepEqual(new Set(api.anatomicalComponentIds(c,[],new Map(a.parts.map(p=>[p.id,p])),a.concepts)),new Set(c.elements));
 assert.ok(c.elements.every(id=>r.byId.get(id)?.sectionSolid));
});

test('Every constructed bladder term selects its actual stored cells in TA98',()=>{
 const a=read(),r=api.urinaryWallRegistry(a),roots=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))));
 const nodes=new Map(),todo=[...roots];while(todo.length){const n=todo.pop();if(nodes.has(n.id))continue;nodes.set(n.id,n);todo.push(...n.children);}
 for(const g of r.marker.termGroups)assert.deepEqual(new Set(nodes.get(g.term)?.directParts.map(p=>p.id)),new Set(g.partIds),g.term);
 for(const term of [...(r.marker.termGroups.some(g=>g.term==='A08.3.01.002'&&g.pointReference)?[]:['A08.3.01.002']),...(a.finiteUrinaryIntegration.apicalAttachment?[]:['A08.3.01.007']),...(hasFeatures(a)?[]:['A08.3.01.006','A08.3.01.025'])]){
  assert.ok(!(nodes.get(term)?.directParts??[]).some(p=>p.id.startsWith('BLADDER:')),'An unfinished landmark must not become a generic wall cell: '+term);
 }
 if(hasFeatures(a))for(const [term,role]of [['A08.3.01.006','mucosa-uvula'],['A08.3.01.025','mucosa-interureteric-crest']]){
  const feature=nodes.get(term)?.directParts??[];assert.ok(feature.length);
  assert.ok(feature.every(p=>p.tissueWallLayer===role));
  assert.ok(feature.every(p=>r.marker.termGroups.find(g=>g.term==='A08.3.01.023').partIds.includes(p.id)));
 }
});

test('The complete bladder keeps its ordinary urine route; a single coat is not a whole organ',()=>{
 const a=read(),c=a.concepts.find(c=>c.id==='ZA:Urinary bladder'),members=a.parts.filter(p=>c.elements.includes(p.id));
 assert.equal(api.urinaryWholeSource(members,a.parts)?.id,'ZA:Urinary bladder');
 assert.equal(api.urinaryWholeSource(members.slice(1),a.parts),undefined);
 const links=api.nativeUrinaryRelations(a.parts).get('ZA:Urinary bladder');assert.ok(links?.some(r=>r.kind==='before'));
 assert.ok(links?.some(r=>r.kind==='after'&&r.targetConcept?.elements.length===9));
 assert.ok(api.anatomicalConceptFunctionalRelations(c,members,a.parts).some(r=>r.kind==='after'));
});

test('Bladder muscle extensions retain their own controls, fixed pubic peers and an exhaustive peel',()=>{
 const a=read(),m=a.finiteUrinaryIntegration.muscleExtensions;if(!m)return;
 const r=api.urinaryWallRegistry(a);assert.ok(r);
 const pubic=m.pubicAttachments.map(e=>e.partId);
 assert.deepEqual(new Set(r.marker.termGroups.find(g=>g.term==='A08.3.01.020').partIds),new Set(pubic));
 assert.deepEqual(r.marker.termGroups.find(g=>g.term==='A08.3.01.021M').partIds,[m.posteriorPartId]);
 for(const e of m.pubicAttachments){
  const p=r.byId.get(e.partId),bone=r.byId.get(e.sourceContact.second);
  assert.equal(bone.system,'skeletal');assert.equal(e.sourceContact.triangles,1);
  assert.equal(api.displayLaterality(p.name).side,e.chart.endsWith('left')?'left':'right');
  assert.ok(p.tissueDrawGroups.some(g=>g.neighbor===bone.id));
 }
 const c=a.concepts.find(c=>c.id==='ZA:Urinary bladder'),steps=api.bowelWallPeelSteps(a,c);
 const muscle=steps.find(s=>s.name==='Outer muscular tissue');
 assert.ok([...pubic,m.posteriorPartId].every(id=>muscle.ids.includes(id)));
 for(const change of [a=>a.finiteUrinaryIntegration.muscleExtensions.pubicAttachments[0].sourceContact.actual_decoded_positions_and_opposed_words_match=false,
  a=>a.finiteUrinaryIntegration.muscleExtensions.pubicAttachments[0].sourceContact.second='ZA:Prostate',
  a=>a.parts.find(p=>p.id===a.finiteUrinaryIntegration.muscleExtensions.posteriorPartId).ta98Term='A08.3.01.020']){
  const broken=read();change(broken);assert.equal(api.urinaryWallRegistry(broken),null);
 }
});

test('Changed bank identities, incomplete seam proofs and broken ownership reject the new assembly',()=>{
 for(const edit of [a=>a.finiteUrinaryIntegration.charts[0].bankSha256='0'.repeat(64),
  a=>a.finiteUrinaryIntegration.externalInterfaces[0].triangles--,
  a=>a.parts.find(p=>p.id===a.finiteUrinaryIntegration.partIds[0]).normalSigns=undefined]){
  const a=read();edit(a);assert.equal(api.urinaryWallRegistry(a),null);
 }
});

test('A named crest or uvula cannot be replaced by a generic coat',()=>{
 const a=read();if(!hasFeatures(a))return;
 const g=a.finiteUrinaryIntegration.termGroups.find(g=>g.term==='A08.3.01.025');
 a.parts.find(p=>p.id===g.partIds[0]).tissueWallLayer='mucosa';
 assert.equal(api.urinaryWallRegistry(a),null);
});

test('Opening controls bind complete clinical rim curves to actual free mucosal faces',()=>{
 const a=read();if(!hasOpenings(a))return;const r=api.urinaryWallRegistry(a);assert.ok(r);
 const groups=r.marker.termGroups.filter(g=>['A08.3.01.026','A08.3.01.027M'].includes(g.term));
 assert.equal(groups.length,2);assert.deepEqual(groups.flatMap(g=>g.openingBoundaries.map(b=>b.name)).sort(),['left','neck','right']);
 for(const g of groups)for(const boundary of g.openingBoundaries){
  assert.ok(boundary.points.length>=12);assert.ok(boundary.maximumSourceOffsetM<=1e-7);
  for(const id of boundary.partIds){
   assert.ok(r.marker.termGroups.find(g=>g.term==='A08.3.01.023').partIds.includes(id));
   assert.equal(r.byId.get(id).ta98Term,'A08.3.01.023','The supporting wall remains mucosa; the witnessed curve is the opening');
   if(boundary.name!=='neck')assert.equal(api.displayLaterality(r.byId.get(id).name).side,boundary.name);
  }
 }
 for(const edit of [b=>b.points[0]=[0,0,0],b=>b.faceIds[0]=-1,b=>b.maximumSourceOffsetM=.001,b=>b.sourcePlanSha256='0'.repeat(64)]){
  const changed=read(),boundary=changed.finiteUrinaryIntegration.termGroups.find(g=>g.term==='A08.3.01.026').openingBoundaries[0];edit(boundary);
  assert.equal(api.urinaryWallRegistry(changed),null);
 }
});

test('Apical attachment keeps a reciprocal wall cap and a separately typed fascial reference contact',()=>{
 const a=read(),attachment=a.finiteUrinaryIntegration.apicalAttachment;if(!attachment)return;
 const r=api.urinaryWallRegistry(a);assert.ok(r);
 const p=r.byId.get(attachment.partId),wall=r.byId.get(attachment.wallPartId),ref=attachment.referenceContact;
 assert.equal(p.ta98Term,'A08.3.01.007');assert.equal(p.tissueWallLayer,'median-umbilical-ligament');
 assert.equal(p.tissueDrawGroups.find(g=>g.neighbor===wall.id).count,wall.tissueDrawGroups.find(g=>g.neighbor===p.id).count);
 assert.equal(ref.representation,'open-reference-sheet-contact');assert.equal(r.chartFor.has(ref.second),false);
 const concept=a.concepts.find(c=>c.id==='ZA:Urinary bladder');
 assert.equal(api.bowelWallPeelSteps(a,concept)[0].ids[0],p.id);
 for(const edit of [m=>m.referenceContact.actual_decoded_positions_and_same_words_match=false,
  m=>m.referenceContact.firstFace--,m=>m.sourcePlanSha256='invalid']){
  const changed=read();edit(changed.finiteUrinaryIntegration.apicalAttachment);assert.equal(api.urinaryWallRegistry(changed),null);
 }
});
