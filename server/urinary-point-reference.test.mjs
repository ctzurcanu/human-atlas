import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/urinary-wall';export * from './app/anatomy-label-plan';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const read=()=>JSON.parse(fs.readFileSync(process.env.ATLAS_URINARY_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
const term='A08.3.01.002',conceptId='URINARY:male:bladder:'+term;
function triangle(atlas,part,face){
 const url=atlas.chunks[part.chunk].url,path=url.startsWith('/local-models/')?'.local-models/'+url.slice('/local-models/'.length):'public'+url;
 const buffer=fs.readFileSync(path);
 return Array.from({length:3},(_,corner)=>{
  const index=buffer.readUInt32LE(part.indices+(face*3+corner)*4);
  return Array.from({length:3},(_,k)=>buffer.readFloatLE(part.positions+(index*3+k)*4)+(part.geometryOrigin?.[k]??0));
 });
}
function fixture(){
 const atlas=read(),marker=atlas.finiteUrinaryIntegration;
 marker.termGroups=marker.termGroups.filter(g=>g.term!==term);atlas.concepts=atlas.concepts.filter(c=>c.id!==conceptId);
 const part=atlas.parts.find(p=>p.id===marker.apicalAttachment.wallPartId),group=part.tissueDrawGroups.find(g=>g.neighbor===marker.apicalAttachment.partId);
 const face=group.start/3,t=triangle(atlas,part,face),point=[0,1,2].map(k=>(t[0][k]+t[1][k]+t[2][k])/3);
 const choice={id:conceptId,name:'Apex of bladder',elements:[part.id],ta98Term:term,ta98Kind:'exact'};
 atlas.concepts.push(choice);marker.termGroups.push({term,conceptId,partIds:[part.id],pointReference:{representation:'point-landmark',partId:part.id,point,face,barycentric:[1/3,1/3,1/3],sourceBankSha256:marker.charts.find(c=>c.partIds.includes(part.id)).bankSha256,sourcePlanSha256:marker.apicalAttachment.sourcePlanSha256}});
 return {atlas,choice,part,point};
}
test('Apex uses its ordinary named label at its real attachment, with no new solid',()=>{
 const before=read(),{atlas,choice,part,point}=fixture();
 assert.equal(atlas.parts.length,before.parts.length);assert.deepEqual(atlas.parts,before.parts);
 assert.deepEqual(api.urinaryPointReference(atlas,choice),{partId:part.id,point});
 const labels=api.selectedAnatomyLabels(atlas,choice.elements,choice,false);
 assert.equal(labels.length,1);assert.equal(labels[0].name,'Apex of bladder');assert.deepEqual(labels[0].anchor,{partId:part.id,point});
 assert.deepEqual(labels[0].components,[]);assert.deepEqual(api.anatomyLabelChoice(labels[0]).elements,[part.id]);
 assert.ok(!api.layerLabelTargets(labels,atlas.parts).has(choice.id));
});
test('Host tissue alone does not acquire an apex label',()=>{
 const {atlas,part}=fixture(),choice=atlas.concepts.find(c=>c.id===part.conceptId);
 assert.equal(api.urinaryPointReference(atlas,choice),undefined);
 assert.ok(api.selectedAnatomyLabels(atlas,[part.id],choice,false).every(item=>!item.anchor));
});
test('A stale bank or a face outside the actual ligament attachment rejects the point',()=>{
 for(const mutate of [ref=>{ref.sourceBankSha256='f'.repeat(64);},(ref,part)=>{ref.face=part.tissueDrawGroups.find(g=>g.neighbor!==undefined&&g.neighbor!==read().finiteUrinaryIntegration.apicalAttachment.partId).start/3;}]){
  const {atlas,choice,part}=fixture(),ref=atlas.finiteUrinaryIntegration.termGroups.find(g=>g.term===term).pointReference;
  mutate(ref,part);assert.equal(api.urinaryWallRegistry(atlas),null);assert.equal(api.urinaryPointReference(atlas,choice),undefined);
 }
});
test('Proposed or installed apex replays on its actual encoded triangle',t=>{
 const atlas=read(),group=atlas.finiteUrinaryIntegration.termGroups.find(g=>g.term===term);
 if(!group?.pointReference){t.skip('Legacy revision has no apex reference');return;}
 const ref=group.pointReference,part=atlas.parts.find(p=>p.id===ref.partId),source=triangle(atlas,part,ref.face);
 const actual=[0,1,2].map(k=>source.reduce((sum,v,i)=>sum+v[k]*ref.barycentric[i],0));
 assert.ok(Math.hypot(...actual.map((v,k)=>v-ref.point[k]))<1e-12);
 assert.ok(api.urinaryWallRegistry(atlas));
});

const readFemale=()=>JSON.parse(fs.readFileSync(process.env.ATLAS_URINARY_FEMALE_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
test('female apex transfers the source point onto its actual exterior facet without adding geometry',t=>{
 const atlas=readFemale(),marker=atlas.finiteUrinaryIntegration,g=marker.termGroups.find(g=>g.term===term),transfer=marker.apexTransfer;
 if(!transfer&&!process.env.ATLAS_URINARY_FEMALE_CANDIDATE){t.skip('Female transfer is staged and not installed');return;}
 assert.ok(g?.pointReference&&transfer);assert.ok(api.urinaryWallRegistry(atlas));
 const ref=g.pointReference,part=atlas.parts.find(p=>p.id===ref.partId),t=triangle(atlas,part,ref.face);
 const point=[0,1,2].map(k=>t.reduce((sum,v,i)=>sum+v[k]*ref.barycentric[i],0));
 assert.ok(Math.hypot(...point.map((v,k)=>v-ref.point[k]))<1e-12);
 const mapped=transfer.matrix.slice(0,3).map(row=>row[3]+transfer.sourcePoint.reduce((sum,v,i)=>sum+v*row[i],0));
 assert.ok(Math.hypot(...mapped.map((v,i)=>v-ref.point[i]))<=transfer.maximumSourceOffsetM+1e-12);
 const choice=atlas.concepts.find(c=>c.id===g.conceptId),labels=api.selectedAnatomyLabels(atlas,choice.elements,choice,false);
 assert.equal(labels.length,1);assert.equal(labels[0].name,'Apex of bladder');assert.deepEqual(labels[0].anchor,{partId:part.id,point:ref.point});assert.deepEqual(labels[0].components,[]);
 assert.ok(!atlas.parts.some(p=>p.ta98Term===term),'The point must not become a separate solid');
 assert.equal(transfer.sourcePartId,'BLADDER:male:body:adventitia');
});
test('female apex rejects stale banks, moved points, reflection, shear and non-exterior facets',t=>{
 if(!readFemale().finiteUrinaryIntegration.apexTransfer&&!process.env.ATLAS_URINARY_FEMALE_CANDIDATE){t.skip('Female transfer is staged and not installed');return;}
 for(const mutate of [
  (m,g)=>g.pointReference.sourceBankSha256='e'.repeat(64),
  (m,g)=>g.pointReference.point[0]+=.0001,
  m=>m.apexTransfer.matrix[0][0]*=-1,
  m=>m.apexTransfer.matrix[0][1]+=.01,
  (m,g,a)=>{const p=a.parts.find(p=>p.id===g.pointReference.partId);g.pointReference.face=p.tissueDrawGroups.find(q=>q.neighbor).start/3;},
  m=>m.apexTransfer.maximumSourceOffsetM=.001,
 ]){
  const a=readFemale(),m=a.finiteUrinaryIntegration,g=m.termGroups.find(g=>g.term===term);mutate(m,g,a);assert.equal(api.urinaryWallRegistry(a),null);
 }
});
