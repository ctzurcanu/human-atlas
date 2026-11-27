import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
const built=await build({entryPoints:['app/carotid-surface-features.ts'],bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const hostId='TA98PART:ZA:Internal carotid artery.l:A12.2.06.002';
const featureId='TA98REFERENCE:male:carotid-sinus:left:medial-proximal-ICA';
const term='A12.2.04.008',bankSha='e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36';
const bankUrl='/local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1/packed-left-v2/cells.bin';
const marker=a=>a.carotidSurfaceFeatureIntegration;
function fixture(){
 const host={id:hostId,conceptId:hostId,name:'Cervical part of internal carotid artery (left)',ta98Term:'A12.2.06.002',sourceId:'Internal carotid artery.l',system:'arterial',pairedTissue:true,sectionSolid:true,sectionAssembly:'left-carotid-lingual-coupled-finite-wall:male',carotidFiniteWallRevision:'left-four-cell-common-bank-v2',carotidFiniteWallBankSha256:bankSha,geometryOrigin:[0,1.625,0],chunk:0,positions:0,normals:195540,normalSigns:806208,uvs:293312,indices:480120,vertexCount:16295,indexCount:7848,bounds:[[.024,1.612,-.008],[.039,1.684,.005]],provenance:{url:'https://github.com/Z-Anatomy/Models-of-human-anatomy'},tissueDrawGroups:[{start:0,count:3810},{start:3810,count:3882},{start:7692,count:156,neighbor:'TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall'}]};
 const group={term,conceptId:featureId,partIds:[hostId],representation:'surface-region',referenceOnly:true,boundaryProvisional:true,surface:{partId:hostId,faces:[1963,1964],sourceBankSha256:bankSha},anchor:{partId:hostId,face:1964,barycentric:[.2,.3,.5],point:[.033,1.614,-.004]}};
 return {sex:'male',parts:[host],concepts:[{id:hostId,name:host.name,elements:[hostId],ta98Term:'A12.2.06.002',ta98Kind:'exact'},{id:featureId,name:api.carotidSurfaceReferenceName,description:api.carotidSurfaceReferenceDescription,elements:[hostId],ta98Term:term,ta98Kind:'exact'}],chunks:[{url:bankUrl,bytes:822503}],carotidSurfaceFeatureIntegration:{schema:'carotid-surface-features/v1',side:'left',hostPartId:hostId,bankSha256:bankSha,bankUrl,sourcePlanSha256:'a'.repeat(64),provisional:true,termGroups:[group]}};
}
test('a qualified surface reference reuses its exact carrier without adding anatomy',()=>{
 const a=fixture(),before=structuredClone(a);assert.ok(api.carotidFeatureRegistry(a));assert.equal(a.parts.length,1);
 assert.equal(api.carotidSurfaceTermParts(a).get(term)[0],a.parts[0]);assert.deepEqual(a,before);
 const choice=a.concepts[1];assert.equal(api.selectedCarotidSurfaceFeature(a,choice),marker(a).termGroups[0]);
 assert.deepEqual(api.carotidSurfaceLabelReference(a,choice),{partId:hostId,point:marker(a).termGroups[0].anchor.point});
 assert.equal(api.selectedCarotidSurfaceFeature(a,a.concepts[0]),undefined);
 assert.equal(api.selectedCarotidSurfaceFeature(a,{...choice,elements:[hostId,hostId]}),undefined);
});
test('both authority occurrences select one qualified left surface and retain their term identity',()=>{
 const a=fixture(),before=structuredClone(a),group=marker(a).termGroups[0];
 const parts=api.carotidSurfaceTermParts(a),references=api.carotidReferenceChoices(a);
 assert.deepEqual([...parts.keys()],['A12.2.04.008','A12.2.06.003']);
 for(const occurrence of parts.keys()){
  assert.deepEqual(parts.get(occurrence),[a.parts[0]]);
  const reference=references.get(occurrence);assert.ok(reference);
  assert.equal(reference.choice.id,featureId);assert.equal(reference.choice.ta98Term,occurrence);
  assert.equal(reference.choice.terminology.ta98,occurrence);
  assert.match(reference.choice.description,/Qualified reference/);
  assert.match(reference.choice.description,/remain unverified/);
  assert.equal(api.carotidSurfaceChoiceDescription(reference.choice),api.carotidSurfaceReferenceDescription);
  assert.equal(api.carotidSurfaceChoiceDescription({...reference.choice,ta98Term:'A12.2.06.002'}),undefined);
  assert.equal(api.carotidSurfaceChoiceDescription(a.concepts[0]),undefined);
  const node={id:`hierarchy:guest:ta98:${occurrence}`,elements:[hostId],ta98Term:occurrence,ta98Kind:'exact'};
  assert.equal(api.selectedCarotidSurfaceFeature(a,node),group);
  assert.deepEqual(api.carotidSurfaceLabelReference(a,node),{partId:hostId,point:group.anchor.point});
  assert.equal(api.carotidReferenceActive(reference,node),true);
  assert.equal(api.toggleCarotidReferenceChoice(reference,node).id,hostId);
 }
 assert.equal(api.selectedCarotidSurfaceFeature(a,{...a.concepts[1],ta98Term:'A12.2.06.004'}),undefined);
 assert.equal(references.has('A12.2.06.004'),false);
 assert.deepEqual(a,before);assert.equal(marker(a).termGroups.length,1);
});
test('strict bank, side, source frame and carrier buffers reject stale metadata',()=>{
 const changes=[a=>a.sex='female',a=>marker(a).side='right',a=>marker(a).schema='carotid-surface-features/v2',a=>marker(a).bankSha256='0'.repeat(64),a=>marker(a).bankUrl='/other.bin',a=>marker(a).sourcePlanSha256='bad',a=>marker(a).provisional=false,a=>a.chunks[0].url='/other.bin',a=>a.chunks[0].bytes++,a=>a.parts[0].id=hostId.replace('.l:', '.r:'),a=>a.parts[0].sourceId='Internal carotid artery.r',a=>a.parts[0].carotidFiniteWallRevision='left-four-cell-common-bank-v3',a=>a.parts[0].carotidFiniteWallBankSha256='0'.repeat(64),a=>a.parts[0].geometryOrigin=[0,0,0],a=>a.parts[0].sourceOffset=[0,0,0],a=>a.parts[0].geometryScale=2,a=>a.parts[0].geometryInstance={scale:1,rotation:[0,0,0,1]},a=>a.parts[0].suppressed=true,a=>a.defaultHidden=[hostId],a=>delete a.parts[0].normalSigns,a=>delete a.parts[0].uvs,a=>a.parts[0].uvs++,a=>a.parts[0].indices++,a=>a.parts.push({...a.parts[0]}),a=>a.concepts[1].ta98Kind='parent',a=>a.concepts[1].ta98Term='A12.2.04.007'];
 for(const change of changes){const a=fixture();change(a);assert.equal(api.carotidFeatureRegistry(a),null,change.toString());assert.equal(api.carotidReferenceChoices(a).size,0);}
});
test('only actual outer facets and valid anchored barycentrics establish selection ink',()=>{
 const changes=[a=>marker(a).termGroups[0].surface.faces=[1269],a=>marker(a).termGroups[0].surface.faces=[2564],a=>marker(a).termGroups[0].surface.faces=[999999],a=>marker(a).termGroups[0].surface.faces=[1964,1964],a=>marker(a).termGroups[0].surface.faces=[1964.5],a=>marker(a).termGroups[0].surface.faces=[],a=>a.parts[0].tissueDrawGroups[1].count-=3,a=>a.parts[0].tissueDrawGroups[0]=null,a=>a.parts[0].tissueDrawGroups[1].neighbor='unrelated',a=>marker(a).termGroups[0].anchor.face=2000,a=>marker(a).termGroups[0].anchor.barycentric=[.2,.3,.6],a=>marker(a).termGroups[0].anchor.barycentric=[-.1,.6,.5],a=>marker(a).termGroups[0].anchor.point=[NaN,1,2],a=>marker(a).termGroups[0].anchor.point=[-.033,1.614,-.004],a=>marker(a).termGroups[0].referenceOnly=false,a=>marker(a).termGroups[0].representation='finite-tissue',a=>marker(a).termGroups=[null]];
 for(const change of changes){const a=fixture();change(a);assert.equal(api.carotidFeatureRegistry(a),null,change.toString());}
});
test('reference choices and Systems aliases keep selection on the existing ICA',()=>{
 const a=fixture(),ref=api.carotidReferenceChoices(a).get(term);assert.ok(ref);assert.match(ref.choice.name,/medial wall reference/);assert.match(ref.choice.description,/Qualified reference/);
 const shown=api.toggleCarotidReferenceChoice(ref,ref.host);assert.equal(shown.id,featureId);
 const cleared=api.toggleCarotidReferenceChoice(ref,shown);assert.equal(cleared.id,hostId);assert.deepEqual(cleared.elements,shown.elements);assert.equal(api.carotidReferenceActive(ref,cleared),false);
 const entries=[{id:hostId,name:a.parts[0].name,parts:[a.parts[0]],system:'arterial',region:'Head & neck',location:['Neck'],ancestry:['Common carotid artery'],terminology:{}}];
 const aliases=api.carotidSurfaceReferenceEntries(a,entries);assert.equal(aliases.length,1);assert.equal(aliases[0].id,featureId);assert.equal(aliases[0].parts[0],a.parts[0]);assert.equal(aliases[0].ta98Term,term);assert.equal(aliases[0].taParentConcept,'ZA:Left common carotid artery');assert.equal(entries.length,1);
 assert.equal(api.carotidSurfaceReferenceEntries(a,[...entries,...aliases]).length,0);
});
const dataDir='.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1/qualified-carotid-sinus-surface-reference-v1';
const dataAvailable=fs.existsSync(dataDir+'/feature.json')&&fs.existsSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json');
test('local actual-bank data replays the outer-face and barycentric anchor binding',{skip:!dataAvailable&&'Local source data is not redistributed.'},()=>{
 const feature=JSON.parse(fs.readFileSync(dataDir+'/feature.json'));const a=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
 a.concepts=a.concepts.filter(c=>c.id!==featureId);a.concepts.push({id:featureId,name:api.carotidSurfaceReferenceName,description:api.carotidSurfaceReferenceDescription,elements:[hostId],ta98Term:term,ta98Kind:'exact'});
 a.carotidSurfaceFeatureIntegration={schema:'carotid-surface-features/v1',side:'left',hostPartId:hostId,bankSha256:bankSha,bankUrl,sourcePlanSha256:createHash('sha256').update(fs.readFileSync(dataDir+'/independent-bounds.json')).digest('hex'),provisional:true,termGroups:[feature.applicationTermGroup]};
 const registry=api.carotidFeatureRegistry(a);assert.ok(registry);const group=registry.termGroups[0],host=a.parts.find(p=>p.id===hostId);const bytes=fs.readFileSync('.local-models/'+bankUrl.split('/local-models/')[1]);assert.equal(createHash('sha256').update(bytes).digest('hex'),bankSha);assert.equal(group.surface.faces.length,58);assert.ok(group.surface.faces.every(f=>f>=1270&&f<=2563));assert.equal(group.anchor.face,1964);
 const points=[0,1,2].map(corner=>{const index=bytes.readUInt32LE(host.indices+group.anchor.face*12+corner*4);return [0,1,2].map(axis=>bytes.readFloatLE(host.positions+index*12+axis*4)+host.geometryOrigin[axis]);});
 const replay=[0,1,2].map(axis=>points.reduce((sum,p,i)=>sum+p[axis]*group.anchor.barycentric[i],0));assert.ok(Math.hypot(...replay.map((x,i)=>x-group.anchor.point[i]))<1e-12);
 const bounds=JSON.parse(fs.readFileSync(dataDir+'/independent-bounds.json'));assert.equal(bounds.passed,true);assert.equal(bounds.sourceBankSha256,bankSha);assert.deepEqual(bounds.actualHostOuterChartFaceRange,[1270,2563]);assert.equal(bounds.anchorActualReconstructionErrorMetres,0);
});
