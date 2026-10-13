import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/bowel-wall-peeling';export * from './app/viewer-state';export * from './app/requested-selection-choice';export * from './app/guest-hierarchy';export * from './app/anatomical-relations';export * from './app/gallbladder-wall';export * from './app/hepatic-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const read=p=>JSON.parse(fs.readFileSync(p));
test('gallbladder coat route remains local, opt-in and female only',()=>{
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('hepatic=coats'),true),'/local-models/ta98-review/gallbladder-coats-preview/female/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('hepatic=coats'),false),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('hepatic=coats'),true),'/local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json');
 for(const q of ['hepatic=coats&bowel=source','hepatic=coats&delivery=1'])assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams(q),true),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('hepatic=coats&bowel=prototype'),true),'/local-models/ta98-review/bowel-wall-main-prototype/female/atlas.json');
});
test('ordinary female TA98 delivery exposes installed coats without a preview flag',{skip:!fs.existsSync('reports/ta98-female/gallbladder-current-coats-install.json')},()=>{
 const install=read('reports/ta98-female/gallbladder-current-coats-install.json'),atlas=read(install.targetCatalogue);
 assert.ok(read('reports/ta98-female/gallbladder-current-coats-installed-validation.json').ordinaryDeliveryGatePass);
 const model=api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams(),true);assert.equal(model,'/local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json');
 const ids=api.resolveSelection(atlas,[atlas.gallbladderWallPreview.sourcePartId]);assert.deepEqual(ids,atlas.gallbladderWallPreview.partIds);
 const choice=api.requestedSelectionChoice(atlas,[atlas.gallbladderWallPreview.sourcePartId],ids);assert.equal(api.bowelWallPeelSteps(atlas,choice).length,4);
 const frame=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5],state=api.readViewUrl(new URLSearchParams({select:choice.id,frame:frame.join(',')}).toString(),atlas,{selected:[],visible:[],view:'three-quarter'});
 assert.deepEqual(state.camera,frame);assert.deepEqual(state.selected,ids);
 assert.equal(new URL(api.viewUrl('http://localhost:3016/?hepatic=coats','local-female-ta98',state,frame,atlas)).searchParams.has('hepatic'),false,'ordinary saved views must stay on ordinary delivery');
 assert.equal(atlas.gallbladderWallIntegration.anatomicallyAccepted,false);
 const artery=atlas.parts.find(p=>p.id==='LOCAL:female:r_hepatic_artery:1719');
 assert.match(api.hepaticDescription(artery).text,/right hepatic artery supplies the liver/);
 assert.equal(api.hepaticDescription({...artery,id:'HRA:foreign'}),undefined);
 assert.equal(api.hepaticDescription({...artery,suppressed:true}),undefined);
 const members=ids.map(id=>atlas.parts.find(p=>p.id===id)),links=api.anatomicalConceptFunctionalRelations(choice,members,atlas.parts);
 for(const [kind,back] of [['arterial','supplies'],['venous','drains'],['innervation','innervates']])for(const row of links.filter(r=>r.kind===kind)){
  const reciprocal=api.anatomicalRelations(row.target,atlas.parts).filter(r=>r.kind===back&&r.targetConcept?.id===choice.id);
  assert.equal(reciprocal.length,1);assert.deepEqual(reciprocal[0].targetConcept.elements,ids);assert.ok(!reciprocal[0].target.suppressed);
  assert.deepEqual(api.resolveSelection(atlas,[reciprocal[0].targetConcept.id]),ids);
 }
});
const path='.local-models/ta98-review/gallbladder-coats-preview/female/atlas.json';
test('actual four coats expose three canonical TA98 terms, whole-organ relations and every peel stage',{skip:!fs.existsSync(path)},()=>{
 const proof=read('reports/ta98-female/gallbladder-current-coats-validation.json');assert.ok(proof.finiteCoatDeliveryGatePass);
 const atlas=read(path),marker=atlas.gallbladderWallPreview,ids=marker.partIds,parts=ids.map(id=>atlas.parts.find(p=>p.id===id));
 const selected=api.resolveSelection(atlas,[marker.sourcePartId]);assert.deepEqual(new Set(selected),new Set(ids));
 const choice=api.requestedSelectionChoice(atlas,[marker.sourcePartId],selected);assert.equal(choice.name,'Gallbladder');
 const steps=api.bowelWallPeelSteps(atlas,choice);assert.equal(steps.length,4);assert.match(steps[0].name,/allocation unresolved/);
 for(let i=0;i<=4;i++)assert.equal(api.wallPeelDepth(steps,steps.slice(0,i).flatMap(s=>s.ids)),i);
 assert.equal(api.bowelWallPeelSteps(atlas,{...choice,elements:selected.slice(1)}),null);
 const nodes=new Map(),walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};api.resolveGuestHierarchy(atlas,api.parseGuestHierarchy(read('public/assets/ta98.json'))).forEach(walk);
 for(const [code,key] of [['A05.8.02.007','subserosa'],['A05.8.02.008','muscular'],['A05.8.02.009','mucosa']])assert.deepEqual(nodes.get(code).directParts.map(p=>p.id),['LOCAL:female:GALLWALL:'+key]);
 assert.ok(!nodes.get('A05.8.02.006').parts.some(p=>ids.includes(p.id)),'unassigned covering cannot become serosa');
 const graph=api.gallbladderWallRelations(atlas.parts),duct=atlas.parts.find(p=>p.id==='LOCAL:female:bile_duct:1031');
 assert.ok(!atlas.materials[duct.material].vertexColors&&duct.colors===undefined,'a material cannot require the removed color attribute');
 for(const part of parts){assert.ok(graph.get(part.id).some(r=>r.kind==='continuous'&&r.target.id===duct.id));assert.ok(graph.get(duct.id).some(r=>r.kind==='continuous'&&r.target.id===part.id));assert.equal(part.positions,parts[0].positions);assert.equal(part.normals,parts[0].normals);assert.equal(part.chunk,parts[0].chunk);assert.equal(part.uvs,undefined);assert.ok(!atlas.materials[part.material].map);}
 for(const part of parts){assert.ok(api.anatomicalRelations(part,atlas.parts).some(r=>r.kind==='continuous'&&r.target.id===duct.id));assert.ok(api.anatomicalRelations(duct,atlas.parts).some(r=>r.kind==='continuous'&&r.target.id===part.id));}
 const relations=api.anatomicalConceptFunctionalRelations(choice,parts,atlas.parts);
 assert.equal(relations.filter(r=>r.kind==='adjacent'&&r.target.id===duct.id).length,1);
 assert.ok(relations.some(r=>r.kind==='arterial'));assert.ok(relations.some(r=>r.kind==='venous'));assert.ok(relations.some(r=>r.kind==='innervation'));
 assert.ok(!api.anatomicalConceptFunctionalRelations({...choice,elements:ids.slice(1)},parts.slice(1),atlas.parts).some(r=>['arterial','venous','innervation'].includes(r.kind)));
 const frame=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5],q=new URLSearchParams({model:'local-female-ta98',hepatic:'coats',frame:frame.join(','),select:marker.sourcePartId});
 const state=api.readViewUrl(q.toString(),atlas,{selected:[],visible:[],view:'three-quarter'});assert.deepEqual(state.camera,frame);assert.deepEqual(new Set(state.selected),new Set(ids));
 assert.equal(new URL(api.viewUrl('http://localhost:3016/','local-female-ta98',state,frame,atlas)).searchParams.get('hepatic'),'coats');
 for(const broken of [{...parts[1],sourceOffset:[0,.01,0]},{...parts[1],bowelWallSource:'HRA:foreign'}]){
  const snapshot=atlas.parts.map(p=>p.id===broken.id?broken:p),selection=parts.map(p=>p.id===broken.id?broken:p);
  assert.equal(api.gallbladderWholeSource(selection,snapshot),undefined);
  assert.equal(api.bowelWallPeelSteps({...atlas,parts:snapshot},choice),null);
 }
});
