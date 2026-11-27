import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import{build}from'esbuild';import{createRequire}from'node:module';
const bundle=await build({stdin:{contents:`export * from './app/organ-models';export * from './app/spinal-cord-controls';export * from './app/hierarchy-navigation';export * from './app/guest-hierarchy';export * from './app/local-ta98-layer-catalogue';export * from './app/preferred-source-representations';`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const uiBundle=await build({stdin:{contents:`import React from'react';import{renderToStaticMarkup}from'react-dom/server';import Links from'./app/separate-organ-model-links';import GuestTree from'./app/guest-tree';export const render=props=>renderToStaticMarkup(React.createElement(Links,props));export const renderTree=props=>renderToStaticMarkup(React.createElement(GuestTree,props));`,resolveDir:process.cwd(),sourcefile:'organ-links-test.tsx',loader:'tsx'},bundle:true,platform:'node',format:'cjs',write:false,define:{'import.meta.env.DEV':'false'}});
const ui={exports:{}};new Function('require','module','exports',uiBundle.outputFiles[0].text)(createRequire(import.meta.url),ui,ui.exports);
const root='.local-models/organ-models/spinal-cord/',read=name=>JSON.parse(fs.readFileSync(root+name));
const detail=read('atlas.json'),guest=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));

test('the detailed cord is a separate valid model with its own complete control inventory',()=>{
 assert.equal(detail.scope,'organ');assert.equal(api.meshCount(detail.parts),3159);assert.ok(api.spinalCordRegistry(detail));
 assert.ok(detail.parts.every(p=>p.id.startsWith('TA98CORD:')));
 assert.equal(detail.chunks.length,1);assert.equal(detail.chunks[0].bytes,85_156_306);
 for(const mode of ['systems','regions','depth','guest:ta98']){
  const hierarchy=api.hierarchyNavigation(detail,mode,mode==='guest:ta98'?guest:undefined);
  assert.equal(api.meshCount(hierarchy.elements.map(id=>detail.parts.find(p=>p.id===id))),3159);
 }
});

for(const sex of ['female','male'])test(`${sex}: body container totals exclude the external model while its parenthesized link shows the separate mesh total`,()=>{
 const a=api.withHiddenReviewAlternatives(read(sex+'-body.json')),link=a.separateOrganModels[0],byId=new Map(a.parts.map(p=>[p.id,p]));
 assert.equal(link.meshCount,api.meshCount(detail.parts));assert.equal(link.overviewPartIds.length,sex==='female'?3:1);
 assert.ok(a.parts.every(p=>!p.id.startsWith('TA98CORD:')),'Detailed cells are absent from the body catalogue, rather than just hidden');
 assert.ok(detail.parts.every(p=>!byId.has(p.id)));
 const before=api.meshCount(a.parts);
 for(const mode of ['systems','regions','depth','guest:ta98']){
  const hierarchy=api.hierarchyNavigation(a,mode,mode==='guest:ta98'?guest:undefined);
  const visited=new Set(),cord=[];
  const visit=c=>{if(visited.has(c))return;visited.add(c);assert.ok(c.elements.every(id=>byId.has(id)));if(api.separateOrganModels(a,{id:c.id,name:c.name,term:c.ta98Term??c.terminology?.ta98}).length)cord.push(c);c.children?.forEach(visit);};visit(hierarchy);
  assert.ok(cord.length>=1,`${mode} exposes the organ's external model link`);
  for(const choice of cord){
   const html=ui.exports.render({atlas:a,id:choice.id,name:choice.name,term:choice.ta98Term??choice.terminology?.ta98,depth:2});
   assert.match(html,/\(Spinal cord · detailed/);assert.match(html,/3,159 meshes\)/);assert.match(html,/href="\?model=local-spinal-cord&amp;tree=systems"/);assert.doesNotMatch(html,/role="checkbox"/);
   assert.ok(api.meshCount(choice.elements.map(id=>byId.get(id)))<10);
  }
  assert.ok(api.meshCount(hierarchy.elements.map(id=>byId.get(id)))<=before);
 }
 assert.equal(api.meshCount(a.parts),before,'Rendering links cannot add the other model to any container');
 assert.equal(api.meshCount([...a.parts,...a.parts]),before,'Term aliases count each mesh once');
 const tree=ui.exports.renderTree({atlas:a,hierarchy:guest,state:{visible:['nervous'],hidden:[],guestQuery:'A14.1.02.001'},setState:()=>{},onChoose:()=>{}});
 assert.match(tree,/class="tree-row tree-leaf separate-organ-model"/);
 assert.match(tree,/3,159 meshes\)/);
 assert.equal(api.localTa98LayerCatalogue(sex==='male'?'local-ta98':'local-female-ta98',new URLSearchParams(),true),(sex==='male'?'/local-models/geometry-workbench/biodigital-male-circulation-v1/hand-current-palm-fit-v981/ordinary-atlas.json':'/local-models/organ-models/spinal-cord/female-body.json'));
});

test('both body catalogues preserve every retained source part and its geometry bank literally',()=>{
 for(const sex of ['female','male']){
  const original=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)),body=read(sex+'-body.json'),old=new Map(original.parts.map(p=>[p.id,p]));
  for(const p of body.parts){
   if(p.id.startsWith('SPINALOVERVIEW:'))continue;
   const q=old.get(p.id),{chunk,...now}=p,{chunk:oldChunk,...before}=q;
   assert.deepEqual(now,before);assert.deepEqual(body.chunks[chunk],original.chunks[oldChunk]);
  }
  if(sex==='male'){
   const inventory=JSON.parse(fs.readFileSync('reports/ta98-spinal-cord-inventory.json')),ids=new Set(body.parts.map(p=>p.id));
   for(const term of inventory.terms)if(term.authority.id.startsWith('A14.1.02.')&&term.authority.id!=='A14.1.02.001')for(const mesh of term.source_meshes.male)assert.ok(!ids.has(mesh.part_id));
   assert.ok(ids.has('LOCAL:male:spinal_cord:1863'));assert.ok(!body.defaultHidden.includes('LOCAL:male:spinal_cord:1863'));
   // The old transverse donor section has no embedded TA98 metadata. Its
   // horns and tracts must not escape into the lightweight cord a second time.
   const section=original.parts.filter(p=>p.id.startsWith('ZA:')&&p.system==='nervous'&&p.bounds[0][1]>1.442&&p.bounds[1][1]<1.449);
   assert.equal(section.length,28);
   for(const p of section)assert.ok(!ids.has(p.id),`${p.name} is a donor detail, not body overview tissue`);
   for(const id of ['ZA:Spinal ganglion.l','ZA:Spinal ganglion.r'])if(old.has(id))assert.ok(ids.has(id),'Spinal ganglia remain in the body');
  }
 }
});

test('the whole-body grey/white/pial charts preserve finite reciprocal interfaces',()=>{
 const a=read('female-body.json'),parts=a.parts.filter(p=>p.id.startsWith('SPINALOVERVIEW:'));assert.equal(parts.length,3);
 const bank=fs.readFileSync(root+'overview.bin'),byId=new Map(parts.map(p=>[p.id,p]));
 for(const p of parts){
  assert.ok(p.pairedTissue&&p.sectionSolid);assert.ok(p.indexCount>0&&p.vertexCount>0);
  let end=0;for(const group of p.tissueDrawGroups){assert.equal(group.start,end);end+=group.count;if(group.neighbor)assert.equal(byId.get(group.neighbor).tissueDrawGroups.filter(g=>g.neighbor===p.id).reduce((n,g)=>n+g.count,0),group.count);}
  assert.equal(end,p.indexCount);assert.ok(p.indices+p.indexCount*4<=bank.byteLength);
 }
 const original=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
 const evidence=read('evidence.json');assert.ok(evidence.sourceCoordinatesAndNormalWordsPreserved);assert.equal(evidence.models.female.sourceMeshes-evidence.models.female.bodyMeshes,3156);
 const old=new Map(original.parts.map(p=>[p.id,p]));for(const p of a.parts){if(p.id.startsWith('SPINALOVERVIEW:'))continue;const q=old.get(p.id),{chunk,...now}=p,{chunk:oldChunk,...before}=q;assert.deepEqual(now,before);assert.deepEqual(a.chunks[chunk],original.chunks[oldChunk]);}
});
