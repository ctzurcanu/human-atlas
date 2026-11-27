import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import{build}from'esbuild';
const compiled=await build({stdin:{contents:`export * from './app/tongue-controls';export * from './app/tissue-wall-material';export * from './app/anatomical-reference-selection';export * from './app/hierarchy-navigation';export * from './app/anatomy-hierarchy';export * from './app/anatomy-label-plan';export * from './app/guest-hierarchy';export * from './app/organ-models';export * from './app/surface-feature-overlay';export {Group,Mesh,BufferGeometry,Float32BufferAttribute,Matrix4} from'three';`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const hash='a'.repeat(64),host='NATIVE:Tongue',muscle='TONGUE:male:genioglossus',whole='TONGUE:male:whole';
const fixture=()=>{
 const part=(id,system)=>({id,conceptId:id,name:id===host?'Tongue':'Genioglossus',system,chunk:0,positions:0,normals:48,indices:72,vertexCount:4,indexCount:12,geometryOrigin:[0,1.6,0],bounds:[[0,1.6,0],[.02,1.62,.02]],material:'m'});
 const anchor=(face=0)=>({partId:host,face,barycentric:[1,0,0],point:[0,1.6,0]});
 const groups=[
  {key:'tongue',name:'Tongue',term:'A05.1.04.001',conceptId:whole,partIds:[host,muscle],kind:'exact'},
  {key:'muscles',parent:'tongue',name:'Muscles of tongue',term:'A04.1.05.001',conceptId:'TONGUE:male:muscles',partIds:[muscle],kind:'exact'},
  {key:'genioglossus',parent:'muscles',name:'Genioglossus',term:'A05.1.04.101',conceptId:muscle,partIds:[muscle],kind:'exact'},
  {key:'dorsum',parent:'tongue',name:'Dorsum of tongue',term:'A05.1.04.004',conceptId:'TONGUE:male:feature:dorsum',partIds:[host],kind:'exact',representation:'surface-region',surface:{partId:host,faces:[0],sourceBankSha256:hash},anchor:anchor()},
  {key:'apex',parent:'tongue',name:'Apex of tongue',term:'A05.1.04.010',conceptId:'TONGUE:male:feature:apex',partIds:[host],kind:'exact',representation:'point-landmark',anchor:anchor()},
  {key:'margin',parent:'tongue',name:'Margin of tongue',term:'A05.1.04.009',conceptId:'TONGUE:male:feature:margin',partIds:[host],kind:'exact',representation:'boundary-curve',anchor:anchor(),curves:[{partId:host,points:[[0,1.6,0],[.02,1.6,0],[0,1.62,0]],closed:false,sourceBankSha256:hash}]},
 ];
 return {version:'fixture; binding semantics only',sex:'male',parts:[part(host,'digestive'),part(muscle,'muscular')],concepts:[...groups.map(g=>({id:g.conceptId,name:g.name,elements:g.partIds,ta98Term:g.term,ta98Kind:g.kind})),{id:host,name:'Tongue',elements:[host],ta98Term:'A05.1.04.001',ta98Kind:'exact'}],chunks:[{url:'/local-models/test-tongue.bin',bytes:120}],materials:{m:{color:[.8,.3,.3]}},triangles:8,finiteTongueIntegration:{schema:'tongue-local/v1',wholeConceptId:whole,partIds:[host,muscle],termGroups:groups,carrierBankHashes:{[host]:hash}}};
};

const coupledFixture=()=>{
 const old=fixture(),base=old.parts[0],mucosa='TA98TONGUE:male:mucosa',apo='TA98TONGUE:male:aponeurosis';
 const ids=[mucosa,apo,'GG:l','GG:r','HG:l','HG:r','SL','IL:l','IL:r','TV','SG:l','SG:r'];
 const parts=ids.map(id=>({...base,id,conceptId:id,name:id,system:[mucosa,apo].includes(id)?'digestive':'muscular'}));
 const mk=(key,name,term,partIds,parent='muscles')=>({key,name,term,partIds,parent,conceptId:'TA98TONGUE:male:control:'+key,kind:'exact'});
 const groups=[{...mk('tongue','Tongue','A05.1.04.001',ids),parent:undefined},
  mk('muscles','Muscles of tongue','A04.1.05.001',ids.slice(2),'tongue'),
  mk('genioglossus','Genioglossus','A05.1.04.101',['GG:l','GG:r']),mk('hyoglossus','Hyoglossus','A05.1.04.102',['HG:l','HG:r']),
  mk('styloglossus','Styloglossus','A05.1.04.105',['SG:l','SG:r']),mk('superior-longitudinal','Superior longitudinal muscle','A05.1.04.106',['SL']),
  mk('inferior-longitudinal','Inferior longitudinal muscle','A05.1.04.107',['IL:l','IL:r']),
  mk('transverse','Transverse muscle of tongue','A05.1.04.108',['TV']),mk('vertical','Vertical muscle of tongue','A05.1.04.109',['TV']),
  mk('mucosa','Mucous membrane of tongue','A05.1.04.011',[mucosa],'tongue'),mk('aponeurosis','Lingual aponeurosis','A05.1.04.025',[apo],'tongue')];
 const oldFeatures=old.finiteTongueIntegration.termGroups.filter(g=>g.representation).map(g=>{
  const n=structuredClone(g);n.partIds=[mucosa];n.anchor.partId=mucosa;
  if(n.surface){n.surface.partId=mucosa;n.surface.sourceBankSha256=hash;}
  n.curves?.forEach(c=>{c.partId=mucosa;c.sourceBankSha256=hash;});return n;
 });
 const dorsal=oldFeatures.find(g=>g.key==='dorsum');
 const extra=[['anterior-dorsum','Anterior part of dorsum','A05.1.04.005'],['posterior-dorsum','Posterior part of dorsum','A05.1.04.006'],['inferior-surface','Inferior surface of tongue','A05.1.04.007']]
  .map(([key,name,term])=>({...structuredClone(dorsal),key,name,term,conceptId:'TA98TONGUE:male:feature:'+key}));
 groups.push(...oldFeatures,...extra);
 parts.push({...base,id:'ZA:Tongue',conceptId:'ZA:Tongue',suppressed:true});
 const concepts=groups.map(g=>({id:g.conceptId,name:g.name,elements:g.partIds,ta98Term:g.term,ta98Kind:'exact'}));
 concepts.push({id:'ZA:Tongue',name:'Tongue',elements:[...ids]});
 return {...old,parts,concepts,finiteTongueIntegration:{schema:'tongue-local/v1',wholeConceptId:groups[0].conceptId,
  partIds:ids,termGroups:groups,carrierBankHashes:{[mucosa]:hash},legacyWholeAlias:true,sourcePartId:'ZA:Tongue'}};
};

test('coupled seventeen-control marker counts twelve active sources and retains the legacy whole alias',()=>{
 const atlas=coupledFixture(),before=JSON.stringify(atlas),registry=api.tongueRegistry(atlas);
 assert.ok(registry);assert.equal(registry.marker.termGroups.length,17);assert.equal(api.meshCount(registry.parts),12);
 assert.ok(!registry.parts.some(p=>p.id==='ZA:Tongue'));
 const legacy=atlas.concepts.find(c=>c.id==='ZA:Tongue');
 assert.equal(api.selectedTongueFeature(atlas,legacy)?.key,'tongue');
 assert.equal(api.tongueComponentChoices(atlas,legacy).length,16);
 const labels=api.selectedAnatomyLabels(atlas,legacy.elements,legacy,true);
 assert.equal(labels.length,16);assert.equal(labels.filter(l=>l.anchor).length,6);
 for(const label of labels){const choice=api.anatomyLabelChoice(label);assert.deepEqual(choice.elements,atlas.concepts.find(c=>c.id===label.id).elements);}
 const tv=registry.marker.termGroups.filter(g=>['A05.1.04.108','A05.1.04.109'].includes(g.term));
 assert.deepEqual(tv[0].partIds,tv[1].partIds);assert.equal(tv[0].partIds.length,1);
 assert.equal(JSON.stringify(atlas),before,'Controls and labels preserve source data and camera-independent state');
});

test('coupled marker rejects a suppressed carrier and the wrong shared covering-bank hash',()=>{
 const hidden=coupledFixture();hidden.parts.find(p=>p.id==='TA98TONGUE:male:mucosa').suppressed=true;
 assert.equal(api.tongueRegistry(hidden),null);
 const stale=coupledFixture();stale.finiteTongueIntegration.carrierBankHashes['TA98TONGUE:male:mucosa']='b'.repeat(64);
 assert.equal(api.tongueRegistry(stale),null);
});

test('legacy whole selection requires its exact declared source, native concept and complete owned assembly',()=>{
 for(const damage of [
  a=>delete a.finiteTongueIntegration.sourcePartId,
  a=>a.finiteTongueIntegration.sourcePartId='other',
  a=>a.finiteTongueIntegration.legacyWholeAlias=false,
  a=>a.parts.find(p=>p.id==='ZA:Tongue').suppressed=false,
  a=>a.parts.find(p=>p.id==='ZA:Tongue').conceptId='other',
  a=>a.concepts.find(c=>c.id==='ZA:Tongue').elements.pop(),
  a=>a.concepts.find(c=>c.id==='ZA:Tongue').ta98Term='A05.1.04.101',
 ]){
  const atlas=coupledFixture();damage(atlas);
  const legacy=atlas.concepts.find(c=>c.id==='ZA:Tongue');
  assert.ok(api.tongueRegistry(atlas),'A malformed optional alias does not become geometry');
  assert.equal(api.selectedTongueFeature(atlas,legacy),undefined);
 }
 const atlas=coupledFixture(),legacy=atlas.concepts.find(c=>c.id==='ZA:Tongue');
 assert.equal(api.selectedTongueFeature(atlas,{...legacy,elements:legacy.elements.slice(1)}),undefined);
 const ordinary=fixture(),native=ordinary.concepts.find(c=>c.id===host);
 native.elements=[...ordinary.finiteTongueIntegration.partIds];
 assert.equal(api.selectedTongueFeature(ordinary,native),undefined,'A native whole with matching TA98 term is not implicitly aliased');
});

test('actual staged coupled source without a TA98 field resolves only the guarded whole alias',()=>{
 const file='.local-models/geometry-workbench/tongue-completion/coupled-male-v3/atlas.json';
 if(!fs.existsSync(file))return;
 const atlas=JSON.parse(fs.readFileSync(file)),before=JSON.stringify(atlas),native=atlas.concepts.find(c=>c.id==='ZA:Tongue');
 assert.equal(native.ta98Term,undefined);
 assert.equal(api.selectedTongueFeature(atlas,native)?.key,'tongue');
 assert.equal(api.tongueComponentChoices(atlas,native).length,16);
 assert.equal(JSON.stringify(atlas),before);
});

test('native selection remains native while explicit tongue components include correctly typed references',()=>{
 const a=fixture(),before=structuredClone(a),registry=api.tongueRegistry(a);assert.ok(registry);
 const root=a.concepts.find(c=>c.id===whole),native=a.concepts.find(c=>c.id===host),dorsum=a.concepts.find(c=>c.id.endsWith(':dorsum'));
 assert.equal(api.selectedTongueFeature(a,native),undefined);
 assert.equal(api.tongueComponentChoices(a,root).length,5);
 assert.deepEqual(api.tongueLabelReference(a,dorsum),{partId:host,point:[0,1.6,0]});
 const labels=api.selectedAnatomyLabels(a,root.elements,root,true);
 assert.equal(labels.length,5);assert.equal(labels.filter(l=>l.anchor).length,3);
 const label=labels.find(l=>l.id===dorsum.id);assert.deepEqual(api.anatomyLabelChoice(label),{...dorsum,elements:[host]});
 assert.equal(api.meshCount(a.parts),2);assert.deepEqual(a,before,'References and component labels cannot add geometry or alter source fields');
 const references=api.anatomicalReferenceChoices(a),ref=references.get('A05.1.04.004');assert.ok(ref);
 assert.equal(api.anatomicalReferenceActive(ref,dorsum),true);
 assert.deepEqual(api.toggleAnatomicalReferenceChoice(ref,dorsum),{id:host,name:'Tongue',elements:[host]});
});

test('Systems contains the complete organ while canonical TA98 muscle parentage stays unchanged',()=>{
 const a=fixture(),tree=api.hierarchyNavigation(a,'systems'),visited=[];const walk=c=>{visited.push(c);c.children?.forEach(walk);};walk(tree);
 const organ=visited.find(c=>c.name==='Tongue'&&c.children?.length);assert.ok(organ);
 const local=[];const collect=c=>{local.push(c);c.children?.forEach(collect);};collect(organ);
 assert.ok(local.some(c=>c.name==='Muscles of tongue'));assert.ok(local.some(c=>c.name==='Genioglossus'));assert.ok(local.some(c=>c.name==='Apex of tongue'));
 assert.deepEqual(new Set(organ.elements),new Set([host,muscle]));assert.equal(api.meshCount(organ.elements.map(id=>a.parts.find(p=>p.id===id))),2);
 const definition=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))),roots=api.resolveGuestHierarchy(a,definition),nodes=new Map(),parents=new Map();
 const visit=n=>{if(nodes.has(n.id))return;nodes.set(n.id,n);for(const c of n.children){const set=parents.get(c.id)??new Set();set.add(n.id);parents.set(c.id,set);visit(c);}};roots.forEach(visit);
 assert.ok(parents.get('A05.1.04.101').has('A04.1.05.001'));assert.ok(parents.get('A04.1.05.001').has('A04.1.00.001'));
 assert.deepEqual(nodes.get('A05.1.04.004').directParts.map(p=>p.id),[host]);
 assert.deepEqual(nodes.get('A05.1.04.001').directParts.map(p=>p.id),[host,muscle]);
});

test('stale, invalid and cyclic tongue references are rejected rather than becoming anatomy',()=>{
 for(const damage of [a=>a.finiteTongueIntegration.termGroups[3].surface.faces.push(4),a=>a.finiteTongueIntegration.termGroups[3].surface.sourceBankSha256='b'.repeat(64),a=>a.finiteTongueIntegration.termGroups[4].anchor.barycentric=[1,1,0],a=>a.finiteTongueIntegration.termGroups[1].parent='genioglossus',a=>a.concepts.find(c=>c.id.endsWith(':dorsum')).elements.push(muscle),a=>a.finiteTongueIntegration.termGroups[5].curves[0].points[1][0]=NaN,a=>a.parts[0].suppressed=true]){
  const a=fixture();damage(a);assert.equal(api.tongueRegistry(a),null);
 }
});

test('existing facet selection ink, point landmarks and curves share their carrier without duplicate tissue',()=>{
 const a=fixture(),parent=new api.Group(),geometry=new api.BufferGeometry().setAttribute('position',new api.Float32BufferAttribute([0,0,0,.02,0,0,0,.02,0,0,0,.02],3)).setIndex([0,1,2,0,3,1,0,2,3,1,3,2]);
 const mesh=new api.Mesh(geometry);mesh.matrixWorld.makeTranslation(1,2.6,3);const overlay=api.createSurfaceFeatureOverlay(parent,[]);
 const features=a.finiteTongueIntegration.termGroups;
 overlay.update(features[3],mesh,a.parts[0]);
 const fill=parent.children.find(o=>o.isMesh),line=parent.children.find(o=>o.isLineSegments),point=parent.children.find(o=>o.isPoints);
 assert.equal(fill.visible,true);assert.equal(fill.geometry.getAttribute('position'),geometry.getAttribute('position'));assert.equal(fill.geometry.getIndex().count,3);assert.equal(point.visible,false);
 overlay.update(features[4],mesh,a.parts[0]);assert.equal(fill.visible,false);assert.equal(line.visible,false);assert.equal(point.visible,true);assert.deepEqual([...point.geometry.getAttribute('position').array],[0,0,0]);assert.deepEqual(point.matrix.elements,mesh.matrixWorld.elements);
 overlay.update(features[5],mesh,a.parts[0]);assert.equal(fill.visible,false);assert.equal(line.visible,true);assert.equal(point.visible,false);assert.equal(line.geometry.getAttribute('position').count,4);assert.deepEqual(line.matrix.elements,mesh.matrixWorld.elements);
 overlay.update(undefined,undefined);assert.ok(parent.children.every(o=>!o.visible));overlay.dispose();assert.equal(parent.children.length,0);assert.equal(geometry.getAttribute('position').count,4);assert.equal(a.parts.length,2);geometry.dispose();
});

test('transverse and vertical fibre controls share one carrier and change only shader weights',()=>{
 const a=fixture();
 for(const [term,key] of [['A05.1.04.108','transverse'],['A05.1.04.109','vertical']]){
  const g={key,parent:'muscles',name:key,term,conceptId:`TONGUE:male:${key}`,partIds:[muscle],kind:'exact'};
  a.finiteTongueIntegration.termGroups.push(g);a.concepts.push({id:g.conceptId,name:key,elements:[muscle],ta98Term:term,ta98Kind:'exact'});
 }
 const before=structuredClone(a),t=a.concepts.find(c=>c.ta98Term==='A05.1.04.108'),v=a.concepts.find(c=>c.ta98Term==='A05.1.04.109');
 assert.deepEqual(api.tongueFibreWeights(a,t),[1,0]);assert.deepEqual(api.tongueFibreWeights(a,v),[0,1]);assert.deepEqual(api.tongueFibreWeights(a,undefined),[.5,.5]);
 const appearance={color:[.7,.3,.3],pattern:'fibrous',scale:.0004,relief:.000015,fibres:{directions:[[1,0,0],[0,1,0]],control:'tongue-tv'}};
 const field=api.tissueWallMaterial(appearance),shader={uniforms:{},vertexShader:'#include <begin_vertex>',fragmentShader:'#include <color_fragment>\n#include <normal_fragment_maps>'};
 field.material.onBeforeCompile(shader,{});assert.equal(shader.uniforms.tissueFibreWeights,field.fibreWeights);
 field.fibreWeights.value.set(...api.tongueFibreWeights(a,v));assert.deepEqual(shader.uniforms.tissueFibreWeights.value.toArray(),[0,1]);
 assert.match(shader.fragmentShader,/tissueDirectionalFibre/);assert.match(shader.fragmentShader,/fwidth\(across/);assert.equal(field.material.map,null);
 assert.deepEqual(a,before);assert.equal(api.meshCount(a.parts),2);field.material.dispose();
 assert.throws(()=>api.tissueWallMaterial({...appearance,fibres:{directions:[[0,0,0],[0,1,0]]}}),/Invalid directional/);
});


test('female legacy whole resolves the actual concept when its source part has a separate identity',()=>{
 const atlas=JSON.parse(JSON.stringify(coupledFixture()).replaceAll(':male:',':female:'));atlas.sex='female';
 const source=atlas.parts.find(p=>p.id==='ZA:Tongue'),legacy=atlas.concepts.find(c=>c.id==='ZA:Tongue');
 source.id='LOCAL:female:tongue:1089';source.conceptId='LOCAL:female:tongue';legacy.id=source.conceptId;
 atlas.finiteTongueIntegration.sex='female';atlas.finiteTongueIntegration.sourcePartId=source.id;
 const before=JSON.stringify(atlas);
 assert.equal(api.selectedTongueFeature(atlas,legacy)?.key,'tongue');
 assert.equal(api.tongueComponentChoices(atlas,legacy).length,16);
 assert.equal(api.selectedAnatomyLabels(atlas,legacy.elements,legacy,true).length,16);
 assert.equal(JSON.stringify(atlas),before,'Female alias selection must not mutate anatomy or camera-independent state');
 const unrelated={...legacy,id:'LOCAL:female:other'};
 assert.equal(api.selectedTongueFeature(atlas,unrelated),undefined,'Matching elements alone do not grant a native whole alias');
});
