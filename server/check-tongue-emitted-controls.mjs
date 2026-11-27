// Verify the actual packed catalogue through the same APIs used by the viewer.
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import fs from 'node:fs';
import path from 'node:path';
import {build} from 'esbuild';

const args=process.argv.slice(2),options={};
for(let i=0;i<args.length;i+=2){assert.ok(['--catalogue','--report'].includes(args[i])&&args[i+1],'Use --catalogue atlas.json [--report new.json]');options[args[i]]=args[i+1];}
assert.ok(options['--catalogue'],'An explicit emitted --catalogue is required');
const cataloguePath=path.resolve(options['--catalogue']),root=process.cwd();
assert.ok(cataloguePath.startsWith(path.join(root,'.local-models/geometry-workbench/tongue-completion/')),'Check a local staged tongue catalogue');
const sha=data=>createHash('sha256').update(data).digest('hex');
const raw=fs.readFileSync(cataloguePath),atlas=JSON.parse(raw),before=sha(JSON.stringify(atlas));
const compiled=await build({stdin:{contents:`export * from './app/tongue-controls';export * from './app/hierarchy-navigation';export * from './app/anatomy-label-plan';export * from './app/guest-hierarchy';export * from './app/organ-models';`,resolveDir:root},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const registry=api.tongueRegistry(atlas);assert.ok(registry,'Actual emitted tongue marker must pass frontend validation');
const {marker}=registry,byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
const coupled=marker.legacyWholeAlias===true,sex=marker.sex??atlas.sex;
assert.ok(['male','female'].includes(sex));
const femalePartial=sex==='female',femaleCoreOnly=femalePartial&&marker.partialMode==='female-core-only';
if(femalePartial){
 assert.equal(coupled,true,'Female partial requires its explicit native whole alias');assert.equal(marker.sex,'female');
 if(femaleCoreOnly){
  assert.equal(marker.activePhysicalMeshes,6,'Core-only female revision declares six physical meshes');
  assert.equal(marker.representedAdultControls,14);assert.equal(marker.missingAdultTerms,21);
  assert.equal(marker.unfinishedExtrinsicSourcesSuppressed,true);
 }else{
  assert.equal(marker.partialMode,undefined,'Unsupported female partial mode');
  assert.equal(marker.activePhysicalMeshes,11,'Diagnostic native-extrinsic female partial declares eleven physical meshes');
 }
}else assert.notEqual(marker.partialMode,'female-core-only','Core-only mode is female-specific');
const controlCount=femaleCoreOnly?14:coupled?17:15,expectedMeshCount=femaleCoreOnly?6:femalePartial?11:coupled?12:11;
const meshCount=marker.activePhysicalMeshes??expectedMeshCount,newMeshCount=femalePartial?6:coupled?8:6;
assert.equal(meshCount,expectedMeshCount,'Declared physical count must match this supported partial revision');
const nativeConceptId=femalePartial?'LOCAL:female:tongue':'ZA:Tongue',nativePartId=femalePartial?'LOCAL:female:tongue:1089':'ZA:Tongue';
const rootChoice=concepts.get(marker.wholeConceptId),native=concepts.get(nativeConceptId);
assert.ok(rootChoice&&native);assert.equal(marker.termGroups.length,controlCount);assert.equal(marker.partIds.length,meshCount);assert.equal(api.meshCount(registry.parts),meshCount);
if(coupled){assert.equal(marker.sourcePartId,nativePartId);assert.deepEqual(native.elements,marker.partIds);assert.equal(byId.get(nativePartId).suppressed,true);assert.ok(!marker.partIds.includes(nativePartId));assert.equal(api.selectedTongueFeature(atlas,native)?.key,'tongue');}
else{assert.deepEqual(native.elements,[nativePartId]);assert.equal(api.selectedTongueFeature(atlas,native),undefined,'Native source selection must remain native');}
assert.equal(api.tongueComponentChoices(atlas,rootChoice).length,controlCount-1);
const labels=api.selectedAnatomyLabels(atlas,rootChoice.elements,rootChoice,true);
assert.equal(labels.length,controlCount-1);assert.equal(new Set(labels.map(l=>l.id)).size,controlCount-1);
if(coupled){
 assert.equal(api.tongueComponentChoices(atlas,native).length,controlCount-1,'Native Tongue selection exposes every available component');
 const legacyLabels=api.selectedAnatomyLabels(atlas,native.elements,native,true);
 assert.equal(legacyLabels.length,controlCount-1);assert.deepEqual(new Set(legacyLabels.map(l=>l.id)),new Set(labels.map(l=>l.id)),'Native and canonical whole labels are identical');
 assert.ok(native.elements.every(id=>id!==nativePartId&&!byId.get(id).suppressed),'Native whole selection cannot expose the old surface');
}
for(const label of labels){const choice=api.anatomyLabelChoice(label);assert.equal(choice.id,label.id);assert.deepEqual(choice.elements,concepts.get(label.id).elements,'Clicking a label selects that specific control');}
const features=marker.termGroups.filter(g=>g.representation);assert.equal(features.length,6);assert.equal(labels.filter(l=>l.anchor).length,6);
const boundBanks=new Map();
const bankFor=part=>{
 const url=atlas.chunks[part.chunk].url;assert.ok(url.startsWith('/local-models/'));
 const file=path.join(root,'.local-models',url.slice('/local-models/'.length));
 if(!boundBanks.has(file)){const data=fs.readFileSync(file);boundBanks.set(file,{data,sha256:sha(data)});}
 return {file,...boundBanks.get(file)};
};
const emittedPosition=(p,view,vertex)=>{
 const scale=p.geometryScale??1,origin=p.geometryOrigin??[0,0,0],offset=p.sourceOffset??[0,0,0];
 assert.ok(Number.isFinite(scale)&&scale>0);assert.ok(origin.length===3&&offset.length===3&&[...origin,...offset].every(Number.isFinite));
 const local=Array.from({length:3},(_,axis)=>view.getFloat32(p.positions+(vertex*3+axis)*4,true)*scale+offset[axis]);
 const instance=p.geometryInstance;
 if(!instance)return local.map((v,i)=>v+origin[i]);
 const q=instance.rotation;assert.ok(Array.isArray(q)&&q.length===4&&q.every(Number.isFinite));
 const length=Math.hypot(...q);assert.ok(Math.abs(length-1)<1e-8&&Number.isFinite(instance.scale)&&instance.scale>0,'Only a proper positive uniform geometry instance is allowed');
 const [x,y,z,w]=q.map(v=>v/length),r=[
  [1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
  [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
  [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)],
 ];
 // Match AtlasReader/viewer order: rotate the scaled bank coordinate plus
 // local source offset, apply instance scale once, then add world origin.
 return r.map((row,i)=>row.reduce((n,v,j)=>n+v*local[j],0)*instance.scale+origin[i]);
};
for(const g of features){
 const p=byId.get(g.anchor.partId),bank=bankFor(p),view=new DataView(bank.data.buffer,bank.data.byteOffset,bank.data.byteLength);
 assert.equal(bank.sha256,marker.carrierBankHashes[g.anchor.partId]);
 const anchor=[0,0,0];
 for(let c=0;c<3;c++){
  const vertex=view.getUint32(p.indices+(g.anchor.face*3+c)*4,true);
  const xyz=emittedPosition(p,view,vertex);for(let axis=0;axis<3;axis++)anchor[axis]+=xyz[axis]*g.anchor.barycentric[c];
 }
 assert.ok(Math.hypot(...anchor.map((v,i)=>v-g.anchor.point[i]))<1e-10,'Reference anchor must match its actual emitted source triangle');
 const choice=concepts.get(g.conceptId);assert.equal(api.selectedTongueSurfaceFeature(atlas,choice)?.key,g.key);assert.equal(api.tongueLabelReference(atlas,choice).partId,p.id);
}
const tv=marker.termGroups.filter(g=>['A05.1.04.108','A05.1.04.109'].includes(g.term));assert.equal(tv.length,2);assert.deepEqual(tv[0].partIds,tv[1].partIds);assert.equal(tv[0].partIds.length,1);
for(const g of tv)assert.deepEqual(api.tongueFibreWeights(atlas,concepts.get(g.conceptId)),g.term.endsWith('108')?[1,0]:[0,1]);
const tvPart=byId.get(tv[0].partIds[0]),tvMaterial=atlas.materials[tvPart.material];
assert.equal(tvMaterial.tissueAppearance.fibres.control,'tongue-tv');assert.equal(tvMaterial.tissueAppearance.relief,.000015);assert.deepEqual(tvMaterial.tissueAppearance.fibres.directions,[[1,0,0],[0,1,0]]);
const prefix='TA98TONGUE:'+sex+':',mucosaId=prefix+'mucosa';
const newIds=atlas.parts.filter(p=>p.id.startsWith(prefix+'intrinsic:')||p.id.startsWith(prefix+'styloglossus:')||coupled&&[mucosaId,prefix+'aponeurosis'].includes(p.id));
assert.equal(newIds.length,newMeshCount);
const retainedNativeIds=['LOCAL:female:genioglossus:239','LOCAL:female:l_hyoglossus:242','LOCAL:female:r_hyoglossus:255','LOCAL:female:l_styloglossus:243','LOCAL:female:r_styloglossus:256'];
if(femalePartial){
 const expectedNew=[mucosaId,prefix+'aponeurosis',...['superior-longitudinal','transverse-vertical-shared','inferior-longitudinal-left','inferior-longitudinal-right'].map(key=>prefix+'intrinsic:'+key)];
 assert.deepEqual(new Set(newIds.map(p=>p.id)),new Set(expectedNew),'Female partial contains exactly its six transferred core/IL meshes');
 const allTerms=['A04.1.05.001',...['001','004','005','006','007','009','010','011','025','101','102','105','106','107','108','109'].map(t=>'A05.1.04.'+t)],extrinsicTerms=['A05.1.04.101','A05.1.04.102','A05.1.04.105'];
 assert.deepEqual(new Set(marker.termGroups.map(g=>g.term)),new Set(allTerms.filter(term=>!femaleCoreOnly||!extrinsicTerms.includes(term))),'Female controls must match the explicit supported term inventory');
 assert.deepEqual(new Set(marker.partIds),new Set(femaleCoreOnly?expectedNew:[...expectedNew,...retainedNativeIds]),'Only the declared revision contributes physical members');
 for(const id of retainedNativeIds){
  const p=byId.get(id);assert.ok(p,'Retained native work geometry must remain available');
  assert.equal(p.geometryInstance,undefined,'Retained native extrinsics are not transferred instances');
  if(femaleCoreOnly){assert.equal(p.suppressed,true,'Unfinished native extrinsics remain suppressed');assert.ok(atlas.defaultHidden?.includes(id));assert.ok(!marker.partIds.includes(id));}
  else assert.ok(!p.suppressed&&!atlas.defaultHidden?.includes(id),'Diagnostic native extrinsics remain active');
 }
 if(femaleCoreOnly){
  assert.deepEqual(new Set(marker.suppressedRetainedWorkPartIds),new Set(retainedNativeIds));
  assert.ok(!marker.termGroups.some(g=>['genioglossus','hyoglossus','styloglossus'].includes(g.key)));
  assert.ok(!labels.some(label=>extrinsicTerms.includes(concepts.get(label.id)?.ta98Term)),'Unfinished extrinsics cannot appear as component labels');
  assert.equal(marker.termGroups.find(g=>g.key==='muscles').partIds.length,4);
 }else{
  assert.deepEqual(marker.termGroups.find(g=>g.key==='genioglossus').partIds,[retainedNativeIds[0]]);
  assert.deepEqual(new Set(marker.termGroups.find(g=>g.key==='hyoglossus').partIds),new Set(retainedNativeIds.slice(1,3)));
  assert.deepEqual(new Set(marker.termGroups.find(g=>g.key==='styloglossus').partIds),new Set(retainedNativeIds.slice(3)));
  assert.equal(marker.termGroups.find(g=>g.key==='muscles').partIds.length,9);
 }
 assert.equal(byId.get('HRA:VH_F_dorsal_tongue').suppressed,true);assert.ok(atlas.defaultHidden.includes('HRA:VH_F_dorsal_tongue'));assert.ok(!marker.partIds.includes('HRA:VH_F_dorsal_tongue'));
}
for(const p of newIds){
 assert.equal(p.uvs,undefined);
 const inShared=coupled&&(marker.sharedBankPartIds??[mucosaId,prefix+'aponeurosis']).includes(p.id);
 const expected=femalePartial?marker.carrierBankHashes?.[p.id]:inShared?bankFor(byId.get(mucosaId)).sha256:marker.bankSha256;
 assert.match(expected??'',/^[a-f0-9]{64}$/,'Every transferred bank has an explicit source binding');assert.equal(bankFor(p).sha256,expected);assert.ok(p.license&&p.provenance);
}
const nav=api.hierarchyNavigation(atlas,'systems'),systemNodes=[];
const walk=(node,seen=new Set())=>{if(seen.has(node))return;seen.add(node);systemNodes.push(node);node.children?.forEach(c=>walk(c,seen));};walk(nav);
const organ=systemNodes.find(n=>n.name==='Tongue'&&n.children?.length&&new Set(n.elements).size===meshCount&&n.elements.every(id=>marker.partIds.includes(id)));
assert.ok(organ,'Systems must expose one complete Tongue assembly');
const descendants=[],visitSystem=node=>{descendants.push(node);node.children?.forEach(visitSystem);};visitSystem(organ);
for(const g of marker.termGroups)assert.ok(descendants.some(n=>n.name===g.name&&g.partIds.every(id=>n.elements.includes(id))),'Systems is missing '+g.name);
assert.equal(api.meshCount(organ.elements.map(id=>byId.get(id))),meshCount,'Systems total must deduplicate surface aliases and shared T/V carrier');
const definition=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json'))),roots=api.resolveGuestHierarchy(atlas,definition),nodes=new Map(),parents=new Map();
const visit=node=>{if(nodes.has(node.id))return;nodes.set(node.id,node);for(const child of node.children){const p=parents.get(child.id)??new Set();p.add(node.id);parents.set(child.id,p);visit(child);}};roots.forEach(visit);
for(const g of marker.termGroups)assert.deepEqual(new Set(nodes.get(g.term)?.directParts.map(p=>p.id)),new Set(g.partIds),'TA98 direct controls are missing or stale for '+g.term);
if(femaleCoreOnly)for(const term of ['A05.1.04.101','A05.1.04.102','A05.1.04.105'])assert.equal(nodes.get(term)?.directParts.length??0,0,'Unfinished native extrinsics do not count as represented TA98 controls');
assert.ok(parents.get('A05.1.04.101').has('A04.1.05.001'));assert.ok(parents.get('A04.1.05.001').has('A04.1.00.001'),'Canonical muscle parents must not be rewritten by Systems grouping');
assert.equal(api.meshCount(nodes.get('A05.1.04.001').parts),meshCount,'TA98 Tongue total must count physical meshes only');
assert.equal(sha(JSON.stringify(atlas)),before,'Hierarchy, label and fibre controls cannot mutate the catalogue');
const result={schema:'tongue-emitted-frontend-controls-check/v1',catalogue:cataloguePath,catalogueSha256:sha(raw),sex,partialMode:femaleCoreOnly?'female-core-only':undefined,controls:controlCount,expandedLabels:controlCount-1,typedSurfaceReferences:6,actualMeshes:meshCount,newMeshes:newMeshCount,retainedNativeMeshes:femalePartial&&!femaleCoreOnly?5:undefined,suppressedRetainedWorkMeshes:femaleCoreOnly?5:undefined,missingAdultTerms:femaleCoreOnly?21:undefined,sharedTVPhysicalMeshes:1,systemsControlsPresent:true,canonicalTA98ParentsPreserved:true,actualBankHashes:[...boundBanks].map(([file,r])=>({path:file,sha256:r.sha256})),sourceDataUnchanged:true,passed:true,anatomicalAcceptance:false};
if(options['--report']){const report=path.resolve(options['--report']);assert.ok(!fs.existsSync(report),'Use a new report file');fs.writeFileSync(report,JSON.stringify(result,null,2)+'\n',{flag:'wx'});}
console.log(JSON.stringify(result));
