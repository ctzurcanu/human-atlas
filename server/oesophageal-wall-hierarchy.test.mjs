import test from 'node:test';
import assert from 'node:assert/strict';
import {existsSync,readFileSync} from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/guest-hierarchy'; export * from './app/bowel-wall-peeling'; export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {resolveGuestHierarchy,bowelWallPeelSteps,anatomicalConceptFunctionalRelations}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const path='.local-models/ta98-review/oesophageal-layer-bank/male/runtime-package/package.json';
test('TA98 oesophagus presents every stored coat, with mucosa and muscularis subtrees',{skip:!existsSync(path)},()=>{
 const packed=JSON.parse(readFileSync(path,'utf8'));
 const keys=packed.parts.slice(0,7).map(p=>p.tissueWallLayer);
 const parts=packed.parts.map(p=>({...p,chunk:0,depth:8}));
 const members=keys=>parts.filter(p=>keys.includes(p.tissueWallLayer)).map(p=>p.id);
 const concepts=[{id:'ZA:Oesophagus',name:'Oesophagus',elements:parts.map(p=>p.id)}];
 for(const key of ['wall-assembly','muscularis','mucosa',...keys]){
  const selected=key==='wall-assembly'?keys:key==='muscularis'?keys.slice(1,3):key==='mucosa'?keys.slice(4):[key];
  concepts.push({id:'OESOPHWALL:male:'+key,name:key,elements:members(selected),...(key==='submucosa'?{ta98Term:'A05.4.01.014',ta98Kind:'exact'}:{})});
 }
 const atlas={parts,concepts,bowelWallPrototype:{partIds:parts.map(p=>p.id)}};
 const hierarchy={schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98',roots:['A05.4.01.001'],nodes:[{id:'A05.4.01.001',name:'oesophagus',matches:['ZA:Oesophagus'],childIds:['A05.4.01.014']},{id:'A05.4.01.014',name:'submucosa'}]};
 const [organ]=resolveGuestHierarchy(atlas,hierarchy);assert.equal(organ.parts.length,14);
 const wall=organ.children.find(n=>n.id==='OESOPHWALL:male:wall-assembly');assert.ok(wall);assert.equal(wall.parts.length,14);
 const muscle=wall.children.find(n=>n.id.endsWith(':muscularis'));assert.equal(muscle.children.length,2);assert.equal(muscle.parts.length,4);
 const mucosa=wall.children.find(n=>n.id.endsWith(':mucosa'));assert.equal(mucosa.children.length,3);assert.equal(mucosa.parts.length,6);
 assert.equal(organ.children.find(n=>n.id==='A05.4.01.014').parts.length,2);
 for(const [node,count] of [[wall,7],[muscle,2],[mucosa,3]])assert.equal(bowelWallPeelSteps(atlas,{id:node.id,elements:node.parts.map(p=>p.id)}).length,count);
});

test('complete oesophageal assembly retains organ-level supply and innervation without assigning them to each coat',{skip:!existsSync(path)},()=>{
 const packed=JSON.parse(readFileSync(path,'utf8'));const selected=packed.parts.map(p=>({...p,chunk:0}));
 const targets=[{id:'ZA:Vagus nerve (X).l',name:'Vagus nerve (X) (left)',system:'nervous'},{id:'ZA:Left gastric artery',name:'Left gastric artery',system:'arterial'},{id:'ZA:Left gastric vein',name:'Left gastric vein',system:'venous'}];
 const concept={id:'ZA:Oesophagus',name:'Oesophagus',elements:selected.map(p=>p.id)};
 const relations=anatomicalConceptFunctionalRelations(concept,selected,[...selected,...targets]);
 for(const kind of ['innervation','arterial','venous'])assert.ok(relations.some(r=>r.kind===kind),kind);
 assert.equal(anatomicalConceptFunctionalRelations(concept,selected.slice(1),[...selected,...targets]).length,0);
});
