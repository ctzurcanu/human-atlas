import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/anatomy-label-plan';export * from './app/label-layout';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
for(const sex of ['male','female']){
 const atlas=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 test(`${sex}: assembled stomach has one expandable label; components retain distinct spaced layers`,()=>{
  const choice=atlas.concepts.find(c=>c.name==='Stomach'&&c.elements.some(id=>id.startsWith('GASTRIC:')));
  assert.ok(choice);
  const collapsed=api.selectedAnatomyLabels(atlas,choice.elements,choice,false);
  assert.equal(collapsed.length,1);assert.equal(collapsed[0].name,'Stomach');assert.ok(collapsed[0].components.length>1);
  const expanded=api.selectedAnatomyLabels(atlas,choice.elements,choice,true),targets=api.layerLabelTargets(expanded,atlas.parts);
  assert.ok(expanded.length>=8);assert.equal(targets.size,expanded.length);
  const coordinates=[...targets.values()].map(t=>t.coordinate);
  assert.equal(new Set(coordinates).size,coordinates.length);
  for(const t of targets.values())assert.ok(t.coordinate>=t.min&&t.coordinate<=t.max);
 });
 test(`${sex}: oesophageal, small-bowel and colonic layers use their own organ interval`,()=>{
  for(const pattern of ['OESOPHWALL:','BOWELWALL:jejunum:','CAECOWALL:ascending:','CAECOWALL:rectum:']){
   const parts=atlas.parts.filter(p=>!p.suppressed&&p.id.includes(pattern));assert.ok(parts.length>=7,pattern);
   const items=api.selectedAnatomyLabels(atlas,parts.map(p=>p.id),undefined,true),targets=api.layerLabelTargets(items,atlas.parts);
   assert.ok(targets.size>=7,pattern);
   const t=[...targets.values()][0],family=api.layerLabelFamily(parts[0]);
   const members=atlas.parts.filter(p=>api.layerLabelFamily(p)===family);
   const low=Math.min(...members.map(p=>p.bounds[0][t.axis])),high=Math.max(...members.map(p=>p.bounds[1][t.axis]));
   assert.ok(Math.abs(t.min-(low+(high-low)*.25))<1e-10);assert.ok(Math.abs(t.max-(low+(high-low)*.75))<1e-10);
  }
 });
 test(`${sex}: cystic label arrow has exactly the inspector's fold component`,()=>{
  const duct=atlas.parts.find(p=>p.id===atlas.cysticFoldIntegration.ductPartId),choice=atlas.concepts.find(c=>c.id===duct.conceptId);
  const items=api.selectedAnatomyLabels(atlas,choice.elements,choice,false);
  assert.equal(items.length,1);assert.deepEqual(items[0].components,[`CYSTICFOLD:${sex}:A05.8.02.012`]);
 });
}
test('anchors choose real retained surface points, using central-half positions before outer ends',()=>{
 const samples=Array.from({length:101},(_,i)=>({partId:'wall',point:[i/100,0,0]}));
 const positions=Array.from({length:8},(_,i)=>api.chooseLayerLabelPoint(samples,{axis:0,min:.25,max:.75,coordinate:.25+i/14},[.5,0,0]).point[0]);
 assert.equal(new Set(positions).size,8);assert.ok(positions.every(p=>p>=.25&&p<=.75));
 const kept=samples.filter(s=>s.point[0]>.9),fallback=api.chooseLayerLabelPoint(kept,{axis:0,min:.25,max:.75,coordinate:.5},[.5,0,0]);assert.ok(kept.includes(fallback));
});
test('a selected component stays a component; ad-hoc selections do not acquire an unselected whole',()=>{
 const parts=[{id:'a',conceptId:'a',name:'A',bounds:[[0,0,0],[1,1,1]]},{id:'b',conceptId:'b',name:'B',bounds:[[0,0,0],[1,1,1]]}];
 const atlas={parts,concepts:[{id:'whole',name:'Whole',elements:['a','b']},{id:'a',name:'A',elements:['a']},{id:'b',name:'B',elements:['b']}]};
 assert.deepEqual(api.selectedAnatomyLabels(atlas,['a'],undefined,false).map(i=>i.name),['A']);
 const grouped=api.selectedAnatomyLabels(atlas,['a','b'],undefined,false);assert.equal(grouped.length,1);assert.equal(grouped[0].name,'Whole');
});
test('group arrow reserves room without overflowing standard label pills',()=>{
 const labels=api.layoutAnatomyLabels([{id:'g',x:400,y:150,text:'Very long anatomical group label',components:['a','b']}],{left:5,right:800,top:5,bottom:500},s=>s.length*7);
 assert.equal(labels.length,1);assert.ok(labels[0].right-labels[0].left>=labels[0].label.length*7+40);assert.ok(labels[0].right<=800);
});
test('native layers without generated-wall metadata use explicit parent membership',()=>{
 const parts=['a','b','c'].map(id=>({id,conceptId:id,name:id,bounds:[[0,0,0],[1,4,1]]}));
 const concepts=[{id:'organ',name:'Organ',elements:['a','b','c']},...parts.map((p,i)=>({id:p.id,name:['Tunica serosa','Tela submucosa','Tunica mucosa'][i],elements:[p.id],partOf:'organ'}))];
 const items=api.selectedAnatomyLabels({parts,concepts},parts.map(p=>p.id),undefined,true),targets=api.layerLabelTargets(items,parts,concepts);
 assert.deepEqual([...targets.values()].map(t=>t.coordinate),[1,2,3]);assert.ok([...targets.values()].every(t=>t.axis===1));
});
test('label selection contains exactly the named item, including a partial layer scope',()=>{
 const item={id:'mucosa',name:'Tunica mucosa',members:['body-mucosa','fundus-mucosa'],components:[]};
 const choice=api.anatomyLabelChoice(item);assert.deepEqual(choice.elements,['body-mucosa','fundus-mucosa']);assert.equal(choice.id,'mucosa');
 item.members.push('neck-mucosa');assert.deepEqual(choice.elements,['body-mucosa','fundus-mucosa'],'selection is an independent snapshot');
 const group={id:'organ',name:'Organ',members:['a','b'],components:['b'],choice:{id:'organ',name:'Organ',elements:['a','b'],children:[]}};
 assert.deepEqual(api.anatomyLabelChoice(group).elements,['a','b']);
});
