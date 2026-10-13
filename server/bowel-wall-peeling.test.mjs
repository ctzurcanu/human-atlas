import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/bowel-wall-peeling';export * from './app/guest-hierarchy';export * from './app/guest-choice';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});const {bowelWallPeelSteps,wallPeelDepth,parseGuestHierarchy,resolveGuestHierarchy,guestNodeChoice,partVisible,sectionPartVisible}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const hierarchy=parseGuestHierarchy(JSON.parse(await fs.readFile('public/assets/ta98.json')));
for(const sex of ['female','male'])test(`${sex}: complete bowel peels in eight tissue stages, with fold cells removed together`,async()=>{
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`));
 for(const organ of ['duodenum','jejunum','ileum']){
  const exemplar=atlas.parts.find(p=>p.bowelWallSource&&p.id.includes(`:BOWELWALL:${organ}:`)),source=atlas.parts.find(p=>p.id===exemplar.bowelWallSource);
  const ids=atlas.parts.filter(p=>p.bowelWallSource===exemplar.bowelWallSource).map(p=>p.id),whole=atlas.concepts.find(c=>c.elements.length===12&&c.elements.every(id=>ids.includes(id))),steps=bowelWallPeelSteps(atlas,whole);assert.ok(steps);assert.equal(steps.length,8);assert.deepEqual(steps.map(s=>s.ids.length),[1,1,1,1,2,2,2,2]);
  assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(whole.elements));
  for(let depth=0;depth<=8;depth++){const hidden=steps.slice(0,depth).flatMap(s=>s.ids);assert.equal(wallPeelDepth(steps,hidden),depth);assert.equal(new Set(hidden).size,hidden.length);}
  assert.equal(wallPeelDepth(steps,[steps[4].ids[0]]),0);
  const fold=atlas.concepts.find(c=>c.id===`SIW:${organ}:circular-folds`),f=bowelWallPeelSteps(atlas,fold);assert.equal(f.length,4);assert.ok(f[0].ids[0].endsWith(':circular-folds'));assert.ok(f[3].ids[0].endsWith(':submucosa'));
  assert.equal(bowelWallPeelSteps(atlas,{...whole,id:'selection-set'}),null);assert.equal(bowelWallPeelSteps(atlas,{...whole,elements:whole.elements.slice(1)}),null);
 }
 assert.equal(bowelWallPeelSteps({...atlas,bowelWallPrototype:undefined},atlas.concepts.find(c=>c.elements.length===12)),null);
});
for(const sex of ['female','male'])test(`${sex}: coat aggregates peel exactly their stored tissues and the hierarchy fold core`,async()=>{
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`));
 const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 for(const organ of ['duodenum','jejunum','ileum'])for(const key of ['muscularis','mucosa']){
  const concept=atlas.concepts.find(c=>c.id===`SIW:${organ}:${key}`),steps=bowelWallPeelSteps(atlas,concept);
  assert.deepEqual(steps.map(s=>s.ids.length),key==='muscularis'?[1,1]:[2,2,2]);assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(concept.elements));
  const choice=guestNodeChoice('ta98',nodes.get(concept.id)),hierarchySteps=bowelWallPeelSteps(atlas,choice);assert.ok(hierarchySteps);
  assert.deepEqual(hierarchySteps.map(s=>s.ids.length),key==='muscularis'?[1,1]:[1,2,2,2]);assert.deepEqual(new Set(hierarchySteps.flatMap(s=>s.ids)),new Set(choice.elements));
  if(key==='mucosa')assert.ok(hierarchySteps[0].ids[0].endsWith(':circular-folds:submucosa'));
  assert.equal(bowelWallPeelSteps(atlas,{...choice,elements:choice.elements.slice(1)}),null);
 }
});

for(const sex of ['female','male'])test(`${sex}: peeled selected tissues actually disappear from surfaces and sections at all eight stages`,async()=>{
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 for(const organ of ['duodenum','jejunum','ileum']){
  const pieces=atlas.parts.filter(p=>p.bowelWallSource&&p.id.includes(`:BOWELWALL:${organ}:`)),ids=pieces.map(p=>p.id),whole=atlas.concepts.find(c=>c.elements.length===12&&c.elements.every(id=>ids.includes(id))),steps=bowelWallPeelSteps(atlas,whole);assert.equal(steps.length,8);
  for(const isolate of [false,true])for(let depth=0;depth<=8;depth++){
   const hidden=steps.slice(0,depth).flatMap(s=>s.ids),state={selected:ids,visible:['digestive'],hidden,depthHidden:[],isolate,region:'all',skinOpacity:0,contextOpacity:0,explode:0};
   assert.deepEqual(new Set(pieces.filter(p=>partVisible(p,state,atlas.parts)).map(p=>p.id)),new Set(ids.filter(id=>!hidden.includes(id))));
   assert.deepEqual(new Set(pieces.filter(p=>sectionPartVisible(p,state)).map(p=>p.id)),new Set(ids.filter(id=>!hidden.includes(id))));
   assert.deepEqual(state.selected,ids,'selection retains the whole assembly and peel controls');
  }
 }
 // Explicit selection still reveals ordinary source anatomy, preserving its
 // existing behaviour. The exception is limited to generated peelable walls.
 const ordinary=atlas.parts.find(p=>!p.bowelWallSource&&!p.suppressed),state={selected:[ordinary.id],visible:[],hidden:[ordinary.id],depthHidden:[],isolate:true,region:'all'};assert.equal(partVisible(ordinary,state),true);assert.equal(sectionPartVisible(ordinary,state),true);
});

test('Ordinary female TA98 jejunum selects its complete peelable wall without an alternate donor exterior',async()=>{
 const atlas=JSON.parse(await fs.readFile('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
 const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 const node=nodes.get('A05.6.03.001'),choice=guestNodeChoice('ta98',node),steps=bowelWallPeelSteps(atlas,choice);
 assert.equal(choice.elements.length,12);assert.equal(steps?.length,8);
 for(const id of ['HRA:VH_F_jejenum','HRA:VH_F_ileum']){const p=atlas.parts.find(p=>p.id===id);assert.equal(p.suppressed,true);assert.equal(p.replacementParts,undefined,'alternate geometry is not source lineage');}
 const source=JSON.parse(await fs.readFile('.local-models/ta98-female-runtime.json'));assert.ok(!source.parts.find(p=>p.id==='HRA:VH_F_jejenum').suppressed,'source comparison remains available');
});

for(const sex of ['female','male'])test(`${sex}: every TA98 bowel organ exposes a complete wall separately from regional and adjoining structures`,async()=>{
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 for(const [organ,code] of [['duodenum','A05.6.02.001'],['jejunum','A05.6.03.001'],['ileum','A05.6.04.001']]){
  const owner=nodes.get(code),wall=owner.children.find(n=>n.id===`SIW:${organ}:wall-assembly`);assert.ok(wall);
  const choice=guestNodeChoice('ta98',wall),steps=bowelWallPeelSteps(atlas,choice);assert.equal(choice.elements.length,12);assert.equal(steps?.length,8);
  assert.ok(wall.parts.every(p=>p.bowelWallSource&&p.id.includes(`:BOWELWALL:${organ}:`)));
  assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(choice.elements));
  for(const part of owner.parts.filter(p=>!p.bowelWallSource))assert.ok(!choice.elements.includes(part.id),'regional or adjoining geometry stays outside peel assembly');
 }
 const source=JSON.parse(await fs.readFile(sex==='male'?'.local-models/ta98-runtime.json':'.local-models/ta98-female-runtime.json'));
 const sourceNodes=new Map();const visit=n=>{sourceNodes.set(n.id,n);n.children.forEach(visit);};resolveGuestHierarchy(source,hierarchy).forEach(visit);
 assert.ok(!sourceNodes.has('SIW:duodenum:wall-assembly'),'no wall assembly is claimed without stored layers');
});
