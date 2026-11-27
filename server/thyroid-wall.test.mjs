import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
import {existsSync,readFileSync} from 'node:fs';

const modules=['thyroid-wall','bowel-wall-peeling','viewer-state','anatomy-label-plan','preferred-source-representations'];
const compiled=await build({stdin:{contents:modules.map(name=>`export * from './app/${name}';`).join('\n'),resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));

// Metadata-only control fixture. Anatomical shape remains a separate gate.
function fixture(sex='male'){
 const prefix=`THYROID:${sex}:`,whole=prefix+'whole',regionNames=['right-lobe','left-lobe','isthmus','pyramidal-lobe'];
 const capsuleParts=regionNames.map(region=>prefix+region+':capsule'),parenchymaParts=regionNames.map(region=>prefix+region+':parenchyma');
 const appearance={color:[.64,.39,.4],pattern:'fibrous',scale:.0008,relief:.000002};
 const cells=['capsule','parenchyma'].flatMap((role,k)=>regionNames.map((region,r)=>({
  id:(k?parenchymaParts:capsuleParts)[r],conceptId:(k?parenchymaParts:capsuleParts)[r],name:`${role} of ${region}`,system:'endocrine',
  chunk:0,positions:0,normals:192,vertexCount:16,indexCount:36,indices:384+(k*4+r)*144,
  normalSigns:1600+k*64+r*16,geometryOrigin:[0,1.55,0],bounds:[[-.03,1.52,0],[.03,1.58,.03]],
  pairedTissue:true,sectionSolid:true,tissueWallAssembly:`thyroid:${sex}`,tissueWallLayer:role,tissueAppearance:appearance,
  tissueDrawGroups:[{start:0,count:12},{start:12,count:12,neighbor:(k?capsuleParts:parenchymaParts)[r]},
   {start:24,count:12,neighbor:(k?parenchymaParts:capsuleParts)[(r+2)%4]}]
 })));
 const marker={schema:'ta98-thyroid-wall-v1',sex,wholeConcept:whole,regionConcepts:regionNames.map(region=>prefix+region),
  capsuleConcept:prefix+'capsule',parenchymaConcept:prefix+'parenchyma',capsuleParts,parenchymaParts,sourceReview:{status:'provisional'}};
 const concepts=[{id:whole,name:'Thyroid gland',elements:cells.map(p=>p.id),ta98Term:'A11.3.00.001',ta98Kind:'exact'},
  ...regionNames.map((region,r)=>({id:marker.regionConcepts[r],name:region,elements:[capsuleParts[r],parenchymaParts[r]],partOf:whole,
   ta98Term:r<2?'A11.3.00.002':r===2?'A11.3.00.003':'A11.3.00.004',ta98Kind:'parent'})),
  {id:marker.capsuleConcept,name:'Fibrous capsule',elements:capsuleParts,partOf:whole,ta98Term:'A11.3.00.006',ta98Kind:'exact'},
  {id:marker.parenchymaConcept,name:'Parenchyma',elements:parenchymaParts,partOf:whole,ta98Term:'A11.3.00.008',ta98Kind:'exact'},
  ...cells.map(p=>({id:p.conceptId,name:p.name,elements:[p.id],partOf:whole}))];
 return {version:'control-test',sex,thyroidWallAssembly:marker,parts:cells,concepts,chunks:[],triangles:0};
}
const whole=a=>a.concepts.find(c=>c.id===a.thyroidWallAssembly.wholeConcept);

const installedFemalePath='.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';
test('Installed female thyroid is visible and peelable after ordinary URL restoration',
 {skip:!existsSync(installedFemalePath)},()=>{
 const source=JSON.parse(readFileSync(installedFemalePath,'utf8'));
 if(!source.femaleThyroidInstallation)return;
 const a=api.withHiddenReviewAlternatives(source),choice=whole(a),registry=api.thyroidWallRegistry(a);
 assert.ok(registry,'Local installation status must preserve the finite tissue registry');
 assert.equal(registry.members.length,11);
 const initial={selected:[],hidden:[],visible:['endocrine'],isolate:false,view:'front',rotate:false,
  reset:0,explode:0,skinOpacity:0,context:0,region:'all'};
 const restored=api.readViewUrl('http://localhost:3016/?select=LOCAL%3Afemale%3Athyroid_gland&skin=0&context=0&isolate=1',a,initial);
 assert.deepEqual(new Set(restored.selected),new Set(choice.elements));
 assert.ok(registry.members.every(p=>api.partVisible(p,restored,a.parts)),
  'A hidden retired source must not hide its selected replacement tissues');
 const steps=api.bowelWallPeelSteps(a,choice);
 assert.deepEqual(steps.map(step=>step.ids.length),[4,7]);
 const peeled={...restored,hidden:[...restored.hidden,...steps[0].ids]};
 assert.equal(registry.members.filter(p=>api.partVisible(p,peeled,a.parts)).length,7);
 assert.ok(a.parts.find(p=>p.id==='LOCAL:female:thyroid_gland:1081').suppressed);
});

for(const sex of ['male','female'])test(`${sex}: whole and complete anatomical regions peel capsule before internal gland tissue`,()=>{
 const a=fixture(sex),r=api.thyroidWallRegistry(a);assert.ok(r);assert.equal(r.marker.sourceReview.status,'provisional');
 assert.equal(a.bowelWallPrototype,undefined,'Thyroid controls do not depend on bowel registry');
 for(const choice of [whole(a),...r.regions]){
  const steps=api.bowelWallPeelSteps(a,choice);assert.deepEqual(steps.map(s=>s.name),['Fibrous capsule','Parenchyma']);
  assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(choice.elements));
  for(let depth=0;depth<=2;depth++){
   const hidden=steps.slice(0,depth).flatMap(s=>s.ids),view={selected:choice.elements,hidden,visible:['endocrine'],isolate:true,view:'front',rotate:false,reset:0,explode:0};
   assert.equal(api.wallPeelDepth(steps,hidden),depth);
   for(const part of a.parts.filter(p=>choice.elements.includes(p.id))){assert.equal(api.partVisible(part,view,a.parts),!hidden.includes(part.id));assert.equal(api.sectionPartVisible(part,view),!hidden.includes(part.id));}
  }
 }
 for(const choice of [r.capsule,r.parenchyma])assert.equal(api.bowelWallPeelSteps(a,choice).length,1);
 assert.equal(api.bowelWallPeelSteps(a,{...whole(a),id:'selection-set'}),null);
 assert.equal(api.bowelWallPeelSteps(a,{...whole(a),elements:whole(a).elements.slice(1)}),null);
 const mixed=[...r.regions[0].elements,...r.regions[2].elements];assert.equal(api.bowelWallPeelSteps(a,{id:'mixed',name:'Mixed tissues',elements:mixed}),null);
 const lobes=r.regions.filter(c=>c.ta98Term==='A11.3.00.002').flatMap(c=>c.elements);
 const collective={id:'hierarchy:guest:ta98:A11.3.00.002',name:'Lobe',elements:lobes};
 assert.deepEqual(api.bowelWallPeelSteps(a,collective).map(s=>s.ids.length),[2,2]);
 assert.equal(api.bowelWallPeelSteps(a,{...collective,elements:[...r.regions[0].elements,...r.regions[2].elements]}),null);
});

test('Unregistered exteriors, wrong roles, incomplete contacts and independent coordinate frames cannot activate thyroid peeling',()=>{
 for(const change of [{pairedTissue:false},{sectionSolid:false},{suppressed:true},{sourceOffset:[0,0,0]},{system:'connective'},
  {tissueWallAssembly:'thyroid:female'},{tissueWallLayer:'parenchyma'},{chunk:3},{positions:999},{normals:999},{vertexCount:12},
  {geometryOrigin:[0,0,0]},{normalSigns:undefined},{tissueAppearance:undefined},{tissueDrawGroups:[]}]){
  const a=fixture();a.parts[0]={...a.parts[0],...change};assert.equal(api.thyroidWallRegistry(a),null,JSON.stringify(change));assert.equal(api.bowelWallPeelSteps(a,whole(a)),null);
 }
 for(const mutate of [a=>delete a.thyroidWallAssembly,a=>a.sex='female',a=>a.thyroidWallAssembly.schema='unknown',
  a=>a.thyroidWallAssembly.sourceReview.status='accepted',a=>a.thyroidWallAssembly.capsuleParts.push(a.thyroidWallAssembly.parenchymaParts[0]),
  a=>a.parts[0].tissueDrawGroups[1].count=9,a=>a.parts[0].tissueDrawGroups[1].neighbor='other-organ',
  a=>a.parts[0].tissueDrawGroups[0].start=3,a=>a.parts[0].indexCount=39,
  a=>a.parts.push({...a.parts[0],id:'unexpected-thyroid-cell'}),
  a=>a.concepts.find(c=>c.id===a.thyroidWallAssembly.capsuleConcept).partOf='other-organ',
  a=>a.concepts.find(c=>c.id===a.thyroidWallAssembly.parenchymaConcept).ta98Term='A11.3.00.007',
  a=>a.concepts.find(c=>c.id===a.thyroidWallAssembly.parenchymaConcept).ta98Kind='parent',
  a=>whole(a).elements.pop(),a=>a.concepts[1].elements.push(a.concepts[2].elements[0]),
  a=>a.thyroidWallAssembly.stromaConcept='unbuilt-stroma',a=>a.thyroidWallAssembly.lobuleParts=['unbuilt-lobule']]){
  const a=fixture();mutate(a);assert.equal(api.thyroidWallRegistry(a),null);
 }
});

test('Registered tissues retain ordinary labels and explicit selection camera frame',()=>{
 const a=fixture(),choice=whole(a),before=JSON.stringify(a),collapsed=api.selectedAnatomyLabels(a,choice.elements,choice,false);
 assert.equal(collapsed.length,1);assert.equal(collapsed[0].name,'Thyroid gland');assert.equal(collapsed[0].components.length,8);
 const expanded=api.selectedAnatomyLabels(a,choice.elements,choice,true),targets=api.layerLabelTargets(expanded,a.parts,a.concepts);
 assert.equal(expanded.length,8);assert.equal(targets.size,8);assert.equal(new Set([...targets.values()].map(t=>t.coordinate)).size,8);
 for(const item of expanded)assert.deepEqual(api.anatomyLabelChoice(item).elements,item.members);
 const frame=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5];
 const view=api.readViewUrl(new URLSearchParams({select:choice.id,frame:frame.join(',')}).toString(),a,{selected:[],visible:[],view:'front'});
 assert.deepEqual(view.camera,frame);assert.deepEqual(view.selected,choice.elements);assert.equal(JSON.stringify(a),before,'Control resolution does not mutate source catalogue or geometry');
});

test('Lobules select existing parenchymal domains without adding another part or peel coat',()=>{
 const a=fixture(),before=JSON.stringify(a.parts),partCount=a.parts.length,indexCount=a.parts.reduce((n,p)=>n+p.indexCount,0);
 const marker=a.thyroidWallAssembly;marker.lobuleParts=marker.parenchymaParts.slice(0,2);marker.lobulesConcept='THYROID:male:lobules';
 a.concepts.push({id:marker.lobulesConcept,name:'Lobules',elements:[...marker.lobuleParts],ta98Term:'A11.3.00.009',ta98Kind:'exact',partOf:marker.wholeConcept});
 const registry=api.thyroidWallRegistry(a);assert.ok(registry);assert.equal(registry.lobules.id,marker.lobulesConcept);
 assert.equal(a.parts.length,partCount);assert.equal(registry.members.length,partCount);assert.equal(JSON.stringify(a.parts),before);
 assert.equal(a.parts.reduce((n,p)=>n+p.indexCount,0),indexCount,'No additional indexed geometry or part draw resources');
 assert.deepEqual(api.bowelWallPeelSteps(a,whole(a)).map(s=>s.name),['Fibrous capsule','Parenchyma']);
 const steps=api.bowelWallPeelSteps(a,registry.lobules);assert.equal(steps.length,1);assert.equal(steps[0].name,'Parenchyma');assert.deepEqual(steps[0].ids,marker.lobuleParts);
 for(const mutate of [b=>b.thyroidWallAssembly.lobuleParts.push(b.thyroidWallAssembly.capsuleParts[0]),
  b=>b.thyroidWallAssembly.lobuleParts.push(b.thyroidWallAssembly.lobuleParts[0]),
  b=>b.parts.push({...b.parts[4],id:'overlapping-lobule',tissueWallLayer:'lobules'}),
  b=>b.parts[4].tissueWallLayer='lobules',
  b=>b.thyroidWallAssembly.parenchymaParts.push(b.thyroidWallAssembly.parenchymaParts[0])]){
  const bad=structuredClone(a);mutate(bad);assert.equal(api.thyroidWallRegistry(bad),null);
 }
});

test('Illustrative parenchymal aliases peel existing cells without claiming a microscopic TA98 identity',()=>{
 const a=fixture(),marker=a.thyroidWallAssembly,before=JSON.stringify(a.parts);
 marker.illustrativeDomainParts=marker.parenchymaParts.slice(0,2);marker.illustrativeDomainConcept='local:parenchymal-territories';
 a.concepts.push({id:marker.illustrativeDomainConcept,name:'Illustrative parenchymal domains',
  elements:[...marker.illustrativeDomainParts],partOf:marker.parenchymaConcept});
 const registry=api.thyroidWallRegistry(a);assert.ok(registry);assert.equal(registry.lobules,undefined);
 assert.equal(registry.domains.ta98Term,undefined);assert.equal(JSON.stringify(a.parts),before);
 assert.deepEqual(api.bowelWallPeelSteps(a,registry.domains).map(s=>s.ids),[marker.illustrativeDomainParts]);
 for(const mutate of [b=>b.concepts.at(-1).ta98Term='A11.3.00.009',
  b=>b.concepts.at(-1).partOf=marker.wholeConcept,
  b=>b.thyroidWallAssembly.illustrativeDomainParts.push(marker.capsuleParts[0])]){
  const bad=structuredClone(a);mutate(bad);assert.equal(api.thyroidWallRegistry(bad),null);
 }
});

test('Distinct stromal geometry joins internal tissue peeling and must have exactly one anatomical region owner',()=>{
 const a=fixture(),marker=a.thyroidWallAssembly,id='THYROID:male:right-lobe:stroma';
 const stroma={...a.parts[4],id,conceptId:id,name:'Stroma',indices:3000,normalSigns:3100,indexCount:24,tissueWallLayer:'stroma',
  tissueDrawGroups:[{start:0,count:12,neighbor:a.parts[0].id},{start:12,count:12,neighbor:a.parts[4].id}]};
 for(const p of [a.parts[0],a.parts[4]]){p.tissueDrawGroups.push({start:p.indexCount,count:12,neighbor:id});p.indexCount+=12;}
 a.parts.push(stroma);marker.stromaParts=[id];marker.stromaConcept='THYROID:male:stroma';
 a.concepts.push({id:marker.stromaConcept,name:'Stroma',elements:[id],ta98Term:'A11.3.00.007',ta98Kind:'exact',partOf:marker.wholeConcept});
 whole(a).elements.push(id);a.concepts[1].elements.push(id);
 const registry=api.thyroidWallRegistry(a);assert.ok(registry);
 const steps=api.bowelWallPeelSteps(a,whole(a));assert.equal(steps.length,2);assert.equal(steps[1].name,'Internal thyroid tissue');assert.ok(steps[1].ids.includes(id));
 assert.equal(api.bowelWallPeelSteps(a,registry.regions[0])[1].ids.length,2);
 const unowned=structuredClone(a);unowned.concepts[1].elements=unowned.concepts[1].elements.filter(p=>p!==id);assert.equal(api.thyroidWallRegistry(unowned),null);
 const twice=structuredClone(a);twice.concepts[2].elements.push(id);assert.equal(api.thyroidWallRegistry(twice),null);
 const overlap=structuredClone(a);overlap.thyroidWallAssembly.stromaParts=[marker.parenchymaParts[0]];assert.equal(api.thyroidWallRegistry(overlap),null);
});

test('Explicit coupled vessel partners preserve endocrine-only peeling and reject unregistered or disconnected walls',()=>{
 const a=fixture(),marker=a.thyroidWallAssembly,id='THYROIDVESSEL:male:arterial-network';
 const vessel={...a.parts[4],id,conceptId:id,name:'Thyroid arterial network',system:'arterial',
  tissueWallAssembly:'thyroid-vascular:male',tissueWallLayer:'vascular-wall',indexCount:24,
  tissueDrawGroups:[{start:0,count:12},{start:12,count:12,neighbor:a.parts[4].id}]};
 a.parts[4].tissueDrawGroups.push({start:a.parts[4].indexCount,count:12,neighbor:id});a.parts[4].indexCount+=12;
 a.parts.push(vessel);marker.vascularParts=[id];
 const r=api.thyroidWallRegistry(a);assert.ok(r);assert.equal(r.members.length,8);
 assert.deepEqual(api.bowelWallPeelSteps(a,whole(a)).map(s=>s.ids.length),[4,4]);
 assert.ok(!api.bowelWallPeelSteps(a,whole(a)).flatMap(s=>s.ids).includes(id),'Blood vessels never become an endocrine peel coat');
 for(const mutate of [b=>delete b.thyroidWallAssembly.vascularParts,
  b=>b.thyroidWallAssembly.vascularParts.push(b.thyroidWallAssembly.capsuleParts[0]),
  b=>b.parts.at(-1).system='nervous',b=>b.parts.at(-1).sectionSolid=false,
  b=>b.parts.at(-1).sourceOffset=[0,0,0],b=>b.parts.at(-1).tissueWallAssembly='thyroid-vascular:female',
  b=>b.parts.at(-1).positions+=4,b=>b.parts.at(-1).tissueDrawGroups[1].count=9,
  b=>{b.parts[4].tissueDrawGroups.pop();b.parts[4].indexCount-=12;b.parts.at(-1).tissueDrawGroups=[{start:0,count:24}];}]){
  const bad=structuredClone(a);mutate(bad);assert.equal(api.thyroidWallRegistry(bad),null);
 }
});
