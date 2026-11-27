import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';

const modules=['vascular-wall','bowel-wall-peeling','anatomical-components','viewer-state','anatomy-label-plan','anatomy-hierarchy','guest-hierarchy','guest-choice','ta98-modeled-parent'];
const bundle=await build({stdin:{contents:modules.map(name=>`export * from './app/${name}';`).join('\n'),resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
// Actual candidate metadata and packed offsets; this test does not certify anatomy.
const candidate=JSON.parse(await readFile(new URL('./fixtures/thoracic-aorta-controls.json',import.meta.url),'utf8'));
const sourceId='ZA:Thoracic aorta',prefix='TA98WALL:male:thoracic-aorta',keys=['externa','media','intima'],ids=keys.map(k=>`${prefix}:${k}`),terms=['A12.0.00.017','A12.0.00.019','A12.0.00.018'];
function fixture(){const a=structuredClone(candidate);a.parts=a.parts.filter(p=>p.id===sourceId||ids.includes(p.id));const concepts=new Set([sourceId,...a.parts.map(p=>p.conceptId)]);a.concepts=a.concepts.filter(c=>concepts.has(c.id));return a;}
const whole=a=>a.concepts.find(c=>c.id===sourceId),byId=a=>new Map(a.parts.map(p=>[p.id,p]));

test('Actual aorta replacement validates its UV bank and rejects corrupted finite declarations',()=>{
 const a=fixture(),r=api.vascularWallIntegration(a);assert.ok(r);assert.deepEqual(r.ids,ids);
 assert.equal(r.configuration.wholeTerm,'A12.2.11.001');
 for(const change of [{uvs:undefined},{uvs:999},{positions:999},{normals:999},{vertexCount:999},{chunk:999},{normalSigns:undefined},{geometryOrigin:[1,0,0]},{pairedTissue:false},{sectionSolid:false},{tissueDrawGroups:[]},{provenance:{url:'https://example.com/donor'}},{sourceOffset:[0,0,0]}]){
  const bad=fixture();const i=bad.parts.findIndex(p=>p.id===ids[0]);bad.parts[i]={...bad.parts[i],...change};assert.equal(api.vascularWallIntegration(bad),null,JSON.stringify(change));
 }
 const bad=fixture();bad.concepts.find(c=>c.id===ids[0]+':concept').partOf='other';assert.equal(api.vascularWallIntegration(bad),null);
});

test('Actual aorta selection, component labels and peeling address all three finite coats',()=>{
 const a=fixture(),choice=whole(a),parts=byId(a);assert.deepEqual(api.resolveSelection(a,[sourceId]),ids);assert.deepEqual(api.anatomicalComponentIds(choice,[],parts,a.concepts),ids);
 const steps=api.bowelWallPeelSteps(a,choice);assert.deepEqual(steps.map(s=>s.ids[0]),ids);
 for(const id of ids){assert.equal(api.nativePartFrame(parts.get(id)),'ZA:');assert.deepEqual(api.storedAssemblyParent(a.concepts.find(c=>c.id===id+':concept'),a.concepts,parts).elements,ids);}
 for(let depth=0;depth<=3;depth++){const view={selected:ids,visible:['arterial'],hidden:ids.slice(0,depth),explode:0,isolate:true,view:'front',rotate:false,reset:0};for(let k=0;k<3;k++){assert.equal(api.partVisible(parts.get(ids[k]),view,a.parts),k>=depth);assert.equal(api.sectionPartVisible(parts.get(ids[k]),view),k>=depth);}}
 const labels=api.selectedAnatomyLabels(a,ids,choice,false);assert.equal(labels.length,1);assert.deepEqual(labels[0].components,ids);
 const expanded=api.selectedAnatomyLabels(a,ids,choice,true);assert.equal(expanded.length,3);for(const label of expanded){assert.equal(label.members.length,1);assert.deepEqual(api.anatomyLabelChoice(label).elements,label.members);}
});

test('Actual aorta has TA98 and systems controls without changing generic coat authority',async()=>{
 const a=fixture(),authority=JSON.parse(await readFile('public/assets/ta98.json','utf8'));
 const all=new Map(),visit=n=>{all.set(n.id,n);n.children.forEach(visit);};api.resolveGuestHierarchy(a,authority).forEach(visit);
 const wall=all.get(prefix+':wall-assembly');assert.ok(wall);assert.deepEqual(wall.parts.map(p=>p.id),ids);assert.deepEqual(wall.children.map(n=>n.parts.map(p=>p.id)),ids.map(id=>[id]));
 for(let k=0;k<3;k++)assert.deepEqual(all.get(terms[k]).directParts.map(p=>p.id),[ids[k]]);
 const entries=api.hierarchyEntries(a).filter(e=>e.parts.length===1&&ids.includes(e.parts[0].id));assert.equal(entries.length,3);for(const e of entries){assert.equal(e.taParentConcept,sourceId);assert.deepEqual(e.ancestry,['Thoracic aorta']);}
});

test('Two explicit vascular assemblies coexist and choices never use the first unrelated registry',()=>{
 const a=fixture(),portal='ZA:Hepatic portal vein',pids=keys.map(k=>`VEINWALL:male:portal:${k}`),source=a.parts.find(p=>p.id===sourceId);
 a.parts.push({...source,id:portal,conceptId:portal,sourceId:'Hepatic portal vein',name:'Hepatic portal vein',system:'venous',replacementParts:pids});
 a.concepts.push({id:portal,name:'Hepatic portal vein',elements:pids,ta98Term:'A12.3.12.001',ta98Kind:'exact'});
 for(let k=0;k<3;k++){const p=structuredClone(a.parts.find(p=>p.id===ids[k]));Object.assign(p,{id:pids[k],conceptId:pids[k],sourceId:'Hepatic portal vein',system:'venous',tissueChart:12345,tissueWallAssembly:'portal:male:native-trunk'});p.tissueDrawGroups=p.tissueDrawGroups.map(g=>({...g,...(g.neighbor?{neighbor:pids[ids.indexOf(g.neighbor)]}:{})}));a.parts.push(p);a.concepts.push({id:p.id,name:p.name,elements:[p.id],partOf:portal,ta98Term:terms[k],ta98Kind:'exact'});}
 a.vascularWallIntegration.push({sex:'male',sourcePartId:portal,wholeConceptId:portal,assemblyId:'portal:male:native-trunk',partIds:pids});
 assert.equal(api.vascularWallIntegrations(a).length,2);assert.deepEqual(api.bowelWallPeelSteps(a,whole(a)).map(s=>s.ids[0]),ids);assert.deepEqual(api.vascularWallComponentIds(whole(a),byId(a)),ids);
 const pchoice=a.concepts.find(c=>c.id===portal);assert.deepEqual(api.bowelWallPeelSteps(a,pchoice).map(s=>s.ids[0]),pids);
});
