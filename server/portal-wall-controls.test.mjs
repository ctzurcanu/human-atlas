import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';

const modules=['vascular-wall','bowel-wall-peeling','anatomical-components','anatomical-relations','viewer-state','anatomy-label-plan','anatomy-hierarchy','guest-hierarchy','guest-choice','ta98-modeled-parent','local-ta98-layer-catalogue','view-url'];
const bundle=await build({stdin:{contents:modules.map(name=>`export * from './app/${name}';`).join('\n'),resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const keys=['externa','media','intima'],terms=['A12.0.00.017','A12.0.00.019','A12.0.00.018'];
const ids=keys.map(key=>`VEINWALL:male:portal:${key}`),sourceId='ZA:Hepatic portal vein',assembly='portal:male:native-trunk';
// Metadata-only control fixture: it makes no claim about anatomical geometry.
function fixture(){
 const common={system:'venous',chunk:7,positions:0,normals:160,tissueChart:320,vertexCount:12,indexCount:36,bounds:[[-.02,1.2,.02],[0,1.3,.04]],sourceId:'Hepatic portal vein',license:'CC BY-SA 4.0',provenance:{label:'Native source',url:'https://github.com/Z-Anatomy/Models-of-human-anatomy',detail:'Control fixture'},groups:['Systemic veins','Trunk']};
 const source={...common,id:sourceId,conceptId:sourceId,name:'Hepatic portal vein',indices:480,suppressed:true,replacementParts:[...ids]};
 const coats=keys.map((key,k)=>({...common,id:ids[k],conceptId:ids[k],name:`Tunica ${key} of hepatic portal vein`,pairedTissue:true,sectionSolid:true,tissueWallAssembly:assembly,tissueWallLayer:key,normalSigns:640+k*12,indices:800+k*144,geometryOrigin:[0,1.25,0],tissueDrawGroups:[{start:0,count:12,...(k>0?{neighbor:ids[k-1]}:{})},{start:12,count:12,...(k<2?{neighbor:ids[k+1]}:{})},{start:24,count:12}]}));
 const liver={...common,id:'ZA:Liver:Organ',conceptId:'ZA:Liver',sourceId:'Liver',name:'Liver · Organ',system:'digestive',indices:1000};
 const smv={...common,id:'ZA:Superior mesenteric vein',conceptId:'ZA:Superior mesenteric vein',sourceId:'Superior mesenteric vein',name:'Superior mesenteric vein',indices:1200};
 const splenic={...common,id:'ZA:Splenic vein',conceptId:'ZA:Splenic vein',sourceId:'Splenic vein',name:'Splenic vein',indices:1400};
 return {version:'control-test',sex:'male',parts:[source,...coats,liver,smv,splenic],concepts:[{id:sourceId,name:source.name,elements:[...ids],ta98Term:'A12.3.12.001',ta98Kind:'exact'},...coats.map((p,k)=>({id:p.id,name:p.name,elements:[p.id],partOf:sourceId,ta98Term:terms[k],ta98Kind:'exact'})),{id:liver.conceptId,name:'Liver',elements:[liver.id],ta98Term:'A05.8.01.001',ta98Kind:'exact'},{id:smv.id,name:smv.name,elements:[smv.id],ta98Term:'A12.3.12.005',ta98Kind:'exact'},{id:splenic.id,name:splenic.name,elements:[splenic.id],ta98Term:'A12.3.12.002',ta98Kind:'exact'}],vascularWallIntegration:[{sex:'male',sourcePartId:sourceId,wholeConceptId:sourceId,assemblyId:assembly,partIds:[...ids]}],chunks:[],triangles:0};
}
const whole=a=>a.concepts.find(c=>c.id===sourceId);
const byId=a=>new Map(a.parts.map(p=>[p.id,p]));
const state={selected:ids,visible:['venous'],hidden:[],explode:0,isolate:true,view:'front',rotate:false,reset:0};

test('Only the explicit native shared-bank finite coats activate the registry',()=>{
 const a=fixture();assert.ok(api.vascularWallIntegration(a));
 for(const change of [{sourceOffset:[0,0,0]},{suppressed:true},{pairedTissue:false},{sectionSolid:false},{chunk:999},{positions:999},{normals:999},{vertexCount:999},{geometryOrigin:[0,0,0]},{tissueChart:undefined},{normalSigns:undefined},{tissueDrawGroups:[]},{sourceId:'Donor portal vein'},{license:'unverified'},{provenance:{url:'https://example.com/donor'}},{tissueWallLayer:'media'}]){
  const bad=fixture();bad.parts[1]={...bad.parts[1],...change};assert.equal(api.vascularWallRegistry(bad.parts),null,JSON.stringify(change));
 }
 for(const change of [{suppressed:false},{sourceOffset:[0,0,0]},{replacementParts:ids.slice(0,2)}]){
  const bad=fixture();bad.parts[0]={...bad.parts[0],...change};assert.equal(api.vascularWallRegistry(bad.parts),null);
 }
 const unmarked=fixture();delete unmarked.vascularWallIntegration;assert.equal(api.vascularWallIntegration(unmarked),null);assert.equal(api.bowelWallPeelSteps(unmarked,whole(unmarked)),null);
 const opposite=fixture();opposite.vascularWallIntegration[0].sex='female';assert.equal(api.vascularWallIntegration(opposite),null);
 for(const change of [{partOf:'unrelated'},{ta98Term:'A12.0.00.021'},{ta98Kind:'parent'},{elements:[]}]){const bad=fixture();bad.concepts[1]={...bad.concepts[1],...change};assert.equal(api.vascularWallIntegration(bad),null);}
 const wrongSex=fixture();wrongSex.sex='female';assert.equal(api.vascularWallIntegration(wrongSex),null);
 const badWhole=fixture();badWhole.concepts[0].elements=ids.slice(0,2);assert.equal(api.vascularWallIntegration(badWhole),null);
});

test('Whole selection expands to the three coats, excluding tributaries from components',()=>{
 const a=fixture(),parts=byId(a),choice=whole(a);
 assert.deepEqual(api.resolveSelection(a,[sourceId]),ids);
 assert.deepEqual(api.resolveSelection(a,['Hepatic portal vein']),ids);
 assert.deepEqual(api.anatomicalComponentIds(choice,[],parts,a.concepts),ids);
 const misleading=[{kind:'contains',target:parts.get('ZA:Splenic vein')},{kind:'contains',target:parts.get('ZA:Superior mesenteric vein')}];
 assert.deepEqual(api.anatomicalComponentRelations(choice,misleading,parts,a.concepts).map(r=>r.target.id),ids);
 for(const id of ids){const child=a.concepts.find(c=>c.id===id);assert.deepEqual(api.storedAssemblyParent(child,a.concepts,parts).elements,ids);assert.equal(api.nativePartFrame(parts.get(id)),'ZA:');}
 assert.equal(api.vascularWallComponentIds({...choice,id:'selection-set'},parts),undefined);
 assert.equal(api.vascularWallComponentIds({...choice,elements:ids.slice(0,2)},parts),undefined);
});

test('Peeling selected coats hides them in both main and section views, then restores them',()=>{
 const a=fixture(),steps=api.bowelWallPeelSteps(a,whole(a));assert.deepEqual(steps.map(s=>s.ids[0]),ids);
 for(let depth=0;depth<=3;depth++){
  const view={...state,hidden:ids.slice(0,depth)};assert.equal(api.wallPeelDepth(steps,view.hidden),depth);
  for(let k=0;k<3;k++){const p=byId(a).get(ids[k]);assert.equal(api.partVisible(p,view,a.parts),k>=depth);assert.equal(api.sectionPartVisible(p,view),k>=depth);}
 }
 assert.equal(api.bowelWallPeelSteps(a,{...whole(a),id:'selection-set'}),null);
 assert.equal(api.bowelWallPeelSteps(a,{...whole(a),elements:ids.slice(0,2)}),null);
 assert.equal(api.partVisible(a.parts[0],state,a.parts),false,'Retired original exterior stays hidden');
});

test('Standard labels group the whole and target separate positions on expanded coats',()=>{
 const a=fixture(),choice=whole(a),grouped=api.selectedAnatomyLabels(a,ids,choice,false);
 assert.equal(grouped.length,1);assert.equal(grouped[0].id,sourceId);assert.deepEqual(grouped[0].components,ids);
 const expanded=api.selectedAnatomyLabels(a,ids,choice,true);assert.equal(expanded.length,3);
 const targets=api.layerLabelTargets(expanded,a.parts,a.concepts);assert.equal(targets.size,3);
 assert.equal(new Set([...targets.values()].map(v=>v.coordinate)).size,3);
 for(const item of expanded){assert.deepEqual(api.anatomyLabelChoice(item).elements,item.members);assert.equal(item.members.length,1);const t=targets.get(item.id);assert.equal(t.axis,1);assert.ok(t.coordinate>=1.225&&t.coordinate<=1.275);}
});

test('Portal coats have contextual hierarchy controls while generic TA98 parents are retained',async()=>{
 const a=fixture(),authority=JSON.parse(await readFile('public/assets/ta98.json','utf8')),parents=JSON.parse(await readFile('app/data/ta98-parents.json','utf8'));
 const nodes=api.resolveGuestHierarchy(a,authority),all=new Map();const visit=n=>{all.set(n.id,n);n.children.forEach(visit);};nodes.forEach(visit);
 const wall=all.get('VEINWALL:male:portal:wall-assembly');assert.ok(wall);assert.equal(wall.children.length,3);assert.deepEqual(wall.parts.map(p=>p.id),ids);
 assert.deepEqual(wall.children.map(n=>n.parts.map(p=>p.id)),ids.map(id=>[id]));
 assert.deepEqual(api.anatomicalComponentIds(api.guestNodeChoice('ta98',wall),[],byId(a),a.concepts),ids);
 for(let k=0;k<3;k++){assert.deepEqual(all.get(terms[k]).directParts.map(p=>p.id),[ids[k]]);assert.equal(parents.parents[terms[k]],'A12.0.00.000');}
 const entries=api.hierarchyEntries(a).filter(e=>ids.includes(e.id));assert.equal(entries.length,3);
 for(const entry of entries){assert.equal(entry.taParentConcept,sourceId);assert.deepEqual(entry.ancestry,['Hepatic portal vein']);assert.equal(api.systemPathFor(entry).at(-1),'Hepatic portal vein');}
});

test('Physiological flow belongs to the whole replacement, with reciprocal whole selection',()=>{
 const a=fixture(),parts=byId(a),members=ids.map(id=>parts.get(id));
 const links=api.anatomicalConceptFunctionalRelations(whole(a),members,a.parts);
 assert.ok(links.some(r=>r.kind==='after'&&r.target.id==='ZA:Liver:Organ'));
 for(const row of links)assert.ok(!ids.includes(row.target.id));
 const liverLinks=api.anatomicalRelations(parts.get('ZA:Liver:Organ'),a.parts);
 const portal=liverLinks.find(r=>r.kind==='before'&&r.targetConcept?.id===sourceId);assert.ok(portal);assert.deepEqual(portal.targetConcept.elements,ids);
 for(const p of members){const rows=api.anatomicalRelations(p,a.parts);assert.ok(!rows.some(r=>['before','after','arterial','venous','innervation'].includes(r.kind)),'No whole-vessel flow or invented wall supply on a coat');}
 assert.deepEqual(api.anatomicalConceptFunctionalRelations({...whole(a),id:'selection-set'},members,a.parts),[]);
 assert.equal(api.vascularWallWholeRelations(whole(a),members.slice(0,2),a.parts,api.anatomicalRelations),null);
 const donor=fixture();donor.parts.find(p=>p.id==='ZA:Liver:Organ').id='HRA:donor-liver';assert.ok(!api.anatomicalConceptFunctionalRelations(whole(donor),ids.map(id=>byId(donor).get(id)),donor.parts).some(r=>r.target.id==='HRA:donor-liver'));
});

test('Portal-wall review routing is opt-in and development-only, including compact copied views',()=>{
 const q=new URLSearchParams('venous=portal-wall');
 assert.equal(api.localTa98LayerCatalogue('local-ta98',q,true),'/local-models/ta98-review/male-portal-wall/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-ta98',q,false),undefined);
 assert.ok(!api.localTa98LayerCatalogue('local-female-ta98',q,true).includes('/male-portal-wall/'));
 for(const [model,dev] of [['local-female-ta98',true],['local-ta98',false],['male',true]]){const u=new URL('http://localhost:3016/?model='+model);api.preserveLocalTa98LayerCatalogueQuery(u,q,model,dev);assert.equal(u.searchParams.has('venous'),false);}
 const url=new URL('http://localhost:3016/?model=local-ta98&state='+Array.from({length:600},(_,i)=>i.toString(36)).join('.'));api.preserveLocalTa98LayerCatalogueQuery(url,q,'local-ta98',true);
 const copied=api.compactViewUrl(url.href);assert.ok(new URL(copied).hash);assert.equal(new URL(copied).searchParams.get('venous'),'portal-wall');assert.equal(api.viewParameters(copied).get('venous'),'portal-wall');assert.equal(api.localTa98LayerCatalogue('local-ta98',new URL(copied).searchParams,true),api.localTa98LayerCatalogue('local-ta98',q,true));
});
