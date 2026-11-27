import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/renal-wall';export * from './app/bowel-wall-peeling';export * from './app/viewer-state';export * from './app/guest-hierarchy';export * from './app/anatomical-relations';export * from './app/anatomical-components';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const read=(sex='male')=>JSON.parse(fs.readFileSync(`${process.env.ATLAS_RENAL_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas'}/${process.env.ATLAS_RENAL_CANDIDATE?sex+'.json':sex+'/atlas.json'}`));
const base={selected:[],hidden:[],visible:['urinary'],depthHidden:[],view:'front',region:'all',contextOpacity:0,skinOpacity:0,explode:0,isolate:false,rotate:false,reset:0,focus:0};
for(const sex of ['male','female']){
test(`${sex}: detailed kidney controls peel actual tissue groups`,()=>{
 const a=read(sex),registry=api.renalWallRegistry(a);assert.ok(registry);const layered=registry.marker.schema==='renal-local-layers-v1';assert.equal(registry.marker.partIds.length,layered?66:54);
 for(const g of registry.assemblies)for(const source of g.definition.sources){
  const c=a.concepts.find(c=>c.id===source.wholeConceptId),steps=api.bowelWallPeelSteps(a,c);
  const components=api.anatomicalComponentIds(c,[],new Map(a.parts.map(p=>[p.id,p])),a.concepts);
  assert.deepEqual(new Set(components),new Set(source.partIds.length>1?source.partIds:[]),'Components must contain only this stored physical assembly');
  if(source.partIds.length!==24){
   if(!layered){assert.equal(steps,null,'The unpartitioned collecting wall must not impersonate three histological coats');continue;}
   assert.ok(steps);assert.deepEqual(steps.map(s=>s.ids.length),source.partIds.length===6?[2,2,2]:[1,1,1]);assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(c.elements));
   assert.equal(api.renalWholeSource([a.parts.find(p=>p.id===c.elements[0])],a.parts),undefined,'One coat is not the complete organ');continue;
  }
  assert.ok(steps,c.id);assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(c.elements));
  assert.deepEqual(steps.map(s=>s.ids.length),[1,1,22]);
  assert.equal(api.bowelWallPeelSteps(a,{...c,id:'selection-set'}),null);
  assert.equal(api.bowelWallPeelSteps(a,{...c,elements:c.elements.slice(1)}),null);
 }
});
test(`${sex}: retired sources stay hidden and new tissues remain selectable`,()=>{
 const a=read(sex),r=api.renalWallRegistry(a);
 for(const g of r.assemblies)for(const source of g.definition.sources){
  const c=a.concepts.find(c=>c.id===source.wholeConceptId),state=api.readViewUrl('?select='+encodeURIComponent(c.id),a,base);
  assert.deepEqual(new Set(state.selected),new Set(c.elements));assert.ok(c.elements.every(id=>!state.hidden.includes(id)));
  if(source.sourcePartId){const hidden=api.readViewUrl('?hide='+encodeURIComponent(source.sourcePartId),a,base);assert.ok(c.elements.every(id=>hidden.hidden.includes(id)));}
 }
});
test(`${sex}: incompatible frames and incomplete surface ownership reject the registry`,()=>{
 for(const change of [a=>a.parts.find(p=>p.id===a.finiteRenalIntegration.partIds[0]).geometryOrigin=[0,0,0],a=>a.parts.find(p=>p.id===a.finiteRenalIntegration.partIds[0]).tissueDrawGroups[0].count-=3]){
  const a=read(sex);change(a);assert.equal(api.renalWallRegistry(a),null);
 }
});
test(`${sex}: actual cortex, medulla, pyramids and papillae resolve in TA98`,()=>{
 const a=read(sex),roots=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')))),nodes=new Map(),todo=[...roots];
 while(todo.length){const n=todo.pop();if(nodes.has(n.id))continue;nodes.set(n.id,n);todo.push(...n.children);}
 for(const term of ['A08.1.01.013','A08.1.01.015','A08.1.01.020','A08.1.01.027','A08.1.01.029'])assert.ok(nodes.get(term)?.directParts.some(p=>a.finiteRenalIntegration.partIds.includes(p.id)),term);
 for(const term of ['A08.1.01.016','A08.1.01.019','A08.1.01.021','A08.1.01.026'])assert.ok(!(nodes.get(term)?.directParts??[]).some(p=>a.finiteRenalIntegration.partIds.includes(p.id)),`${term} must not acquire unrelated finite geometry`);
 for(const term of ['A08.1.05.009','A08.1.05.010','A08.1.05.011','A08.2.01.005','A08.2.01.006','A08.2.01.007'])assert.equal((nodes.get(term)?.directParts??[]).some(p=>a.finiteRenalIntegration.partIds.includes(p.id)),a.finiteRenalIntegration.schema==='renal-local-layers-v1',term);
});
test(`${sex}: whole kidneys retain native supply and ipsilateral suprarenal links`,()=>{
 const a=read(sex),r=api.renalWallRegistry(a);
 for(const g of r.assemblies){
  const source=g.definition.sources.find(s=>s.partIds.length===24),c=a.concepts.find(c=>c.id===source.wholeConceptId),selected=source.partIds.map(id=>a.parts.find(p=>p.id===id));
  const links=api.anatomicalConceptFunctionalRelations(c,selected,a.parts);assert.ok(links.some(l=>l.kind==='arterial'),c.id);assert.ok(links.some(l=>l.kind==='venous'),c.id);
  assert.equal(api.renalWholeSource(selected.slice(1),a.parts),undefined);
  const adjacent=selected.flatMap(p=>api.anatomicalRelations(p,a.parts)).filter(l=>l.kind==='adjacent');assert.ok(adjacent.length);
  const adrenal=adjacent.find(l=>l.targetConcept?.id.includes(g.definition.side[0]+'_adrenal')||l.targetConcept?.id===`ZA:Suprarenal gland.${g.definition.side[0]}`);assert.ok(adrenal,c.id);
 }
});
test(`${sex}: reciprocal collecting joins and whole urine routes remain available`,()=>{
 const a=read(sex),r=api.renalWallRegistry(a),graph=api.renalTissueRelations(a.parts);
 for(const g of r.assemblies){
  const layered=r.marker.schema==='renal-local-layers-v1',source=g.definition.sources.find(s=>s.partIds.length===(layered?6:2)),c=a.concepts.find(c=>c.id===source.wholeConceptId),selected=source.partIds.map(id=>a.parts.find(p=>p.id===id));
  const routes=api.anatomicalConceptFunctionalRelations(c,selected,a.parts);assert.ok(routes.some(l=>l.kind==='after'&&l.targetConcept?.elements.length===(layered?3:1)));assert.ok(routes.some(l=>l.kind==='before'&&l.targetConcept?.elements.length===24));
  for(const p of selected){const next=graph.get(p.id)?.find(l=>l.kind==='after');assert.ok(next,p.id);assert.equal(next.target.tissueWallLayer,p.tissueWallLayer);assert.ok(graph.get(next.target.id)?.some(l=>l.kind==='before'&&l.target.id===p.id));assert.ok(api.anatomicalRelations(p,a.parts).some(l=>l.kind==='after'&&l.target.id===next.target.id));}
 }
});
}
