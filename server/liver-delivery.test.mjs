import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:`
 export * from './app/hierarchy-navigation';
 export * from './app/liver-controls';
 export * from './app/preferred-source-representations';
 export * from './app/reviewed-catalogue-identities';
 export * from './app/ta98-modeled-parent';
 export * from './app/hepatic-relations';
 export * from './app/bowel-wall-peeling';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const path='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
const source=fs.existsSync(path)?JSON.parse(fs.readFileSync(path)):null;
const options={skip:!source?.liverIntegration&&'Local review geometry is not redistributed.'};
const runtime=a=>api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(a)));
const walk=n=>[n,...(n.children??[]).flatMap(walk)];
const code=n=>n.ta98Term??n.terminology?.ta98;

test('Systems contains one complete liver assembly with its four lobes and capsule',options,()=>{
 const a=runtime(source),r=api.liverControlRegistry(a);assert.ok(r);
 const nodes=walk(api.hierarchyNavigation(a,'systems'));
 const liver=nodes.filter(n=>n.name==='Liver'&&code(n)==='A05.8.01.001');assert.equal(liver.length,1);
 for(const id of r.marker.partIds)assert.ok(liver[0].elements.includes(id),id);
 const children=walk(liver[0]);
 for(const term of ['A05.8.01.026','A05.8.01.027','A05.8.01.029','A05.8.01.030','A05.8.01.054'])
  assert.equal(children.filter(n=>code(n)===term).length,1,term);
 for(const id of source.liverIntegration.sourceParts){assert.equal(a.parts.find(p=>p.id===id)?.suppressed,true);assert.ok(!liver[0].elements.includes(id));}
});

test('TA98, direct links and Show components use the same stored liver without moving its meshes',options,()=>{
 const a=runtime(source),before=JSON.stringify(a),authority=JSON.parse(fs.readFileSync('public/assets/ta98.json'));
 const nodes=walk(api.hierarchyNavigation(a,'guest:ta98',authority));
 const liver=nodes.find(n=>code(n)==='A05.8.01.001');assert.ok(liver);
 for(const id of a.liverIntegration.partIds)assert.ok(liver.elements.includes(id),id);
 const direct=api.liverTermParts(a);
 const parts=new Map(a.parts.map(p=>[p.id,p]));const whole=a.concepts.find(c=>c.id==='ZA:Liver');
 const expectedComponents=[...a.liverIntegration.partIds];
 if(a.hepaticConnectionIntegration)expectedComponents.push('TA98HEPATIC:male:biliary','TA98HEPATIC:male:region:A05.8.01.061','TA98PART:FJ3079:A05.8.01.061');
 assert.deepEqual(new Set(api.ta98ModeledComponentIds(whole,a.concepts,parts)),new Set(expectedComponents));
 for(const g of a.liverIntegration.termGroups){
  assert.deepEqual(new Set(direct.get(g.term).map(p=>p.id)),new Set(g.parts));
  for(const id of g.parts)assert.ok(nodes.find(n=>code(n)===g.term)?.elements.includes(id),g.term);
  assert.deepEqual(new Set(a.concepts.find(c=>c.id===g.conceptId).elements),new Set(g.parts));
 }
 assert.equal(JSON.stringify(a),before,'hierarchy construction preserves every literal mesh and transform');
});

test('lobe blood-flow associations stay available and never certify an unbuilt lumen junction',options,()=>{
 const a=runtime(source),graph=api.hepaticRelations(a.parts);
 for(const p of a.parts.filter(api.isInstalledLiverLobe)){
  const links=graph.get(p.id)??[];assert.ok(links.some(r=>r.target.id===(a.hepaticConnectionIntegration?'TA98HEPATIC:male:portal':'ZA:Hepatic portal vein')));
  assert.ok(links.some(r=>r.target.id===(a.hepaticConnectionIntegration?'TA98HEPATIC:male:hepatic_vein':'ZA:Hepatic veins')));
  assert.match(api.hepaticSourceLimitation(p),a.liverCompletionIntegration?/Regional extents remain provisional/:/interfaces remain unfinished/);
 }
 const wrong=structuredClone(a);wrong.liverIntegration.partIds.pop();assert.equal(api.liverControlRegistry(wrong),null);
 const hidden=structuredClone(a);hidden.parts.find(p=>p.id===hidden.liverIntegration.partIds[0]).suppressed=true;assert.equal(api.liverControlRegistry(hidden),null);
});

test('liver peeling removes the capsule before all four parenchymal lobes together',options,()=>{
 const a=runtime(source),whole=a.concepts.find(c=>c.id==='ZA:Liver');
 const steps=api.bowelWallPeelSteps(a,whole);
 assert.equal(steps.length,2);
 assert.deepEqual(steps[0].ids,['TA98PART:ZA:Liver:capsule']);
 assert.deepEqual(new Set(steps[1].ids),new Set(a.parts.filter(api.isInstalledLiverLobe).map(p=>p.id)));
 for(const concept of a.concepts.filter(c=>c.id.startsWith('TA98PART:ZA:Liver:')))
  assert.equal(api.liverTissuePeelSteps(a,concept),null,'a selected lobe is not an invented wall stack');
 assert.equal(api.liverTissuePeelSteps(a,{...whole,elements:[...whole.elements,'ZA:Hepatic veins']}),null);
});
