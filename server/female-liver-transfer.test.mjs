import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:`export * from './app/hierarchy-navigation';export * from './app/preferred-source-representations';export * from './app/reviewed-catalogue-identities';export * from './app/liver-controls';export * from './app/liver-completion-controls';export * from './app/liver-surface-features';export * from './app/hepatic-connection-controls';export * from './app/liver-section-volumes';export * from './app/bowel-wall-peeling';export * from './app/viewer-state';`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const path='.local-models/geometry-workbench/liver-completion/female-transfer-v206/atlas.json';
const raw=fs.existsSync(path)?JSON.parse(fs.readFileSync(path)):null,options={skip:!raw&&'Local geometry is not redistributed.'};
const walk=n=>[n,...(n.children??[]).flatMap(walk)];
test('female transfer exposes one complete Systems liver with the same 63/69 controls',options,()=>{
 const a=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(raw)));
 assert.ok(api.liverControlRegistry(a));assert.ok(api.hepaticConnectionRegistry(a));assert.ok(api.liverSectionRegistry(a));
 const m=api.liverCompletionRegistry(a);assert.ok(m);assert.equal(m.partIds.length,16);
 const roots=walk(api.hierarchyNavigation(a,'systems')).filter(n=>n.name==='Liver'&&n.terminology?.ta98==='A05.8.01.001');assert.equal(roots.length,1);
 const terms=walk(roots[0]).filter(n=>n.terminology?.ta98?.startsWith('A05.8.01.'));assert.equal(new Set(terms.map(n=>n.terminology.ta98)).size,69);assert.equal(terms.filter(n=>n.elements.length).length,63);
 for(const id of [...m.functionalPartIds,...m.microscopicPartIds])assert.ok(!roots[0].elements.includes(id));
 for(const id of raw.liverTransferIntegration.retiredPartIds)assert.equal(a.parts.find(p=>p.id===id).suppressed,true);
});
test('female direct links, TA98, peeling and functional visibility use the transferred assembly',options,()=>{
 const authority=JSON.parse(fs.readFileSync('public/assets/ta98.json'));
 const tree=walk(api.hierarchyNavigation(raw,'guest:ta98',authority));
 for(const marker of ['liverIntegration','liverFeatureIntegration','liverCompletionIntegration'])for(const g of raw[marker].termGroups){const n=tree.find(n=>n.id==='hierarchy:guest:ta98:'+g.term);assert.ok(n,g.term);for(const id of g.partIds??g.parts)assert.ok(n.elements.includes(id),g.term+' '+id);}
 const whole=raw.concepts.find(c=>c.id==='LOCAL:female:liver');assert.equal(api.bowelWallPeelSteps(raw,whole).length,2);
 const m=api.liverCompletionRegistry(raw),p=raw.parts.find(p=>p.id==='TA98PART:ZA:Liver:right');
 const state={selected:[m.functionalPartIds[0]],hidden:[],visible:['digestive'],depthHidden:[],isolate:false,region:'all',skinOpacity:0};
 assert.equal(api.partVisible(p,state,raw.parts),false);assert.equal(api.sectionPartVisible(p,state),false);
 for(const g of m.termGroups.filter(g=>g.detailModel))assert.match(api.liverMicroscopicDetailUrl(raw.concepts.find(c=>c.id===g.conceptId)),/model=liver-unit/);
});
test('female live section geometry applies its lossless source scale exactly once',options,()=>{
 const p=raw.parts.find(p=>p.id==='TA98PART:ZA:Liver:right'),row=raw.liverSectionIntegration.volumes.find(v=>v.partId===p.id);
 const data=fs.readFileSync(raw.chunks[p.chunk].url.replace(/^\/local-models\//,'.local-models/')),buffer=Uint8Array.from(data).buffer,cache=new Map();
 const shared=(key,create)=>{if(!cache.has(key))cache.set(key,create());return cache.get(key);};
 const v=api.liverSectionGeometry(raw,p,buffer,0,shared);assert.ok(v);const points=v.geometry.getAttribute('position');
 for(let i=0;i<Math.min(row.vertexCount,50);i++)for(let k=0;k<3;k++)assert.equal(points.array[i*3+k],Math.fround(data.readFloatLE(row.positions+(i*3+k)*4)*p.geometryScale));
 for(const g of [v.geometry,...v.exclusions])g.dispose();assert.equal(raw.liverTransferIntegration.newGeometryBytes,0);
});
