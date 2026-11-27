import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/ta98-term-bindings';export * from './app/guest-hierarchy';export * from './app/local-ta98-layer-catalogue';export * from './app/anatomy-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const term='A06.2.03.010',pid='CRICOTRACHEAL:female:ligament';
const fixture=()=>({parts:[{id:pid,conceptId:pid,name:'Cricotracheal ligament',system:'respiratory',sectionSolid:true,ta98Term:term,coverageLimitation:'Attachment remains under review.',bounds:[[0,0,0],[1,1,1]]}],concepts:[{id:pid,name:'Cricotracheal ligament',elements:[pid],ta98Term:term,ta98Kind:'exact'}],ta98TermBindings:[{term,conceptId:pid,partIds:[pid]}]});
test('a declared finite tissue gets a canonical TA98 control without a stale static selector',()=>{
 const a=fixture(),hierarchy={schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98',roots:[term],nodes:[{id:term,name:'cricotracheal ligament',matches:[],description:'TA98: '+term+'\nNo direct reviewed mesh correspondence.'}]};
 const [node]=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(hierarchy));
 assert.deepEqual(node.directParts.map(p=>p.id),[pid]);assert.deepEqual(node.parts.map(p=>p.id),[pid]);
 assert.match(node.description,/Attachment remains under review/);assert.doesNotMatch(node.description,/No direct reviewed mesh|duct/);
});
test('absent, hidden, suppressed, mismatched or nonfinite records cannot acquire a binding',()=>{
 const mutations=[a=>a.parts.pop(),a=>a.parts[0].suppressed=true,a=>a.parts[0].defaultHidden=true,a=>a.defaultHidden=[pid],a=>a.parts[0].sectionSolid=false,a=>a.parts[0].ta98Term='A06.2.03.009',a=>a.parts[0].conceptId='unrelated',a=>a.concepts[0].ta98Kind='parent',a=>a.concepts[0].elements.push('unrelated'),a=>a.ta98TermBindings[0].partIds.push(pid),a=>a.ta98TermBindings[0]=null];
 for(const mutate of mutations){const a=fixture();mutate(a);assert.equal(api.catalogueTa98TermParts(a).size,0);}
});
test('catalogues without explicit bindings retain their existing hierarchy behavior',()=>{
 const a=fixture();delete a.ta98TermBindings;assert.equal(api.catalogueTa98TermParts(a).size,0);
});
test('an explicit source identity supplies an ordinary TA98 control without claiming finite tissue',()=>{
 const a=fixture();a.parts[0].sectionSolid=false;
 a.ta98TermBindings[0].representation='source-identity';
 a.concepts[0].coverageLimitation='Named source vessel; target body fit remains under review.';
 a.parts[0].coverageLimitation=a.concepts[0].coverageLimitation;
 const root='A12.2.00.001';
 const hierarchy={schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98',roots:[root],nodes:[
  {id:root,name:'arteries',childIds:[term]},
  {id:term,name:'source vessel',matches:[],description:'TA98: '+term+'\nNo direct reviewed mesh correspondence.'},
 ]};
 const [parent]=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(hierarchy));
 assert.deepEqual(parent.parts.map(p=>p.id),[pid]);
 assert.deepEqual(parent.children[0].directParts.map(p=>p.id),[pid]);
 assert.match(parent.children[0].description,/identity is established separately from body fit/i);
 assert.doesNotMatch(parent.children[0].description,/Stored finite tissue|No direct reviewed mesh/);
});
test('source identity controls require an exact agreement and an explicit qualification',()=>{
 const mutations=[a=>delete a.concepts[0].coverageLimitation,a=>a.concepts[0].ta98Kind='parent',a=>a.parts[0].ta98Term='A06.2.03.009',a=>a.parts[0].suppressed=true,a=>a.parts[0].defaultHidden=true,a=>a.defaultHidden=[pid],a=>a.concepts[0].elements.push('unrelated')];
 for(const mutate of mutations){
  const a=fixture();a.parts[0].sectionSolid=false;a.ta98TermBindings[0].representation='source-identity';a.concepts[0].coverageLimitation='Named source; geometry fit is separate.';
  mutate(a);assert.equal(api.catalogueTa98TermParts(a).size,0);
 }
});
test('blood vessels retain Circulatory membership across organ and sensory TA98 chapters',()=>{
 for(const system of ['arterial','venous'])for(const ta98Term of ['A08.1.03.003','A15.2.05.008']){
  const entry={name:'Named vessel',system,ta98Term,parts:[{id:'vessel',conceptId:'vessel',name:'Named vessel',system}]};
  assert.equal(api.majorSystemFor(entry),'circulatory');
  assert.equal(api.majorSystemFor({...entry,assemblyMajorSystem:'urinary'}),'urinary','A reviewed organ-assembly override remains authoritative');
 }
});
test('an explicitly qualified partial occurrence retains its parent-term classification',()=>{
 const a=fixture();a.concepts[0].ta98Kind='parent';a.ta98TermBindings[0].kind='parent';
 assert.equal(api.catalogueTa98TermParts(a).size,0);
 a.concepts[0].coverageLimitation='Only the first occurrence is represented; the remaining collection is missing.';
 assert.deepEqual(api.catalogueTa98TermParts(a).get(term).map(p=>p.id),[pid]);
});
test('ligament and partial-cartilage reviews remain local and survive share URLs',()=>{
 for(const [choice,path]of [['continuous-rim-v102','female-cricotracheal-continuous-rim-v102/atlas-controls.json'],['first-ring-v122','cricotracheal-first-cartilage-v122/atlas.json']]){
  const q=new URLSearchParams({cricotracheal:choice}),url=new URL('http://localhost:3016/?model=local-female-ta98');
  assert.ok(api.localTa98LayerCatalogue('local-female-ta98',q,true).endsWith(path));
  for(const model of ['local-female-ta98','local-ta98','female'])assert.equal(api.localTa98LayerCatalogue(model,q,false),undefined);
  assert.doesNotMatch(api.localTa98LayerCatalogue('local-ta98',q,true),/cricotracheal/);
  api.preserveLocalTa98LayerCatalogueQuery(url,q,'local-female-ta98',true);assert.equal(url.searchParams.get('cricotracheal'),choice);
 }
});
