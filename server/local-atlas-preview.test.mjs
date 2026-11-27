import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from './atlas-test-build.mjs';

const bundle=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/view-url';",resolveDir:process.cwd()},bundle:true,write:false,platform:'node',format:'esm'});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const source='/local-models/geometry-workbench/example/candidate.json';

test('an explicit candidate catalogue overrides only its local development review URL',()=>{
 const normal=new URLSearchParams('blood=biodigital-review');
 const preview=new URLSearchParams(normal);preview.set('atlasSource',source);
 const installed=api.localTa98LayerCatalogue('local-ta98',normal,true);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',preview,true),source);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',normal,true),installed);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',preview,false),undefined);
 assert.equal(api.localTa98LayerCatalogue('male-detail',preview,true),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',preview,true),source);
});

test('external, traversal, encoded and non-JSON catalogue URLs cannot become local preview sources',()=>{
 for(const value of ['https://example.com/a.json','//example.com/a.json','/local-models/../secret.json','/local-models/a/../../secret.json','/local-models/%2e%2e/secret.json','/local-models/a//b.json','/local-models/a/./b.json','/local-models/a.json?token=x','/local-models/a.bin','/local-models/a\\b.json']){
  const query=new URLSearchParams({atlasSource:value});
  assert.equal(api.localAtlasPreviewSource(query),undefined,value);
 }
});

test('compressed comparison links retain the exact same camera, selection and source',()=>{
 const frame='0.5,0.51,0.9999,0.5,0.5,0.5,0.5,0.44,0.27,0.5,1,0.5';
 const selection=Array.from({length:30},(_,i)=>'BDM:male:orbital-vein-'+i);
 const q=new URLSearchParams({model:'local-ta98',atlasSource:source,frame,isolate:'1',context:'0.08',skin:'0',labels:'0'});
 selection.forEach(id=>q.append('select',id));
 const url=api.compactViewUrl('http://127.0.0.1:3016/?'+q);
 assert.ok(new URL(url).hash.startsWith('#atlas-view=1.'));
 assert.equal(new URL(url).searchParams.get('atlasSource'),source);
 const restored=api.viewParameters(url);
 assert.equal(restored.get('frame'),frame);assert.deepEqual(restored.getAll('select'),selection);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',new URL(url).searchParams,true),source);
 const copied=new URL('http://127.0.0.1:3016/?model=local-ta98');
 api.preserveLocalTa98LayerCatalogueQuery(copied,restored,'local-ta98',true);
 assert.equal(copied.searchParams.get('atlasSource'),source);
 const publicCopy=new URL('https://example.com/?model=local-ta98');
 api.preserveLocalTa98LayerCatalogueQuery(publicCopy,restored,'local-ta98',false);
 assert.equal(publicCopy.searchParams.get('atlasSource'),null);
});
