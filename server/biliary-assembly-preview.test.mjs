import test from 'node:test';import assert from 'node:assert/strict';import {readFileSync,existsSync} from 'node:fs';import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/local-ta98-layer-catalogue';export * from './app/guest-hierarchy';export * from './app/digestive-duct-relations';export * from './app/viewer-state';export * from './app/view-part-codes';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));const read=p=>JSON.parse(readFileSync(p));
test('paired biliary preview stays opt-in and never changes public, male, source or prototype routing',()=>{
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('hepatic=assembly'),true),'/local-models/ta98-review/hepatic-biliary-preview/female/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('hepatic=assembly'),true),'/local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json');
 for(const q of ['hepatic=assembly&bowel=source','hepatic=assembly&delivery=0'])assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams(q),true),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('hepatic=assembly&bowel=prototype'),true),'/local-models/ta98-review/bowel-wall-main-prototype/female/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('hepatic=assembly'),false),undefined);
});
test('actual finite walls preserve both source identities, TA98 controls, reciprocal neck rims, relations and saved camera',{skip:!existsSync('reports/ta98-female/biliary-assembly-preview-validation.json')},()=>{
 const proof=read('reports/ta98-female/biliary-assembly-preview-validation.json');assert.equal(proof.deliveryGatePass,true);
 const a=read('.local-models/ta98-review/hepatic-boundary-preview/female/atlas.json'),b=read('.local-models/ta98-review/hepatic-biliary-preview/female/atlas.json');
 const ids=['LOCAL:female:gallbladder:1030','LOCAL:female:bile_duct:1031'];
 assert.deepEqual(a.concepts,b.concepts);assert.deepEqual(a.defaultHidden,b.defaultHidden);assert.deepEqual(a.parts.map(p=>p.id),b.parts.map(p=>p.id));
 assert.deepEqual(a.parts.filter(p=>!ids.includes(p.id)),b.parts.filter(p=>!ids.includes(p.id)));
 const parts=ids.map(id=>b.parts.find(p=>p.id===id));assert.ok(parts.every(p=>p.pairedTissue&&p.sectionSolid&&!p.suppressed));
 for(const p of parts)assert.deepEqual(b.materials[p.material].color,p.tissueAppearance.color.map(c=>Math.round(c*255)),'outer and section material colors share the same authored tissue field');
 assert.equal(parts[0].positions,parts[1].positions);assert.equal(parts[0].normals,parts[1].normals);assert.equal(parts[0].chunk,parts[1].chunk);
 for(const p of parts){assert.deepEqual(api.resolveSelection(b,[p.id]),[p.id]);const rim=p.tissueDrawGroups.filter(g=>g.neighbor);assert.equal(rim.length,1);assert.equal(rim[0].count,42);assert.equal(rim[0].neighbor,ids.find(id=>id!==p.id));}
 const guest=api.parseGuestHierarchy(read('public/assets/ta98.json'));const flat=nodes=>nodes.flatMap(n=>[{id:n.id,parts:n.parts.map(p=>p.id)},...flat(n.children)]);
 assert.deepEqual(flat(api.resolveGuestHierarchy(a,guest)),flat(api.resolveGuestHierarchy(b,guest)));
 const links=graph=>[...graph].map(([id,rows])=>[id,rows.map(r=>({kind:r.kind,target:r.target.id,note:r.note,via:r.via}))]);assert.deepEqual(links(api.digestiveDuctRelations(a.parts)),links(api.digestiveDuctRelations(b.parts)));
 const frame=[.4,.5,.9,.5,.5,.5,.5,.45,.3,.5,1,.5],q=new URLSearchParams({model:'local-female-ta98',hepatic:'assembly',frame:frame.join(','),isolate:'1'});ids.forEach(id=>q.append('select',id));
 const state=api.readViewUrl(q.toString(),b,{selected:[],visible:[],view:'three-quarter'});assert.deepEqual(state.selected,ids);assert.deepEqual(state.camera,frame);
 assert.equal(new URL(api.viewUrl('http://localhost:3016/','local-female-ta98',state,frame,b)).searchParams.get('hepatic'),'assembly');
 assert.equal(api.partFingerprint(a),api.partFingerprint(b));
});
