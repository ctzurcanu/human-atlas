import test from 'node:test';import assert from 'node:assert/strict';import {readFile} from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const kinds=new Set(['lymphaticDrainage','lymphaticTributaries']);
for(const sex of ['female','male'])test(`${sex}: whole bowel and all compatible finite cells have reciprocal mesenteric lymph routes`,async()=>{
 const ordinary=JSON.parse(await readFile(sex==='female'?'.local-models/ta98-female-runtime.json':'.local-models/ta98-runtime.json'));
 const review=JSON.parse(await readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`));
 const expected=sex==='female'?['LOCAL:female:superior_mesenteric_lymph_nodes:174']:['ZA:Juxta-intestinal mesenteric nodes','ZA:Central superior mesenteric nodes'];
 const whole=sex==='female'?['LOCAL:female:small_intestine_jejunum:182','LOCAL:female:small_intestine_ileum:181']:['ZA:Jejunum'];
 const check=(source,parts)=>{const links=anatomicalRelations(source,parts).filter(r=>r.kind==='lymphaticDrainage');assert.deepEqual(new Set(links.map(r=>r.target.id)),new Set(expected));for(const link of links){assert.equal(link.target.system,'lymphatic');assert.ok(link.via.includes('Intestinal wall lymphatic plexuses'));assert.match(link.note,/not attach.*directly/);assert.match(link.note,/unverified/);assert.ok(anatomicalRelations(link.target,parts).some(r=>r.kind==='lymphaticTributaries'&&r.target.id===source.id&&r.note===link.note));}};
 for(const id of whole)check(ordinary.parts.find(p=>p.id===id),ordinary.parts);
 const cells=review.parts.filter(p=>p.bowelWallSource&&/:BOWELWALL:(jejunum|ileum):/.test(p.id));assert.equal(cells.length,24);
 for(const p of cells){if(sex==='male'&&p.id.includes(':ileum:'))assert.equal(anatomicalRelations(p,review.parts).filter(r=>kinds.has(r.kind)).length,0);else check(p,review.parts);}
 // A grouped bowel collector must not silently include paracolic nodes.
 for(const p of cells)assert.ok(!anatomicalRelations(p,review.parts).some(r=>r.kind==='lymphaticDrainage'&&r.target.name.includes('Paracolic')));
});
test('lymph routes reject donor, wrong source, suppressed nodes and mismatched frames in both directions',async()=>{
 const atlas=JSON.parse(await readFile('.local-models/ta98-review/bowel-wall-main-prototype/female/atlas.json')),source=atlas.parts.find(p=>p.id.endsWith(':BOWELWALL:ileum:epithelium')),node=atlas.parts.find(p=>p.id==='LOCAL:female:superior_mesenteric_lymph_nodes:174');
 for(const bad of [{...source,id:'HRA:donor-ileum'},{...source,sourceId:'wrong'},{...source,pairedTissue:undefined},{...source,suppressed:true}]){const parts=atlas.parts.map(p=>p.id===source.id?bad:p);assert.ok(!anatomicalRelations(bad,parts).some(r=>kinds.has(r.kind)));assert.ok(!anatomicalRelations(node,parts).some(r=>kinds.has(r.kind)&&r.target.id===bad.id));}
 for(const bad of [{...node,id:'ZA:Superior mesenteric nodes'},{...node,suppressed:true},{...node,system:'venous'}]){const parts=atlas.parts.map(p=>p.id===node.id?bad:p);assert.ok(!anatomicalRelations(source,parts).some(r=>kinds.has(r.kind)));}
});
