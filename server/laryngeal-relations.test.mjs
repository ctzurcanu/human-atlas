import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/laryngeal-relations';export * from './app/reviewed-catalogue-identities';export * from './app/anatomical-relations';export * from './app/anatomy-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const read=sex=>JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
test('male correction changes only seven tissue-system fields and retains every bank, material and coordinate descriptor',()=>{
 const a=read('male'),corrected=api.withReviewedCatalogueIdentities(a);
 const changes=corrected.parts.filter((p,i)=>p!==a.parts[i]);assert.equal(changes.length,7);
 assert.strictEqual(corrected.chunks,a.chunks);assert.strictEqual(corrected.materials,a.materials);
 for(const p of changes){const old=a.parts.find(q=>q.id===p.id);assert.equal(p.system,'respiratory');assert.deepEqual({...p,system:old.system},old);}
 assert.strictEqual(api.withLaryngealSystems(corrected),corrected,'classification is idempotent');
 const female=read('female');assert.strictEqual(api.withLaryngealSystems(female),female,'already correct female sources retain their exact catalogue object');
});
test('both ordinary models expose reciprocal cricothyroid, cricoarytenoid and same-side corniculate references',()=>{
 for(const sex of ['male','female']){
  const a=api.withReviewedCatalogueIdentities(read(sex)),graph=api.laryngealRelations(a.parts);
  assert.equal(graph.size,7);assert.equal([...graph.values()].flat().length,14);
  for(const [id,links] of graph){const p=a.parts.find(p=>p.id===id);
   for(const r of links){assert.ok(api.anatomicalRelations(p,a.parts).some(q=>q.target.id===r.target.id&&q.kind===r.kind));
    const inverse=r.kind==='connects'?'connectedBy':r.kind==='connectedBy'?'connects':r.kind;
    assert.ok(graph.get(r.target.id).some(q=>q.target.id===id&&q.kind===inverse));
    assert.ok(!r.target.id.startsWith('HRA:'));
    if(/Corniculate/.test(p.name))assert.equal(/left/.test(p.name),/left/.test(r.target.name));
   }
  }
  const entries=api.hierarchyEntries(a);for(const id of graph.keys()){const entry=entries.find(e=>e.parts.some(p=>p.id===id));assert.ok(entry);assert.equal(api.majorSystemFor(entry),'respiratory');assert.ok(api.laryngealDescription(a.parts.find(p=>p.id===id)));}
 }
});
test('displaced, suppressed, donor, wrong-system and misidentified copies cannot become laryngeal attachments',()=>{
 const a=api.withLaryngealSystems(read('male')),thyroid=a.parts.find(p=>p.id==='ZA:Thyroid cartilage'),cricoid=a.parts.find(p=>p.id==='ZA:Cricoid cartilage');
 for(const delta of [{suppressed:true},{sourceOffset:[0,0,0]},{id:'HRA:thyroid'},{sourceId:'other'},{conceptId:'other'},{name:'Thyroid gland'},{system:'endocrine'}])assert.equal(api.laryngealRelations([{...thyroid,...delta},cricoid]).size,0);
 assert.equal(api.laryngealRelations([...api.withLaryngealSystems(read('male')).parts.filter(p=>p.id==='ZA:Thyroid cartilage'),...read('female').parts.filter(p=>p.id==='LOCAL:female:cricoid_cartilage:1072')]).size,0);
});
