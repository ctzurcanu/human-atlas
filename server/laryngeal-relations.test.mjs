import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/laryngeal-relations';export * from './app/reviewed-catalogue-identities';export * from './app/anatomical-relations';export * from './app/anatomy-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const read=sex=>JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
test('male classification changes only active native tissue-system fields and retains every geometry descriptor',()=>{
 const a=read('male'),corrected=api.withReviewedCatalogueIdentities(a);
 const native=['Thyroid cartilage','Cricoid cartilage','Arytenoid cartilage.l','Arytenoid cartilage.r','Corniculate cartilage.l','Corniculate cartilage.r','Median cricothyroid ligament'].map(id=>'ZA:'+id);
 const expected=a.parts.filter(p=>native.includes(p.id)&&!p.suppressed&&p.system!=='respiratory').map(p=>p.id).sort();
 const changes=corrected.parts.filter((p,i)=>p!==a.parts[i]);assert.deepEqual(changes.map(p=>p.id).sort(),expected);
 assert.strictEqual(corrected.chunks,a.chunks);assert.strictEqual(corrected.materials,a.materials);
 for(const p of changes){const old=a.parts.find(q=>q.id===p.id);assert.equal(p.system,'respiratory');assert.deepEqual({...p,system:old.system},old);}
 assert.strictEqual(api.withLaryngealSystems(corrected),corrected,'classification is idempotent');
 const female=read('female');assert.strictEqual(api.withLaryngealSystems(female),female,'already correct female sources retain their exact catalogue object');
});
test('both ordinary models expose reciprocal cricothyroid, cricoarytenoid and same-side corniculate references',()=>{
 for(const sex of ['male','female']){
  const a=api.withReviewedCatalogueIdentities(read(sex)),graph=api.laryngealRelations(a.parts);
  const replaced=a.parts.some(p=>p.id===`CRICOID:${sex}:arch`&&!p.suppressed);
  assert.equal(graph.size,replaced?11:7);assert.equal([...graph.values()].flat().length,replaced?16:14);
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
 const a=api.withLaryngealSystems(read('male')),thyroid=a.parts.find(p=>p.id==='ZA:Thyroid cartilage'),cricoid={...a.parts.find(p=>p.id==='ZA:Cricoid cartilage'),suppressed:false,system:'respiratory'};
 assert.equal(api.laryngealRelations([thyroid,cricoid]).size,2,'The literal retained native identities are a positive control');
 for(const delta of [{suppressed:true},{sourceOffset:[0,0,0]},{id:'HRA:thyroid'},{sourceId:'other'},{conceptId:'other'},{name:'Thyroid gland'},{system:'endocrine'}])assert.equal(api.laryngealRelations([{...thyroid,...delta},cricoid]).size,0);
 assert.equal(api.laryngealRelations([...api.withLaryngealSystems(read('male')).parts.filter(p=>p.id==='ZA:Thyroid cartilage'),...read('female').parts.filter(p=>p.id==='LOCAL:female:cricoid_cartilage:1072')]).size,0);
});

test('actual staged male regions retain reciprocal source-specific joint and ligament links',()=>{
 const a=api.withReviewedCatalogueIdentities(JSON.parse(fs.readFileSync('.local-models/ta98-review/male-cricoid-packed-v2/male/atlas.json'))),graph=api.laryngealRelations(a.parts);
 assert.equal(graph.size,11);assert.equal([...graph.values()].flat().length,16);assert.equal(graph.has('ZA:Cricoid cartilage'),false);
 assert.equal(graph.get('CRICOID:male:arch')[0].target.id,'ZA:Median cricothyroid ligament');
 for(const side of ['left','right']){
  const id=`CRICOID:male:arytenoid-articular:${side}`,ary=`ZA:Arytenoid cartilage.${side==='left'?'l':'r'}`;
  assert.equal(graph.get(id)[0].target.id,ary);assert.ok(graph.get(ary).some(r=>r.target.id===id));
 }
});

test('collective cricoid descriptions refer to the selected whole or surface pair',()=>{
 const a=api.withReviewedCatalogueIdentities(JSON.parse(fs.readFileSync('.local-models/ta98-review/male-cricoid-packed-v2/male/atlas.json')));
 const ids=a.cricoidPartition.partIds,members=ids.map(id=>a.parts.find(p=>p.id===id));
 assert.match(api.laryngealSelectionDescription(members,a.parts).text,/complete ring/);
 assert.match(api.laryngealSelectionDescription(members.slice(1),a.parts).text,/all four articular/);
 assert.match(api.laryngealSelectionDescription(members.slice(2,4),a.parts).text,/paired superior/);
 assert.match(api.laryngealSelectionDescription(members.slice(4,6),a.parts).text,/paired lateral/);
 assert.equal(api.laryngealSelectionDescription(members.slice(0,2),a.parts),undefined);
 assert.equal(api.laryngealSelectionDescription(members,a.parts.filter(p=>p.id!==ids[2])),undefined);
 assert.equal(api.laryngealSelectionDescription([...members.slice(0,5),members[0]],a.parts),undefined);
});

function partitionedCricoid(sex){
 const a=api.withLaryngealSystems(read(sex)),sourceId=sex==='female'?'LOCAL:female:cricoid_cartilage:1072':'ZA:Cricoid cartilage';
 const source=a.parts.find(p=>p.id===sourceId),prefix=`CRICOID:${sex}:`;
 const suffixes=['arch','lamina-remainder','arytenoid-articular:left','arytenoid-articular:right','thyroid-articular:left','thyroid-articular:right'];
 const ids=suffixes.map(s=>prefix+s),laminaIds=ids.slice(1);
 const patches=suffixes.map((suffix,index)=>({...source,id:ids[index],conceptId:ids[index],name:suffix,system:'respiratory',suppressed:false,defaultHidden:false,
  cricoidSourcePartId:sourceId,cricoidRegion:index===0?'arch':'lamina',
  cricoidSurface:index===0?undefined:index===1?'remainder':suffix.replace(':','-'),
  sectionSolid:index===0,sectionAssembly:`TA98:cricoid:${sex}:${index===0?'arch':'lamina'}`,
  sectionClosureParts:index===0?undefined:[...laminaIds],geometryOrigin:[0,1.5,0]}));
 const parts=a.parts.filter(p=>!p.id.startsWith(prefix)).map(p=>p.id===sourceId?{...p,suppressed:true,replacementParts:ids}:p);
 return {parts:[...parts,...patches],patches,sourceId};
}
test('complete cricoid partitions inherit only the joint or ligament of their actual surface region',()=>{
 for(const sex of ['female','male']){
  const {parts,patches,sourceId}=partitionedCricoid(sex),graph=api.laryngealRelations(parts);
  assert.equal(graph.has(sourceId),false,'The retired whole is not another visible joint target');
  assert.equal(graph.has(patches[1].id),false,'Lamina remainder does not stand in for either articular surface');
  const arch=graph.get(patches[0].id);assert.equal(arch.length,1);assert.equal(arch[0].kind,'connectedBy');
  assert.match(arch[0].target.name,/cricothyroid ligament/i);
  for(const index of [2,3,4,5]){
   const links=graph.get(patches[index].id);assert.equal(links.length,1);assert.equal(links[0].kind,'articulates');
   if(index<4)assert.match(links[0].target.name,index===2?/Arytenoid.*left/i:/Arytenoid.*right/i);
   else assert.equal(links[0].target.name,'Thyroid cartilage');
   assert.ok(graph.get(links[0].target.id).some(r=>r.target.id===patches[index].id&&r.kind==='articulates'));
   assert.ok(api.isReviewedLaryngealRelation(patches[index],links[0].target,'articulates'));
   assert.match(links[0].note,/provisional/);
  }
  assert.match(api.laryngealDescription(patches[0]).text,/anterior arch/);
  assert.match(api.laryngealDescription(patches[2]).text,/surface region, not a separate solid/);
 }
});
test('missing, displaced, mixed-frame or inconsistent cricoid patches cannot inherit whole-source interfaces',()=>{
 const check=mutate=>{
  const f=partitionedCricoid('female');mutate(f);const graph=api.laryngealRelations(f.parts);
  assert.ok(f.patches.every(p=>!graph.has(p.id)));
 };
 check(f=>f.parts=f.parts.filter(p=>p.id!==f.patches[2].id));
 check(f=>f.patches[2].sourceOffset=[0,0,0]);
 check(f=>f.patches[2].cricoidSourcePartId='ZA:Cricoid cartilage');
 check(f=>f.patches[2].cricoidSurface='arytenoid-articular-right');
 check(f=>f.patches[2].normals+=6);
 check(f=>f.patches[2].sectionClosureParts.pop());
 check(f=>f.patches[2].sectionSolid=true);
 check(f=>f.parts.find(p=>p.id===f.sourceId).sourceId='wrong-source');
 check(f=>f.parts.push({...f.patches[2],id:'unexpected-cricoid-patch'}));
});
