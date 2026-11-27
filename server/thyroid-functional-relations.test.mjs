import assert from 'node:assert/strict';
import test from 'node:test';
import fs from 'node:fs';
import {build} from 'esbuild';

const compiled=await build({stdin:{contents:"export * from './app/thyroid-functional-relations';export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const load=sex=>JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));

for(const sex of ['male','female'])test(`${sex}: actual thyroid details expose scoped innervation, lymph drainage and reciprocal whole-gland links`,()=>{
 const a=load(sex),before=JSON.stringify(a),id=sex==='male'?'ZA:Thyroid gland':'LOCAL:female:thyroid_gland';
 const c=a.concepts.find(c=>c.id===id),members=a.parts.filter(p=>c.elements.includes(p.id));
 const links=members.length>1?api.anatomicalConceptFunctionalRelations(c,members,a.parts):api.anatomicalRelations(members[0],a.parts);
 const nerves=links.filter(r=>r.kind==='innervation'),nodes=links.filter(r=>r.kind==='lymphaticDrainage');
 assert.equal(nerves.length,4);assert.equal(nodes.length,sex==='male'?7:3);
 assert.equal(links.filter(r=>r.kind==='arterial').length,4);assert.equal(links.filter(r=>r.kind==='venous').length,4);
 const sympathetic=nerves.filter(r=>/sympathetic/i.test(r.target.name));assert.equal(sympathetic.length,2);
 for(const r of sympathetic){
  assert.equal(r.viaRoutes.length,3);
  assert.ok(r.viaRoutes.every(route=>route.via.length===2&&route.via[1]==='Periarterial thyroid plexus'&&!route.modeled));
  assert.deepEqual(r.viaRoutes.map(route=>route.via[0]),['Superior cervical ganglion','Middle cervical ganglion','Inferior cervical ganglion']);
 }
 assert.ok(nerves.filter(r=>/vagus/i.test(r.target.name)).every(r=>r.note.includes('Variable')&&r.note.includes('thyroid-specific branches')));
 for(const r of [...nerves,...nodes]){
  assert.ok(!r.target.suppressed&&!r.target.sourceOffset&&!r.target.id.startsWith('HRA:'));
  for(const targetId of r.targetConcept?.elements??[r.target.id]){
   const target=a.parts.find(p=>p.id===targetId),reverse=api.anatomicalRelations(target,a.parts);
   const back=reverse.find(x=>x.kind===(r.kind==='innervation'?'innervates':'lymphaticTributaries')&&x.targetConcept?.id===id);
   assert.ok(back,targetId);assert.deepEqual(new Set(back.targetConcept.elements),new Set(c.elements));
  }
 }
 if(sex==='male'){
  assert.equal(nodes.filter(r=>r.target.name.startsWith('Thyroid nodes')).length,2);
  assert.ok(nodes.filter(r=>r.target.name.startsWith('Pretracheal')).every(r=>r.via[0]==='Prelaryngeal lymph nodes'&&!r.viaModeled&&r.note.includes('possible')));
  const hierarchy={...c,id:'hierarchy:guest:ta98:A11.3.00.001',children:[]};
  assert.deepEqual(api.anatomicalConceptFunctionalRelations(hierarchy,members,a.parts),links);
  assert.deepEqual(api.anatomicalConceptFunctionalRelations(c,members.slice(1),a.parts),[]);
 }else{
  const group=nodes.find(r=>r.target.name==='Paratracheal lymph nodes');assert.equal(group.targetConcept.elements.length,2);
  assert.ok(nodes.filter(r=>r.target.name.includes('Superior deep')).every(r=>r.note.includes('Upper/lateral')&&r.note.includes('Lower deep cervical nodes')));
 }
 assert.equal(JSON.stringify(a),before,'Functional details never modify geometry, materials or camera state');
});

for(const sex of ['male','female'])test(`${sex}: independent donors, shifted sources, suppressed parts and corrupted identities cannot supply thyroid routes`,()=>{
 const original=load(sex),c=original.concepts.find(c=>c.id===(sex==='male'?'ZA:Thyroid gland':'LOCAL:female:thyroid_gland'));
 const members=original.parts.filter(p=>c.elements.includes(p.id));const links=api.thyroidWholeFunctionalRelations(c,members,original.parts);
 for(const relation of links)for(const change of [{sourceOffset:[0,0,0]},{suppressed:true},{sourceId:'foreign'},{provenance:{mirrored:true}},{id:'HRA:foreign-copy'}]){
  const changed=original.parts.map(p=>p.id===relation.target.id?{...p,...change}:p);
  const selected=changed.filter(p=>c.elements.includes(p.id));
  assert.ok(!api.thyroidWholeFunctionalRelations(c,selected,changed).some(r=>r.target.id===relation.target.id&&r.kind===relation.kind));
 }
 const wrongSex=original.parts.map(p=>p.id===links[0].target.id?{...p,id:sex==='male'?'LOCAL:female:foreign:1':'ZA:foreign',sourceId:'foreign'}:p);
 assert.ok(!api.thyroidWholeFunctionalRelations(c,members,wrongSex).some(r=>r.target.id.includes('foreign')));
});

test('Installed male blood entries require reciprocal shared geometry before claiming a junction',()=>{
 const a=load('male'),c=a.concepts.find(c=>c.id==='ZA:Thyroid gland'),members=a.parts.filter(p=>c.elements.includes(p.id));
 if(!a.thyroidWallAssembly.vascularParts)return;
 const links=api.thyroidWholeFunctionalRelations(c,members,a.parts);
 const joined=links.filter(r=>['arterial','venous'].includes(r.kind)&&r.note.startsWith('Connected'));
 assert.equal(joined.length,6);
 for(const relation of joined){
  const changed=a.parts.map(p=>p.id===relation.target.id?{...p,tissueDrawGroups:p.tissueDrawGroups.filter(g=>!g.neighbor)}:p);
  const selected=changed.filter(p=>c.elements.includes(p.id));
  const r=api.thyroidWholeFunctionalRelations(c,selected,changed).find(r=>r.target.id===relation.target.id&&r.kind===relation.kind);
  assert.ok(r);assert.ok(r.note.includes('not yet installed'),'A missing tissue interface cannot claim a physical join');
 }
});

function femaleFiniteFixture(){
 const male=load('male'),female=load('female');
 const maleWhole=male.concepts.find(c=>c.id==='ZA:Thyroid gland');
 const remap=id=>id.replace('THYROID:male:','THYROID:female:');
 const members=male.parts.filter(p=>maleWhole.elements.includes(p.id)).map(p=>({...p,
  id:remap(p.id),conceptId:remap(p.conceptId),tissueWallAssembly:'thyroid:female',sectionAssembly:'TA98:thyroid-finite:female',
  tissueDrawGroups:p.tissueDrawGroups.map(g=>({...g,...(g.neighbor?{neighbor:remap(g.neighbor)}:{})}))}));
 assert.equal(members.length,11,'Exercise the actual shared eleven-cell tissue records');
 const concept={...female.concepts.find(c=>c.id==='LOCAL:female:thyroid_gland'),elements:members.map(p=>p.id)};
 const parts=[...female.parts.filter(p=>!p.id.startsWith('THYROID:female:')&&p.id!=='LOCAL:female:thyroid_gland:1081'),...members];
 return {concept,members,parts};
}

test('Female finite eleven-cell thyroid retains native functional targets and reciprocal whole-gland routes without claiming vessel joins',()=>{
 const {concept,members,parts}=femaleFiniteFixture(),before=JSON.stringify(parts);
 const links=api.anatomicalConceptFunctionalRelations(concept,members,parts);
 assert.equal(links.filter(r=>r.kind==='innervation').length,4);
 assert.equal(links.filter(r=>r.kind==='lymphaticDrainage').length,3);
 assert.equal(links.filter(r=>r.kind==='arterial').length,4);
 assert.equal(links.filter(r=>r.kind==='venous').length,4);
 const graph=api.thyroidFunctionalRelations(parts);
 for(const relation of links){
  assert.ok(relation.target.id.startsWith('LOCAL:female:'));
  for(const id of relation.targetConcept?.elements??[relation.target.id]){
   const target=parts.find(p=>p.id===id),back=graph.get(id).find(r=>r.targetConcept?.id===concept.id);
   assert.ok(back);assert.deepEqual(new Set(back.targetConcept.elements),new Set(concept.elements));
   assert.ok(api.isReviewedThyroidFunctionalRelation(members[0],target,relation.kind));
   assert.ok(api.isReviewedThyroidFunctionalRelation(target,members[0],back.kind));
  }
 }
 assert.ok(links.filter(r=>['arterial','venous'].includes(r.kind)).every(r=>r.note==='The thyroid–vessel mesh junction is not yet installed.'));
 const hierarchy={...concept,id:'hierarchy:guest:ta98:A11.3.00.001',children:[]};
 assert.deepEqual(api.anatomicalConceptFunctionalRelations(hierarchy,members,parts),links);
 assert.equal(JSON.stringify(parts),before,'No source anatomy, materials or vessel records are modified');
});

test('Female finite routes reject incomplete cells, mixed sex, independent offsets and broken shared interfaces',()=>{
 const {concept,members,parts}=femaleFiniteFixture();
 for(const missing of members){
  const reduced=parts.filter(p=>p!==missing),selected=members.filter(p=>p!==missing);
  assert.deepEqual(api.thyroidWholeFunctionalRelations({...concept,elements:selected.map(p=>p.id)},selected,reduced),[],missing.id);
 }
 for(const change of [{tissueWallAssembly:'thyroid:male'},{sectionAssembly:'TA98:thyroid-finite:male'},
  {sourceOffset:[0,0,0]},{sourceOffset:[.001,0,0]},{suppressed:true},{provenance:{mirrored:true}},
  {normalSigns:undefined},{geometryOrigin:[0,1.56,.001]},{chunk:members[0].chunk+1},
  {pairedTissue:false},{tissueDrawGroups:members[0].tissueDrawGroups.map(g=>({...g,neighbor:undefined}))}]){
  const changed=parts.map(p=>p===members[0]?{...p,...change}:p),selected=changed.filter(p=>concept.elements.includes(p.id));
  assert.deepEqual(api.thyroidWholeFunctionalRelations(concept,selected,changed),[],JSON.stringify(change));
 }
 const mixed=parts.map(p=>p===members[0]?{...p,id:p.id.replace('THYROID:female:','THYROID:male:')}:p);
 const selected=mixed.filter(p=>p.id.startsWith('THYROID:'));
 assert.deepEqual(api.thyroidWholeFunctionalRelations({...concept,elements:selected.map(p=>p.id)},selected,mixed),[]);
 const male=load('male'),maleWhole=male.concepts.find(c=>c.id==='ZA:Thyroid gland');
 const maleMembers=male.parts.filter(p=>maleWhole.elements.includes(p.id));
 assert.deepEqual(api.thyroidWholeFunctionalRelations({...concept,elements:maleWhole.elements},maleMembers,male.parts),[],
  'A known female whole identity cannot activate the complete male assembly');
});

test('Native female whole-gland fallback remains active only when no finite tissue set is present',()=>{
 const a=load('female'),source=a.parts.find(p=>p.id==='LOCAL:female:thyroid_gland:1081');
 assert.ok(source,'The original female source record remains available');
 // Recreate the native-only display state even after a finite set supersedes it.
 const native={...source,suppressed:false};
 const parts=a.parts.filter(p=>!p.id.startsWith('THYROID:female:')).map(p=>p.id===native.id?native:p);
 const concept={id:native.conceptId,name:native.name,ta98Term:'A11.3.00.001',elements:[native.id]};
 const links=api.thyroidWholeFunctionalRelations(concept,[native],parts);
 assert.equal(links.filter(r=>r.kind==='innervation').length,4);
 assert.equal(links.filter(r=>r.kind==='lymphaticDrainage').length,3);
 const {members}=femaleFiniteFixture();
 assert.deepEqual(api.thyroidWholeFunctionalRelations(concept,[native],[...parts,members[0]]),[]);
});
