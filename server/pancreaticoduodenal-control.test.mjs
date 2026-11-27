import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const built=await build({stdin:{contents:"export * from './app/ta98-term-bindings'; export * from './app/anatomy-hierarchy'; export * from './app/anatomical-components';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const term='A12.3.12.023',id='control:pancreaticoduodenal',anchor='original:pancreatic-organizer',parent='source:smv';
function fixture(){
 const members=['anterior-superior','anterior-inferior','posterior-inferior'];
 const parts=[...members.map(n=>({id:n,sourceId:`original:${n}`,conceptId:n,sourceGroups:[anchor],name:n,system:'venous',chunk:0,positions:0,normals:36,indices:54,vertexCount:3,indexCount:3,bounds:[[0,1,0],[.01,1.01,.01]]})),
  {id:parent,sourceId:'original:smv',conceptId:parent,ta98Term:'A12.3.12.018',name:'Superior mesenteric vein',system:'venous',chunk:0,positions:0,normals:36,indices:54,vertexCount:3,indexCount:3,bounds:[[0,1,0],[.01,1.01,.01]]},
  {id:'posterior-superior',sourceId:'original:posterior-superior',conceptId:'posterior-superior',sourceGroups:[anchor],ta98Term:'A12.3.12.014',name:'Posterior superior pancreaticoduodenal vein',system:'venous',chunk:0,positions:0,normals:36,indices:54,vertexCount:3,indexCount:3,bounds:[[0,1,0],[.01,1.01,.01]]}];
 const concepts=parts.map(p=>({id:p.conceptId,name:p.name,elements:[p.id],...(p.ta98Term?{ta98Term:p.ta98Term,ta98Kind:'exact'}:{})}));
 concepts.push({id:anchor,name:'Original pancreatic venous source organizer',elements:[...members,'posterior-superior']});
 concepts.push({id,name:'Pancreaticoduodenal veins',elements:members,componentElements:members,partOf:parent,ta98Term:term,ta98Kind:'exact',sourceCollectionId:id,coverageLimitation:'Named source collection; fit, walls and drainage continuity remain unverified.'});
 return {sex:'male',parts,concepts,chunks:[{url:'source.bin',bytes:66}],materials:{},defaultHidden:[],
  sourceCollectionIdentities:[{id,term,sourceId:'original:pancreatic-organizer',sourceAnchorConceptId:anchor,taParentTerm:'A12.3.12.018',taParentConceptId:parent,identityOnly:true,members:parts.filter(p=>members.includes(p.id)).map(p=>({partId:p.id,sourceId:p.sourceId,conceptId:p.conceptId}))}],
  ta98TermBindings:[{term,conceptId:id,partIds:members,kind:'exact',representation:'source-collection-identity',sourceCollectionId:id}]};
}
test('anatomical TA98 parent is separate from immutable source provenance, with exactly three pancreaticoduodenal members',()=>{
 const a=fixture(),before=structuredClone(a),parts=api.catalogueTa98TermParts(a).get(term);
 assert.deepEqual(parts.map(p=>p.id),['anterior-superior','anterior-inferior','posterior-inferior']);
 assert.ok(!parts.some(p=>p.id==='posterior-superior'));assert.deepEqual(a,before);
 const entry=api.hierarchyEntries(a).find(e=>e.id===id);assert.ok(entry);
 assert.equal(entry.taParentConcept,parent);assert.ok(entry.ancestry.some(n=>n.toLowerCase()==='superior mesenteric vein'));
 const c=a.concepts.find(c=>c.id===id);
 assert.deepEqual(api.anatomicalComponentIds(c,[],new Map(a.parts.map(p=>[p.id,p])),a.concepts),c.elements);
 assert.ok(a.parts.filter(p=>c.elements.includes(p.id)).every(p=>p.sourceGroups.includes(anchor)&&!p.sourceGroups.includes(parent)));
});
test('separate-parent collections reject forged lineage, wrong blood system, hidden parents, provenance edits and cycles',()=>{
 const mutations=[
  a=>a.sourceCollectionIdentities[0].taParentTerm='A12.3.12.014',
  a=>delete a.sourceCollectionIdentities[0].taParentConceptId,
  a=>a.sourceCollectionIdentities[0].taParentConceptId=anchor,
  a=>a.concepts.find(c=>c.id===id).partOf=anchor,
  a=>a.concepts.find(c=>c.id===parent).ta98Kind='parent',
  a=>a.parts.find(p=>p.id===parent).system='arterial',
  a=>a.defaultHidden.push(parent),
  a=>a.parts[0].sourceGroups=[],
  a=>a.concepts.find(c=>c.id===parent).partOf=id,
 ];
 for(const mutate of mutations){const a=fixture();mutate(a);assert.equal(api.catalogueTa98TermParts(a).has(term),false);}
});
const candidate='.local-models/geometry-workbench/biodigital-male-circulation-v1/blood-pancreaticoduodenal-identities-v2435/circulation-atlas.json';
test('actual local candidate preserves all physical anatomy and exposes the reviewed three-mesh control',{skip:!fs.existsSync(candidate)},()=>{
 const a=JSON.parse(fs.readFileSync(candidate,'utf8')),old=JSON.parse(fs.readFileSync('.local-models/geometry-workbench/biodigital-male-circulation-v1/blood-temporal-skin-course-v2432/circulation-atlas.json','utf8'));
 assert.deepEqual(a.parts,old.parts);assert.deepEqual(a.chunks,old.chunks);assert.deepEqual(a.materials,old.materials);
 const control='TA98BLOOD:male:A12.3.12.023',parts=api.catalogueTa98TermParts(a).get(term);assert.equal(parts.length,3);
 assert.ok(parts.every(p=>/^(?:Anterior superior|Anterior inferior|Posterior inferior) pancreaticoduodenal vein$/.test(p.name)));
 const entry=api.hierarchyEntries(a).find(e=>e.id===control);assert.equal(entry.parts.length,3);
 assert.equal(entry.taParentConcept,'BDM:male:human_20_male_cardiovascular_system-superior_mesenteric_vein_ID');
 assert.ok(entry.ancestry.some(n=>n.toLowerCase()==='superior mesenteric vein'));
});
