import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';

const bundled=await build({stdin:{contents:"export * from './app/blood-reference-identities';export * from './app/reviewed-catalogue-identities';export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundled.outputFiles[0].contents).toString('base64'));
const metadata=JSON.parse(readFileSync('app/data/ta98-metadata.json'));
const vertebral='A12.2.08.002';
const thalamicIds=['FJ1675','FJ1720'];
const nativeIds=['ZA:Vertebral artery.l','ZA:Vertebral artery.r'];

function fixture(){
 const parts=[
  {id:'FJ1675',name:'Right thalamogeniculate artery',conceptId:'FMA50670',system:'arterial',ta98Term:'A12.2.07.090'},
  {id:'FJ1720',name:'Right thalamoperforating artery',conceptId:'FMA50664',system:'arterial',ta98Term:'A12.2.07.086'},
  {id:nativeIds[0],name:'Vertebral artery (left)',conceptId:nativeIds[0],system:'arterial'},
  {id:nativeIds[1],name:'Vertebral artery (right)',conceptId:nativeIds[1],system:'arterial',ta98Term:vertebral},
 ];
 const concepts=[
  {id:'FMA3958',name:'right vertebral artery',elements:[...thalamicIds,...nativeIds],componentElements:[...thalamicIds,...nativeIds]},
  {id:'FMA4066',name:'left vertebral artery',elements:[...thalamicIds,...nativeIds]},
  {id:'FMA50670',name:'Right thalamogeniculate artery',elements:['FJ1675']},
  {id:'FMA50664',name:'Right thalamoperforating artery',elements:['FJ1720']},
  ...nativeIds.map(id=>({id,name:id.slice(3),elements:[id]})),
 ];
 return {id:'local-ta98',name:'Minimal reviewed blood identities',parts,concepts,chunks:[],materials:{}};
}
function freeze(value){
 if(value&&typeof value==='object'){Object.values(value).forEach(freeze);Object.freeze(value);}
 return value;
}

test('real stale vertebral aggregates drop explicitly identified thalamic carriers and retain native vertebral sources',()=>{
 const atlas=fixture();
 for(const id of ['FMA3958','FMA4066'])assert.deepEqual(metadata.byConceptMatch[id],{term:vertebral,kind:'exact',method:'name'});
 assert.equal(metadata.byConceptMatch.FMA50670.term,'A12.2.07.090');
 assert.equal(metadata.byConceptMatch.FMA50664.term,'A12.2.07.086');
 const result=api.bloodReferenceConcepts(atlas,atlas.concepts);
 for(const id of ['FMA3958','FMA4066'])assert.deepEqual(result.find(c=>c.id===id).elements,nativeIds,id);
 assert.deepEqual(result[0].componentElements,nativeIds,'component selection cannot retain the contradicted carriers');
 assert.equal(result[1].componentElements,undefined,'optional component list is not invented');
 for(let i=2;i<atlas.concepts.length;i++)assert.strictEqual(result[i],atlas.concepts[i],'unrelated/source concepts preserve their record identity');
 assert.strictEqual(atlas.parts[2].ta98Term,undefined,'the true untagged native source is retained without authoring a new tag');
});

test('an explicit authored exact mapping and a bounded authored carrier keep their contracts',()=>{
 const atlas=fixture();
 const exact={...atlas.concepts[0],ta98Term:vertebral,ta98Kind:'exact'};
 const carrier={id:'BOUNDED:vertebral:surface',name:'Bounded authored surface',conceptId:'BOUNDED:vertebral:surface',system:'arterial',ta98Term:'A12.2.07.090',surfaceOnly:true,sourceFaceIds:[4,5]};
 const bounded={...atlas.concepts[1],ta98Term:vertebral,ta98Kind:'exact',elements:[carrier.id],componentElements:[carrier.id]};
 const input={...atlas,parts:[...atlas.parts,carrier],concepts:[exact,bounded]};
 const result=api.bloodReferenceConcepts(input,input.concepts);
 assert.strictEqual(result,input.concepts);
 assert.strictEqual(result[0],exact);assert.strictEqual(result[1],bounded);
 assert.deepEqual(result[1].elements,[carrier.id]);assert.deepEqual(carrier.sourceFaceIds,[4,5]);
});

test('a genuine authority descendant remains in the parent blood assembly',()=>{
 const atlas=fixture(),id='ZA:Anterior spinal artery.l';
 assert.ok(metadata.byConceptPath[id].includes(vertebral),'fixture uses an actual authority descendant');
 assert.equal(metadata.byConceptMatch[id].kind,'exact');
 const descendant={id,name:'Anterior spinal artery (left)',conceptId:id,system:'arterial',ta98Term:metadata.byConceptMatch[id].term};
 const concept={...atlas.concepts[0],elements:['FJ1675',id],componentElements:['FJ1675',id]};
 const result=api.bloodReferenceConcepts({...atlas,parts:[...atlas.parts,descendant]},[concept]);
 assert.deepEqual(result[0].elements,[id]);assert.deepEqual(result[0].componentElements,[id]);
});

test('missing parts, missing source identity and conflicting explicit evidence are preserved',()=>{
 const atlas=fixture();
 const missing={id:'UNVERIFIED',name:'Unverified vessel',conceptId:'UNREGISTERED:source',system:'arterial',ta98Term:'A12.2.07.090'};
 const conflicting={id:'CONFLICTING',name:'Conflicting source identity',conceptId:'FMA50670',system:'arterial',ta98Term:'A12.2.07.086'};
 const ids=['ABSENT-PART',missing.id,conflicting.id];
 const concept={...atlas.concepts[0],elements:ids,componentElements:[...ids]};
 const result=api.bloodReferenceConcepts({...atlas,parts:[...atlas.parts,missing,conflicting]},[concept]);
 assert.strictEqual(result[0],concept);assert.deepEqual(result[0].elements,ids);
});

test('identity cleanup does not mutate inputs and is idempotent',()=>{
 const atlas=freeze(fixture()),snapshot=structuredClone(atlas);
 const once=api.bloodReferenceConcepts(atlas,atlas.concepts);
 assert.deepEqual(atlas,snapshot);assert.notStrictEqual(once,atlas.concepts);
 assert.strictEqual(api.bloodReferenceConcepts(atlas,once),once,'a second pass does not allocate or alter records');
});

test('reviewed-catalogue integration removes thalamic parts from direct vertebral selection',()=>{
 const atlas=freeze(fixture()),snapshot=structuredClone(atlas),effective=api.withReviewedCatalogueIdentities(atlas);
 const hierarchy=api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json')));
 const roots=api.resolveGuestHierarchy(effective,hierarchy),nodes=new Map();
 const walk=node=>{if(nodes.has(node.id))return;nodes.set(node.id,node);node.children.forEach(walk);};roots.forEach(walk);
 const direct=nodes.get(vertebral)?.directParts.map(p=>p.id).sort();
 assert.deepEqual(direct,[...nativeIds].sort());
 for(const id of thalamicIds)assert.ok(!direct.includes(id),id);
 assert.strictEqual(effective.parts,atlas.parts,'only contradicted aggregate references change');
 assert.deepEqual(atlas,snapshot);
 assert.deepEqual(api.withReviewedCatalogueIdentities(effective),effective);
});
