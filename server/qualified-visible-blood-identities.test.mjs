import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {existsSync,readFileSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';

const definitions=JSON.parse(readFileSync('public/assets/ta98.json'));
const metadata=JSON.parse(readFileSync('app/data/ta98-metadata.json'));
const rules=JSON.parse(readFileSync('scripts/ta98-reviewed-mappings.json'));
const groups=new Map([
 ['A12.2.09.043',['O3M:Anterior ulnar recurrent artery.l','O3M:Anterior ulnar recurrent artery.r']],
 ['A12.2.09.044',['O3M:Posterior ulnar recurrent artery.l','O3M:Posterior ulnar recurrent artery.r']],
 ['A12.2.12.055',['ZA:Anterior inferior pancreaticoduodenal artery']],
 ['A12.2.01.212',['ZA:Inferior lingular artery of left lung']],
 ['A12.2.01.213',['ZA:Superior lingular artery of left lung']],
]);
const expectedParts=new Map([...groups].map(([code,ids])=>[code,code.endsWith('.212')?
 [...ids,'TA98CIRCULATION:male:inferior-lingular-artery:source-component-1','TA98CIRCULATION:male:inferior-lingular-artery:proximal']:
 code.endsWith('.213')?[...ids,'TA98CIRCULATION:male:superior-lingular-artery:source-component-1','TA98CIRCULATION:male:superior-lingular-artery:source-component-2','TA98CIRCULATION:male:superior-lingular-artery:proximal']:ids]));
const bundled=await build({stdin:{contents:`
 export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities';
 export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations';
 export {parseGuestHierarchy,resolveGuestHierarchy} from './app/guest-hierarchy';
 `,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundled.outputFiles[0].contents).toString('base64'));
const definition=api.parseGuestHierarchy(definitions);
const sha=value=>createHash('sha256').update(value).digest('hex');
function nodes(atlas){
 const result=new Map();
 const visit=node=>{if(result.has(node.id))return;result.set(node.id,node);node.children.forEach(visit);};
 api.resolveGuestHierarchy(atlas,definition).forEach(visit);return result;
}

test('reviewed mappings are exact scoped source identities with preserved qualifications',()=>{
 for(const [term,ids] of groups){
  for(const id of ids){
   assert.deepEqual(metadata.reviewedCatalogueIdentities[id],{term,kind:'exact'});
   assert.deepEqual(metadata.byConceptMatch[id],{term,kind:'exact',method:'reviewed'});
   assert.equal(metadata.byConcept[id].ta98,term);
   assert.ok(metadata.byConceptPath[id].length);
   const rule=rules.find(rule=>rule.conceptIds?.includes(id));
   assert.equal(rule.term,term);assert.equal(rule.kind,'exact');
   assert.match(rule.note,/unverified|not.*finite|not.*accepted/i);
   const node=definitions.nodes.find(node=>node.id===term);
   assert.ok(node.matches.includes(id));assert.match(node.description,/source identity|lingular assembly/);
  }
 }
 assert.equal(metadata.byConceptMatch['Unqualified anterior branch'],undefined);
 assert.equal(definitions.nodes.find(node=>node.id==='A12.2.09.043').name,'anterior branch');
 assert.equal(definitions.nodes.find(node=>node.id==='A12.2.12.055').name,'anterior branch');
});

test('source scoped identities preserve physical ownership and never bind a generic homonym',()=>{
 const concepts=[...groups].flatMap(([code,ids])=>ids.map(id=>({id,name:id.slice(id.indexOf(':')+1),elements:expectedParts.get(code).filter(part=>ids.length===1||part===id)})));
 const parts=concepts.flatMap(concept=>concept.elements.map(id=>({id,name:id,conceptId:concept.id,system:'arterial'})));
 concepts.push({id:'Unqualified anterior branch',name:'Anterior branch',elements:['GENERIC']});
 parts.push({id:'GENERIC',name:'Anterior branch',conceptId:'Unqualified anterior branch',system:'arterial'});
 const raw={id:'local-ta98',name:'Scoped source fixture',parts,concepts,chunks:[],materials:{}};
 const before=structuredClone(raw),effective=api.withReviewedCatalogueIdentities(raw),tree=nodes(effective);
 for(const [term,ids] of expectedParts)assert.deepEqual(tree.get(term).directParts.map(part=>part.id).sort(),[...ids].sort());
 assert.deepEqual(raw,before);assert.strictEqual(effective.parts,raw.parts);
 assert.ok(effective.concepts.every(concept=>!concept.partOf),'no physical parent is invented');
 assert.equal(effective.concepts.find(concept=>concept.id==='Unqualified anterior branch').ta98Term,undefined);
});

const current='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
const oldSnapshot='.local-models/geometry-workbench/circulation-census-v1/ordinary-view-v1/resolved-male.json';
const oldCensus='.local-models/geometry-workbench/circulation-census-v1/ordinary-view-v1/census-male.json';
test('actual ordinary male resolver adds exactly five controls without changing any source geometry',{
 skip:!existsSync(current)||!existsSync(oldSnapshot)||!existsSync(oldCensus)?'Private current catalogue/census unavailable; standalone scoped tests still run.':false,
},()=>{
 const bytes=readFileSync(current);assert.equal(sha(bytes),'346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6');
 const raw=JSON.parse(bytes),rawBefore=JSON.stringify(raw);
 const prepared=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(raw)),false);
 const tree=nodes(prepared),baseline=JSON.parse(readFileSync(oldSnapshot)),census=JSON.parse(readFileSync(oldCensus));
 assert.equal(baseline.sourceSha256,sha(bytes));
 const maleCodes=census.terms.filter(term=>term.male_applicable).map(term=>term.id);
 const oldAvailable=new Set(baseline.resolved.filter(node=>node.directPartIds.length).map(node=>node.id));
 const nowAvailable=new Set([...tree.values()].filter(node=>node.directParts.length).map(node=>node.id));
 const gained=maleCodes.filter(code=>nowAvailable.has(code)&&!oldAvailable.has(code));
 const lost=maleCodes.filter(code=>oldAvailable.has(code)&&!nowAvailable.has(code));
 assert.deepEqual(gained.sort(),[...groups.keys()].sort());assert.deepEqual(lost,[]);
 for(const [term,ids] of expectedParts)assert.deepEqual(tree.get(term).directParts.map(part=>part.id).sort(),[...ids].sort());
 assert.equal(tree.get('A12.2.12.056')?.directParts.length??0,0,'posterior pancreaticoduodenal branch stays pending');
 for(const id of ['O3M:Anterior ulnar recurrent artery.l','O3M:Posterior ulnar recurrent artery.l'])assert.equal(raw.parts.find(part=>part.id===id).provenance.mirrored,true);
 assert.equal(JSON.stringify(raw),rawBefore,'ordinary resolver never mutates raw anatomy');
 assert.equal(sha(readFileSync(current)),sha(bytes));
 const female='.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';
 if(existsSync(female))assert.equal(sha(readFileSync(female)),'8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e');
 console.log(JSON.stringify({ordinaryNewControls:gained,ordinaryLostControls:lost,sourceGeometryChanged:false}));
});
