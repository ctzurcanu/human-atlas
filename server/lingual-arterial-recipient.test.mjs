import assert from 'node:assert/strict';
import test from 'node:test';
import {existsSync,readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const body='.local-models/organ-models/spinal-cord/male-body.json';
if(!existsSync(body)){
 test('Lingual recipient checks require the installed private anatomy',{skip:'Local geometry is intentionally absent'},()=>{});
}else{
const bytes=readFileSync(body),atlas=JSON.parse(bytes),parts=atlas.parts;
const bundle=await build({stdin:{contents:"export * from './app/anatomical-relations';export * from './app/lingual-arterial-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const whole=atlas.concepts.find(c=>c.id==='TA98TONGUE:male:whole');
const selected=whole.elements.map(id=>parts.find(p=>p.id===id));
const id=(role,side='left')=>`TA98ARTERIAL:male:${role}:${side}:finite-wall`;
const byId=new Map(parts.map(p=>[p.id,p]));
const replace=(p,bad)=>parts.map(item=>item===p?bad:item);

test('Both installed finite lingual banks and reciprocal six-branch interfaces are the reviewed geometry',()=>{
 for(const [side,hash] of [['left','6baee9921ce60801c347f7772457783ccc4cc5d9af982dfbc92f55477129c5b4'],['right','1abcbeb01df57e36ab6fc740e768f21b317594ef532c9e537fb22541dcd46b9d']]){
  const family=api.registeredLingualFamily(parts,side);assert.equal(family?.size,6);
  const main=family.get('lingual'),bank=readFileSync(atlas.chunks[main.chunk].url.replace(/^\/local-models\//,'.local-models/'));
  const expected=side==='left'&&main.carotidFiniteWallRevision==='left-four-cell-common-bank-v2'
   ?'e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36':hash;
  assert.equal(createHash('sha256').update(bank).digest('hex'),expected);
  const graph=api.lingualArterialRelations(parts),parent=byId.get(`ZA:External carotid artery.${side[0]}`);
  assert.ok(graph.get(parent.id).some(r=>r.kind==='before'&&r.target===main));
  for(const child of family.values()){
   const upstream=child===main?parent:main;
   assert.ok(graph.get(child.id).some(r=>r.kind==='after'&&r.target===upstream));
   assert.ok(graph.get(upstream.id).some(r=>r.kind==='before'&&r.target===child));
  }
 }
});

test('Whole tongue has eight selectable arterial references without inventing coat-specific perfusion',()=>{
 const rows=api.anatomicalConceptFunctionalRelations(whole,selected,parts);
 assert.equal(rows.filter(r=>r.kind==='arterial').length,8);
 for(const side of ['left','right'])for(const role of ['lingual','dorsal-lingual-0','dorsal-lingual-1','deep-lingual']){
  const artery=byId.get(id(role,side)),row=rows.find(r=>r.target===artery);assert.ok(row);
  assert.match(row.note,/Whole-tongue arterial reference/);assert.match(row.note,/does not establish entry into an individual/);
  const reverse=api.anatomicalRelations(artery,parts).find(r=>r.kind==='supplies'&&r.targetConcept?.id===whole.id);
  assert.ok(reverse);assert.deepEqual([...reverse.targetConcept.elements].sort(),[...whole.elements].sort());
 }
 assert.ok(!api.lingualArterialRelations(parts).get('TA98TONGUE:male:mucosa')?.some(r=>r.kind==='arterial'));
});

test('Tested terminal annuli expose only their ipsilateral muscle and gland recipients',()=>{
 for(const side of ['left','right'])for(const [role,name] of [['suprahyoid','Hyoglossus muscle'],['sublingual','Sublingual gland']]){
  const artery=byId.get(id(role,side)),recipient=byId.get(`ZA:${name}.${side[0]}`);
  const forward=api.anatomicalRelations(recipient,parts).find(r=>r.kind==='arterial'&&r.target===artery);assert.ok(forward);
  const reverse=api.anatomicalRelations(artery,parts).find(r=>r.kind==='supplies'&&r.target===recipient);assert.ok(reverse);
  assert.equal(reverse.note,forward.note);assert.match(reverse.note,/terminal annulus enters the ipsilateral/);
  assert.ok(!api.lingualArterialRelations(parts).get(artery.id).some(r=>r.kind==='supplies'&&r.target.id.endsWith(`.${side==='left'?'r':'l'}`)));
 }
});

test('Wrong sources, transformations, revisions and broken shared rims cannot borrow an installed arterial interface',()=>{
 const artery=byId.get(id('lingual'));
 const variants=[{...artery,sourceId:'other'},{...artery,conceptId:'other'},{...artery,ta98Term:'A12.2.05.020'},
  {...artery,lingualArterialWallRevision:'unreviewed'},{...artery,sectionAssembly:'right-lingual-arterial-family:male'},
  {...artery,geometryOrigin:[0,0,0]},{...artery,sourceOffset:[0,0,.001]},{...artery,geometryScale:.9},
  {...artery,geometryInstance:{rotation:[0,0,0,1],scale:1}},{...artery,suppressed:true},
  {...artery,normalSigns:undefined},{...artery,pairedTissue:false},
  {...artery,tissueDrawGroups:artery.tissueDrawGroups.filter(g=>g.neighbor!==byId.get(id('sublingual')).id)},
  {...artery,tissueDrawGroups:artery.tissueDrawGroups.map(g=>g.neighbor?{...g,start:artery.indexCount}:g)}];
 for(const bad of variants){
  const mixed=replace(artery,bad);assert.equal(api.registeredLingualFamily(mixed,'left'),undefined);
  assert.equal(api.registeredLingualFamily(mixed,'right')?.size,6);
  assert.ok(!api.lingualArterialRelations(mixed).get(bad.id)?.some(r=>['before','after','supplies'].includes(r.kind)));
 }
 const parent=byId.get('ZA:External carotid artery.l');
 assert.equal(api.registeredLingualFamily(replace(parent,{...parent,tissueDrawGroups:[]}), 'left'),undefined);
});

test('Incomplete, moved, unrelated and female recipient assemblies cannot acquire whole-tongue links',()=>{
 assert.deepEqual(api.anatomicalConceptFunctionalRelations(whole,selected.slice(1),parts),[]);
 assert.deepEqual(api.lingualWholeFunctionalRelations({...whole,id:'TA98TONGUE:female:whole'},selected,parts),undefined);
 for(const changed of ['TA98TONGUE:male:mucosa','ZA:Hyoglossus muscle.l']){
  const original=byId.get(changed),bad={...original,geometryOrigin:[0,0,0]},mixed=replace(original,bad);
  assert.deepEqual(api.lingualWholeFunctionalRelations(whole,selected.map(p=>p===original?bad:p),mixed),[]);
 }
 const gland=byId.get('ZA:Sublingual gland.l'),mixed=replace(gland,{...gland,sourceOffset:[0,0,.001]});
 assert.ok(!api.lingualArterialRelations(mixed).get(id('sublingual'))?.some(r=>r.kind==='supplies'));
});

test('Functional controls preserve the installed ordinary-body file',()=>assert.deepEqual(readFileSync(body),bytes));
}
