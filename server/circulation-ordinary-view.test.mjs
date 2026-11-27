import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
import {circulationCatalogueSnapshot,ordinaryCirculationCatalogue} from './circulation-effective-catalogue.mjs';

const source='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
const available=fs.existsSync(source);
const skip={skip:!available&&'Private local source is not redistributed.'};
const built=await build({stdin:{contents:`
 export * from './app/guest-hierarchy';
 export * from './app/reviewed-catalogue-identities';
 export * from './app/preferred-source-representations';
 export * from './app/carotid-surface-features';
 export * from './app/anatomy-hierarchy';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const hash=x=>createHash('sha256').update(x).digest('hex');

test('ordinary census excludes hidden donors while retaining installed finite and qualified controls',skip,()=>{
 const raw=fs.readFileSync(source),atlas=JSON.parse(raw),before=hash(JSON.stringify(atlas));
 const metadata=JSON.parse(fs.readFileSync('app/data/ta98-metadata.json'));
 const definition=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
 const snapshot=circulationCatalogueSnapshot(atlas,metadata,api,definition);
 const visibility=new Map(snapshot.effectivePartVisibility.map(p=>[p.id,p]));
 assert.equal(visibility.size,atlas.parts.length);
 for(const id of ['FJ2884','FJ2885']){
  assert.equal(visibility.get(id).rawSuppressed,false);
  assert.equal(visibility.get(id).rawDefaultHidden,true);
  assert.equal(visibility.get(id).suppressed,true);
  assert.equal(visibility.get(id).available,false);
  assert.ok(snapshot.resolved.every(n=>![...n.directPartIds,...n.descendantPartIds].includes(id)));
 }
 const donorSource=snapshot.reviewedConceptBindings.find(c=>c.conceptId==='FMA8644');
 assert.equal(donorSource.term,'A12.2.01.219');
 assert.deepEqual(donorSource.partIds,['FJ2884','FJ2885']);
 assert.ok(snapshot.rawConceptBindings.some(c=>c.conceptId==='FMA8644'));
 assert.equal(snapshot.conceptBindings.some(c=>c.conceptId==='FMA8644'),false,'Retained donor bindings became effective installed bindings');
 const nodes=new Map(snapshot.resolved.map(n=>[n.id,n]));
 assert.ok(nodes.get('A12.2.01.219').directPartIds.every(id=>!['FJ2884','FJ2885'].includes(id)));
 assert.deepEqual(nodes.get('A12.2.04.007').directPartIds,['TA98BODY:male:carotid-body:left:macro']);
 assert.deepEqual(new Set(nodes.get('A12.2.04.009').directPartIds),new Set(['TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall','TA98ARTERIAL:male:carotid-bifurcation:right:finite-wall']));
 const host='TA98PART:ZA:Internal carotid artery.l:A12.2.06.002';
 for(const code of ['A12.2.04.008','A12.2.06.003'])assert.deepEqual(nodes.get(code).directPartIds,[host]);
 const effective=ordinaryCirculationCatalogue(atlas,api).effective;
 const registry=api.carotidFeatureRegistry(effective);assert.ok(registry);
 assert.equal(registry.termGroups.length,1);assert.equal(registry.termGroups[0].referenceOnly,true);
 assert.equal(registry.termGroups[0].representation,'surface-region');
 assert.equal(api.hierarchyEntries(effective).filter(e=>e.id==='TA98REFERENCE:male:carotid-sinus:left:medial-proximal-ICA').length,1);
 assert.ok(snapshot.resolved.every(n=>[...n.directPartIds,...n.descendantPartIds].every(id=>visibility.get(id)?.available)));
 assert.equal(hash(JSON.stringify(atlas)),before);assert.equal(hash(fs.readFileSync(source)),hash(raw));
});
