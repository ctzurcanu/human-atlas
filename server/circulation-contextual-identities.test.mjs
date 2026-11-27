import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {existsSync, readFileSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';

const source='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';

test('male vascular and node homonyms resolve by explicit source and authoritative region', {skip:!existsSync(source)}, async()=>{
 const raw=readFileSync(source),sha=value=>createHash('sha256').update(value).digest('hex');
 const catalogue=JSON.parse(raw),metadata=JSON.parse(readFileSync('app/data/ta98-metadata.json'));
 const bundled=await build({stdin:{contents:"export * from './app/guest-hierarchy';export * from './app/reviewed-catalogue-identities';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
 const api=await import('data:text/javascript;base64,'+Buffer.from(bundled.outputFiles[0].contents).toString('base64'));
 const effective=api.withReviewedCatalogueIdentities(catalogue);
 const roots=api.resolveGuestHierarchy(effective,api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json'))));
 const nodes=new Map(),walk=node=>{if(nodes.has(node.id))return;nodes.set(node.id,node);node.children.forEach(walk);};roots.forEach(walk);
 const cases=[
  ['ZA:Dorsal digital arteries of hand.l','A12.2.09.035','A12.2.16.053'],
  ['ZA:Dorsal digital arteries of foot.r','A12.2.16.053','A12.2.09.035'],
  ['ZA:Superficial anterior cervical nodes','A13.3.00.018','A13.3.00.027'],
  ['ZA:Superficial lateral cervical nodes.l','A13.3.00.027','A13.3.00.018'],
  ['ZA:Paratracheal cervical nodes','A13.3.00.024','A13.3.02.011'],
  ['ZA:Paratracheal thoracic nodes','A13.3.02.011','A13.3.00.024'],
  ['ZA:Superior pancreatic nodes','A13.3.03.024','A13.3.03.028'],
  ['ZA:Superior pancreaticoduodenal nodes','A13.3.03.028','A13.3.03.024'],
  ['ZA:Medial common iliac nodes','A13.3.04.004','A13.3.04.010'],
  ['ZA:Medial external iliac nodes.r','A13.3.04.010','A13.3.04.004'],
  ['ZA:Inferior superficial inguinal nodes.l','A13.3.05.006','A13.3.05.013'],
  ['ZA:(Proximal deep inguinal node).l','A13.3.05.008','A13.3.05.009'],
  ['ZA:(Intermediate deep inguinal node).r','A13.3.05.009','A13.3.05.008'],
  ['ZA:Deep popliteal nodes.r','A13.3.05.013','A13.3.05.007'],
  ['ZA:Medial segmental artery of right lung','A12.2.01.111','A12.2.01.118'],
  ['ZA:Medial basal segmental artery of right lung','A12.2.01.118','A12.2.01.111'],
  ['ZA:Medial vein of right lung','A12.3.02.113','A12.3.02.112'],
  ['ZA:Lateral vein of right lung','A12.3.02.112','A12.3.02.113'],
  ['ZA:Anterior vein of right lung','A12.3.02.105','A12.3.02.305'],
  ['ZA:Anterior vein of left lung','A12.3.02.305','A12.3.02.105'],
  ['ZA:Apical vein of right lung','A12.3.02.102','A12.3.02.202'],
  ['ZA:Superior vein of right lung','A12.3.02.202','A12.3.02.102'],
  ['ZA:Superior basal vein of right lung','A12.3.02.206','A12.3.02.406'],
  ['ZA:Superior basal vein of left lung','A12.3.02.406','A12.3.02.206'],
  ['ZA:Inferior basal vein of right lung','A12.3.02.210','A12.3.02.410'],
  ['ZA:Inferior basal vein of left lung','A12.3.02.410','A12.3.02.210'],
  ['ZA:Posterior lateral nasal branches of sphenopalatine artery..l','A12.2.05.089','A12.2.05.090'],
  ['ZA:Posterior lateral nasal branches of sphenopalatine artery..r','A12.2.05.089','A12.2.05.090'],
  ['ZA:Insular branches of middle cerebral artery (M2).r','A12.2.07.053','A12.2.07.054'],
  ['ZA:Insular branches of middle cerebral artery (M2-segment).l','A12.2.07.053','A12.2.07.054'],
  ['ZA:Right superior phrenic vein','A12.3.07.010','A12.3.09.002'],
  ['ZA:Parieto-occipital artery.l','A12.2.07.101','A12.2.07.045'],
  ['ZA:Parieto-occipital artery.r','A12.2.07.101','A12.2.07.045'],
  ['ZA:Calcaneal branches of fibular artery.l','A12.2.16.075','A12.2.16.059'],
  ['ZA:Calcaneal branches of fibular artery.r','A12.2.16.075','A12.2.16.059'],
 ];
 for(const [id,term,other] of cases){
  const concept=effective.concepts.find(c=>c.id===id);assert.ok(concept,id);
  assert.equal(concept.ta98Term,term,id);assert.equal(concept.ta98Kind,'exact',id);
  assert.equal(metadata.byConceptMatch[id]?.term,term,id);
  for(const partId of concept.elements){
   assert.ok(nodes.get(term)?.directParts.some(p=>p.id===partId),`${id} directly selectable at ${term}`);
   assert.ok(!nodes.get(other)?.directParts.some(p=>p.id===partId),`${id} must not resolve at homonymous ${other}`);
  }
 }
 // Reviewed rules must not let a generic name or duplicate supplement inherit a specific regional identity.
 const fake={id:'UNREGISTERED:generic-lateral-node',name:'Lateral nodes',elements:[],ta98Kind:undefined};
 const probe=api.withReviewedCatalogueIdentities({...catalogue,concepts:[...catalogue.concepts,fake]});
 assert.equal(probe.concepts.at(-1).ta98Term,undefined);
 for(const id of ['O3M:Dorsal digital arteries.l','O3M:Dorsal digital arteries.r']){
  assert.equal(metadata.reviewedCatalogueIdentities?.[id],undefined,id);
  assert.equal(effective.concepts.find(c=>c.id===id)?.ta98Term,undefined,id);
 }
 assert.equal(sha(readFileSync(source)),sha(raw),'identity resolution cannot mutate installed source');
 assert.deepEqual(catalogue,JSON.parse(raw),'identity projection cannot mutate source arrays or records');
});
