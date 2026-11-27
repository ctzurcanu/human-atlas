import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';

const result=await build({stdin:{contents:"export * from './app/reviewed-catalogue-identities';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {withReviewedCatalogueIdentities}=await import('data:text/javascript;base64,'+Buffer.from(result.outputFiles[0].contents).toString('base64'));
const metadata=JSON.parse(readFileSync('app/data/ta98-metadata.json'));
const catalogue=JSON.parse(readFileSync('.local-models/geometry-workbench/biodigital-male-circulation-v1/hand-current-palm-fit-v981/ordinary-atlas.json'));

test('both installed imported hand veins lose the false exact foot identity without losing meshes',()=>{
 assert.match(metadata.byCode['A12.3.11.014'].latin,/pedis/);
 assert.match(metadata.byCode['A12.3.08.023'].latin,/manus/);
 const current=withReviewedCatalogueIdentities(catalogue);
 for(const side of ['left','right']){
  const id=`BDM:male:human_20_male_cardiovascular_system-${side}_dorsal_digital_veins_ID`;
  const part=current.parts.find(p=>p.id===id),before=catalogue.parts.find(p=>p.id===id);
  assert.equal(part.ta98Term,undefined);
  assert.deepEqual({...part,ta98Term:before.ta98Term},before);
  const c=current.concepts.find(c=>c.id===id);
  assert.equal(c.ta98Term,'A12.3.08.023');assert.equal(c.ta98Kind,'parent');
  assert.ok(!current.ta98TermBindings.some(b=>b.conceptId===id&&b.term==='A12.3.11.014'));
 }
 assert.equal(current.parts.length,catalogue.parts.length);
 assert.deepEqual(current.chunks,catalogue.chunks);
 const sourceFoot=catalogue.parts.filter(p=>p.ta98Term==='A12.3.11.014'&&!p.id.includes('dorsal_digital_veins_ID'));
 for(const p of sourceFoot)assert.deepEqual(current.parts.find(q=>q.id===p.id),p);
});
