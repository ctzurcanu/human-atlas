import assert from 'node:assert/strict';
import {existsSync,readFileSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';
const base='.local-models/geometry-workbench/biodigital-male-circulation-v1/';
const profile=base+'atlas-cardiovascular-review-v15/atlas.json',local=existsSync(profile);
let api,atlas,old;
if(local){
 const out=await build({stdin:{contents:"export * from './app/biodigital-heart';export * from './app/cardiac-label-sites';export * from './app/ta98-term-bindings';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
 api=await import('data:text/javascript;base64,'+Buffer.from(out.outputFiles[0].contents).toString('base64'));
 atlas=JSON.parse(readFileSync(profile));old=JSON.parse(readFileSync(base+'atlas-cardiovascular-review-v14/atlas.json'));
}
const svc='BDM:male:human_20_male_cardiovascular_system-superior_vena_cava_ID',ra='BDH:male:human_20_male_cardiovascular_system-right_atrium_ID';
test('installed cardiac labels and the untouched v14 rollback both validate',{skip:!local},()=>{
 assert.equal(api.sourceCardiacIntegration(atlas)?.termGroups.length,46);
 assert.equal(api.sourceCardiacIntegration(old)?.termGroups.length,43);
 assert.equal(api.cardiacComponentChoices(atlas,atlas.concepts.find(c=>c.id==='BDM:male:human_20_male_cardiovascular_system-heart_ID')).length,70);
 assert.equal([...api.catalogueTa98TermParts(atlas)].filter(([,ps])=>ps.some(p=>p.id.startsWith('BDM:male:'))).length,679);
 assert.deepEqual(atlas.sourceCardiacIntegration.carrierPoses,old.sourceCardiacIntegration.carrierPoses);
 const prior=new Map(old.parts.map(p=>[p.id,p]));
 assert.equal(atlas.parts.length,old.parts.length);
 for(const p of atlas.parts){
  if(p.id===svc)continue;
  if(p.id===ra){const {tissueDrawGroups,...record}=p;assert.deepEqual(record,prior.get(p.id));}
  else assert.deepEqual(p,prior.get(p.id));
 }
});
test('the revised SVC reference rejects stale sites, altered caps and forged current carriers',{skip:!local},()=>{
 const site=a=>a.sourceCardiacIntegration.referenceSites.find(s=>s.term==='A12.1.01.011');
 for(const change of [
  a=>{a.parts.find(p=>p.id===svc).geometryOrigin[0]+=.001;},
  a=>{a.parts.find(p=>p.id===svc).tissueDrawGroups.find(g=>g.neighbor).neighbor='other';},
  a=>{a.parts.find(p=>p.id===ra).tissueDrawGroups[0].count-=3;},
  a=>{const p=a.parts.find(p=>p.id===svc);a.chunks[p.chunk].sha256='f'.repeat(64);},
  a=>{site(a).geometry.paths[0].points[0][0]+=.001;},
  a=>{site(a).geometry.filled=true;},
  a=>{a.sourceCardiacIntegration.referenceSites=a.sourceCardiacIntegration.referenceSites.map(s=>s.term==='A12.1.01.011'?structuredClone(site(old)):s);},
 ]){const a=structuredClone(atlas);change(a);assert.equal(api.sourceCardiacIntegration(a),null);}
});
