import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {build} from 'esbuild';

const built=await build({entryPoints:['app/cardiac-label-sites.ts'],bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const temporary=fs.mkdtempSync(path.join(os.tmpdir(),'atlas-heart-compaction-'));
const stage='.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1';
const read=p=>JSON.parse(fs.readFileSync(p));
const child=spawnSync('.local-models/geometry-tools-env/bin/python',['-B','-c',`
import json,runpy,sys
from pathlib import Path
compact=runpy.run_path('server/build-separate-spinal-model.py')['compact_catalogue']
for sex,source in [('male','${stage}/atlas.json'),('female','.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json')]:
 a=json.load(open(source)); active={p['chunk'] for p in a['parts'] if not p.get('suppressed')}
 # Whole-body delivery removes unavailable/detail-only bank references.
 a['parts']=[p for p in a['parts'] if p['chunk'] in active]
 compact(a);Path(sys.argv[1],sex+'.json').write_text(json.dumps(a))
`,temporary],{encoding:'utf8',timeout:10000});
assert.equal(child.status,0,child.stderr);

try{
 for(const sex of ['male','female'])test(`${sex}: production compaction preserves 113 heart controls and rejects stale mapping or pose`,()=>{
  const raw=read(path.join(temporary,sex+'.json')),bundle=read(sex==='male'?`${stage}/male/cardiac-label-sites.json`:`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/cardiac-label-sites.json`);
  const atlas=api.withCardiacLabelBundle(raw,bundle);assert.ok(api.cardiacLabelRegistry(atlas));assert.equal(api.cardiacLabeledTermParts(atlas).size,113);
  const carrier=Object.keys(bundle.integration.carrierPoses)[0],chunk=raw.parts.find(p=>p.id===carrier).chunk;
  assert.equal(raw.sourceChunkIndices[chunk],bundle.integration.carrierPoses[carrier].chunk);
  for(const mutate of [
   a=>a.sourceChunkIndices.pop(),
   a=>a.sourceChunkIndices[chunk]=-1,
   a=>a.sourceChunkIndices[chunk]=.5,
   a=>a.sourceChunkIndices[chunk]=a.sourceChunkCount,
   a=>a.sourceChunkIndices[chunk]=a.sourceChunkIndices[chunk===0?1:0],
   a=>a.sourceChunkCount=a.chunks.length-1,
   a=>{const other=chunk===0?1:0;[a.sourceChunkIndices[chunk],a.sourceChunkIndices[other]]=[a.sourceChunkIndices[other],a.sourceChunkIndices[chunk]];},
   a=>a.parts.find(p=>p.id===carrier).geometryOrigin[0]+=.001,
   a=>a.parts.find(p=>p.id===carrier).positions+=4,
   a=>{delete a.sourceChunkIndices;delete a.sourceChunkCount;},
  ]){const wrong=structuredClone(raw);mutate(wrong);assert.equal(api.withCardiacLabelBundle(wrong,bundle),wrong);}
 });
 for(const sex of ['male','female'])test(`${sex}: refreshed ordinary body loads its actual installed 113 controls`,t=>{
  const raw=read(`.local-models/organ-models/spinal-cord/${sex}-body.json`);
  if(!raw.sourceChunkIndices){t.skip('Ordinary body has not yet been refreshed with source chunk mapping');return;}
  const bundle=read(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/cardiac-label-sites.json`);
  const loaded=api.withCardiacLabelBundle(raw,bundle);
  assert.ok(api.cardiacLabelRegistry(loaded));assert.equal(api.cardiacLabeledTermParts(loaded).size,113);
 });
}finally{
 // Test callbacks retain their decoded values, so remove only after the suite.
 process.on('exit',()=>fs.rmSync(temporary,{recursive:true,force:true}));
}
