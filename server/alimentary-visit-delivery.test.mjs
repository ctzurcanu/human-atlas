import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import crypto from 'node:crypto';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/alimentary-visit';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {matchesAlimentaryVisit}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
for(const sex of ['male','female'])test(`${sex}: digestive visit remains available after oesophageal and rectal wall installation`,()=>{
 const raw=fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`),atlas=JSON.parse(raw),route=JSON.parse(fs.readFileSync(`.local-models/ta98-review/visits/${sex}.json`));
 assert.equal(route.sourceCatalogueSha256,crypto.createHash('sha256').update(raw).digest('hex'));
 assert.ok(matchesAlimentaryVisit(atlas,route));
 for(const label of ['Esophagus','Rectum']){const s=route.segments.find(s=>s.name===label);assert.ok(s?.parts.length);assert.ok(s.points.length>20);assert.ok(s.parts.every(id=>{const p=atlas.parts.find(p=>p.id===id);return p&&!p.suppressed&&p.pairedTissue&&(p.tissueWallLayer==='epithelium'||p.caecocolicLayer==='epithelium');}));}
 const binding=route.bindings[0],part=atlas.parts.find(p=>p.id===binding.id);
 assert.equal(matchesAlimentaryVisit({...atlas,parts:atlas.parts.map(p=>p===part?{...p,suppressed:true}:p)},route),false);
 assert.equal(matchesAlimentaryVisit({...atlas,parts:atlas.parts.map(p=>p===part?{...p,indices:p.indices+4}:p)},route),false);
});
