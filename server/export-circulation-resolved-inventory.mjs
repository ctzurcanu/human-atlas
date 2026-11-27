// Read the same identity/binding/hierarchy APIs as the viewer; no geometry edits.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
import {circulationCatalogueSnapshot} from './circulation-effective-catalogue.mjs';
const root=process.cwd(),sha=b=>createHash('sha256').update(b).digest('hex');
const arguments_=process.argv.slice(2),sex=arguments_.shift()??'male';assert.ok(['male','female'].includes(sex));
let output,sourceOption,expectedSource;
while(arguments_.length){
 const option=arguments_.shift();
 if(option==='--output'){assert.ok(!output,'--output may be supplied once');output=arguments_.shift();assert.ok(output,'--output needs a new JSON path');}
 else if(option==='--source'){assert.ok(!sourceOption,'--source may be supplied once');sourceOption=arguments_.shift();assert.ok(sourceOption,'--source needs a catalogue path');}
 else if(option==='--expected-source-sha256'){assert.ok(!expectedSource,'--expected-source-sha256 may be supplied once');expectedSource=arguments_.shift();assert.match(expectedSource??'',/^[a-f0-9]{64}$/,'An explicit catalogue SHA256 is required');}
 else throw new Error(`Unknown option: ${option}`);
}
assert.ok(!sourceOption||expectedSource,'A nondefault catalogue requires --expected-source-sha256');
const source=path.resolve(root,sourceOption??`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`),raw=fs.readFileSync(source);
assert.equal(sha(raw),expectedSource??(sex==='male'?'1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39':'8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'),'Explicit expected catalogue revision changed');
const atlas=JSON.parse(raw),before=sha(JSON.stringify(atlas)),metadata=JSON.parse(fs.readFileSync('app/data/ta98-metadata.json'));
const compiled=await build({stdin:{contents:`export * from './app/guest-hierarchy';export * from './app/reviewed-catalogue-identities';export * from './app/preferred-source-representations';`,resolveDir:root},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const definition=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
const snapshot=circulationCatalogueSnapshot(atlas,metadata,api,definition);
assert.equal(sha(JSON.stringify(atlas)),before);assert.equal(sha(fs.readFileSync(source)),sha(raw));
const out=output?path.resolve(root,output):path.join(root,`.local-models/geometry-workbench/circulation-census-v1/resolved-${sex}.json`);assert.ok(!fs.existsSync(out));fs.mkdirSync(path.dirname(out),{recursive:true});
fs.writeFileSync(out,JSON.stringify({...snapshot,sex,source,sourceSha256:sha(raw),hierarchySha256:sha(fs.readFileSync('public/assets/ta98.json')),metadataSha256:sha(fs.readFileSync('app/data/ta98-metadata.json')),anatomicalAcceptance:false},null,2)+'\n',{flag:'wx'});
console.log(JSON.stringify({output:out,sourceSha256:sha(raw),conceptBindings:snapshot.conceptBindings.length,resolvedNodes:snapshot.resolved.length,availableParts:snapshot.effectivePartVisibility.filter(p=>p.available).length,sourceUnchanged:true}));
