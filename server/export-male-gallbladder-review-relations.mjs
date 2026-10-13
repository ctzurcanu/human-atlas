import {build} from 'esbuild';
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
const path='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json',raw=readFileSync(path),atlas=JSON.parse(raw);
const compiled=await build({stdin:{contents:"export * from './app/gallbladder-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
// The complete original source supplies organ-level references for review.
// No retired record is reactivated in the real catalogue or clinical junction.
const parts=atlas.parts.map(p=>p.id==='ZA:Gallbladder'?{...p,suppressed:false}:p);
const references=(api.gallbladderRelations(parts).get('ZA:Gallbladder')??[]).map(r=>({kind:r.kind,id:r.target.id,name:r.target.name,via:r.via??[],note:r.note}));
if(!references.some(r=>r.kind==='arterial')||!references.some(r=>r.kind==='venous')||!references.some(r=>r.kind==='innervation'))throw new Error('Incomplete native organ reference routes');
writeFileSync('.local-models/ta98-review/male-gallbladder-native-outlet/functional-references.json',JSON.stringify({sourceCatalogueSha256:createHash('sha256').update(raw).digest('hex'),references,wholeOrganReferencesOnly:true,physicalJunctionsVerified:false},null,2)+'\n');
console.log(`${references.length} source-qualified whole-gallbladder reference routes; no catalogue or geometry changes.`);
