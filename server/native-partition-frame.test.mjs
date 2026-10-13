import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({entryPoints:['app/ta98-modeled-parent.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {nativePartFrame}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const atlas=JSON.parse(readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
const partitions=atlas.parts.filter(p=>!p.suppressed&&p.id.startsWith('TA98PART:ZA:'));
test('actual native TA98 partitions retain the source frame after subdivision',()=>{
 assert.ok(partitions.length>=292);
 for(const p of partitions){assert.equal(nativePartFrame(p),'ZA:',p.id);assert.equal(nativePartFrame(p),'ZA:');}
 for(const p of atlas.parts.filter(p=>p.id.startsWith('ZA:')))assert.equal(nativePartFrame(p),'ZA:');
});
test('partition frame cannot be inferred from terminology, donor provenance or a conflicting source ID',()=>{
 const p=partitions.find(p=>p.sourceId==='Gallbladder');assert.ok(p);
 for(const bad of [{...p,id:'HRA:VH_F_gallbladder'}, {...p,id:'TA98PART:FJ3080:neck'}, {...p,sourceId:'Other'}, {...p,sourceId:undefined}, {...p,provenance:undefined}, {...p,provenance:{...p.provenance,url:'https://example.org/donor'}}, {...p,license:'unverified'}])assert.equal(nativePartFrame(bad),undefined);
 assert.equal(nativePartFrame({...p,id:'LOCAL:female:gallbladder:1030'}),'LOCAL:female:');
 assert.equal(nativePartFrame({...p,id:'LOCAL:male:gallbladder:1030'}),'LOCAL:male:');
});
