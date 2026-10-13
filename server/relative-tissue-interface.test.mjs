import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export {relativeTissueInterface} from './app/relative-tissue-interface';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {relativeTissueInterface}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
for(const sex of ['female','male'])test(`${sex}: browser scalar arithmetic reproduces the independently checked packed interface exactly`,async()=>{
 const bank='.local-models/ta98-review/colon-integrated-wall/',candidate='.local-models/ta98-review/colon-relative-muscle-thickness/';
 const m=JSON.parse(await readFile(bank+sex+'-manifest.json','utf8')),c=JSON.parse(await readFile(candidate+sex+'-manifest.json','utf8'));
 const data=await readFile(bank+sex+'-wall.bin'),raw=data.buffer.slice(data.byteOffset,data.byteOffset+data.byteLength);
 const bytes=await readFile(candidate+sex+'-weights.bin'),weights=new Float32Array(bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength));
 const a=new Float32Array(raw,m.interfaces+2*m.vertices*12,m.vertices*3),b=new Float32Array(raw,m.interfaces+3*m.vertices*12,m.vertices*3),before=Buffer.from(a.buffer,a.byteOffset,a.byteLength).slice();
 const result=relativeTissueInterface(a,b,weights);
 assert.equal(createHash('sha256').update(Buffer.from(result.buffer)).digest('hex'),c.derivedInterfaceSha256);
 assert.deepEqual(Buffer.from(a.buffer,a.byteOffset,a.byteLength),before,'original shared bank stays immutable');
 assert.notEqual(result.buffer,a.buffer);
 for(let i=0;i<weights.length;i++)if(weights[i]===0)assert.deepEqual(result.slice(i*3,i*3+3),a.slice(i*3,i*3+3));
});
test('Invalid scalar counts, extrapolation and nonfinite controls are rejected',()=>{
 assert.throws(()=>relativeTissueInterface(new Float32Array(3),new Float32Array(3),new Float32Array(2)));
 for(const value of [NaN,Infinity,-.1,.6])assert.throws(()=>relativeTissueInterface(new Float32Array(3),new Float32Array(3),Float32Array.of(value)));
});
