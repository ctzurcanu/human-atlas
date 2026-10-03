import test from 'node:test';
import assert from 'node:assert/strict';
import {gzipSync} from 'node:zlib';
import {build} from 'esbuild';
const compiled=await build({entryPoints:['app/model-download.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {decodeModelResponse}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
test('gzip model buffers load with native and fallback decoding, without double decoding',async()=>{
 const bytes=Buffer.from('original anatomical vertices and indices'),gzip=gzipSync(bytes),native=globalThis.DecompressionStream;
 for(const decoder of [native,undefined]){try{globalThis.DecompressionStream=decoder;
  assert.deepEqual(Buffer.from(await decodeModelResponse(new Response(gzip),bytes.length,true)),bytes);
  assert.deepEqual(Buffer.from(await decodeModelResponse(new Response(bytes),bytes.length,true)),bytes);
  await assert.rejects(decodeModelResponse(new Response(gzip),bytes.length+1,true),/incomplete/);
 }finally{globalThis.DecompressionStream=native;}}
});
