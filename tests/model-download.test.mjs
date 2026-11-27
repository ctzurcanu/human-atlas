import test from 'node:test';
import assert from 'node:assert/strict';
import {gzipSync} from 'node:zlib';
import {build} from 'esbuild';
const compiled=await build({entryPoints:['app/model-download.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {decodeModelResponse,loadModelBuffer}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
test('gzip model buffers load with native and fallback decoding, without double decoding',async()=>{
 const bytes=Buffer.from('original anatomical vertices and indices'),gzip=gzipSync(bytes),native=globalThis.DecompressionStream;
 for(const decoder of [native,undefined]){try{globalThis.DecompressionStream=decoder;
  assert.deepEqual(Buffer.from(await decodeModelResponse(new Response(gzip),bytes.length,true)),bytes);
  assert.deepEqual(Buffer.from(await decodeModelResponse(new Response(bytes),bytes.length,true)),bytes);
  await assert.rejects(decodeModelResponse(new Response(gzip),bytes.length+1,true),/incomplete/);
 }finally{globalThis.DecompressionStream=native;}}
});
test('a transient heart geometry fetch retries once and returns the complete buffer',async()=>{
 const original=globalThis.fetch,bytes=Buffer.from('heart geometry'),requests=[];
 try{
  globalThis.fetch=async(url,init)=>{requests.push({url,init});if(requests.length===1)throw new TypeError('Failed to fetch');return new Response(bytes);};
  assert.deepEqual(Buffer.from(await loadModelBuffer('/local-models/heart/cells.bin',bytes.length,false,new AbortController().signal)),bytes);
  assert.equal(requests.length,2);assert.equal(requests[1].init.cache,'reload');
 }finally{globalThis.fetch=original;}
});
test('a heart geometry body interrupted after headers also retries once',async()=>{
 const original=globalThis.fetch,bytes=Buffer.from('female heart geometry');let requests=0;
 try{
  globalThis.fetch=async()=>++requests===1?new Response(new ReadableStream({start(controller){controller.error(new TypeError('Failed to fetch'));}})):new Response(bytes);
  assert.deepEqual(Buffer.from(await loadModelBuffer('/local-models/female-heart/cells.bin',bytes.length,false,new AbortController().signal)),bytes);
  assert.equal(requests,2);
 }finally{globalThis.fetch=original;}
});
test('persistent fetch failure stops after two attempts and identifies the failing asset',async()=>{
 const original=globalThis.fetch;let requests=0;
 try{
  globalThis.fetch=async()=>{requests++;throw new TypeError('Failed to fetch');};
  await assert.rejects(loadModelBuffer('/local-models/heart/cells.bin',16,false,new AbortController().signal),error=>error.message.includes('/local-models/heart/cells.bin')&&error.cause instanceof TypeError);
  assert.equal(requests,2);
 }finally{globalThis.fetch=original;}
});
test('cancelled anatomy fetches and missing files do not retry',async()=>{
 const original=globalThis.fetch,controller=new AbortController();let requests=0;
 try{
  globalThis.fetch=async()=>{requests++;controller.abort();throw controller.signal.reason;};
  await assert.rejects(loadModelBuffer('/cancelled.bin',16,false,controller.signal),{name:'AbortError'});assert.equal(requests,1);
  requests=0;globalThis.fetch=async()=>{requests++;return new Response('',{status:404});};
  await assert.rejects(loadModelBuffer('/missing.bin',16,false,new AbortController().signal),/missing\.bin/);assert.equal(requests,1);
 }finally{globalThis.fetch=original;}
});
