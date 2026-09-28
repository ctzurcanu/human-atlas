import assert from 'node:assert/strict';
import {test} from 'node:test';
import {build} from 'esbuild';
import {IDBFactory,IDBObjectStore} from 'fake-indexeddb';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';
import {mkdtemp,mkdir,writeFile,rm} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {assetManifest} from '../server/asset-manifest.mjs';

const bundle=await build({stdin:{contents:"export * from './app/asset-cache'; export * from './app/model-download';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {createAssetCache,initializeAssetCache,fetchAsset,fetchAssetJson,loadModelBuffer,decodeModelResponse}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);
const base='http://atlas.test/human-atlas/';
function fixture(initial={'models/body.bin':'body-v1'}){
 const files=new Map(Object.entries(initial)),calls=new Map(),values=new Map(),factory=new IDBFactory();
 const storage={getItem:key=>values.get(key)??null,setItem:(key,value)=>values.set(key,value),removeItem:key=>values.delete(key)};
 const fetcher=async(input,init={})=>{
  if(init.signal?.aborted)throw init.signal.reason;
  const key=new URL(input,base).pathname.slice('/human-atlas/'.length);calls.set(key,(calls.get(key)??0)+1);
  if(key==='asset-manifest.json')return Response.json({schema:1,files:Object.fromEntries([...files].map(([key,body])=>[key,{version:createHash('sha256').update(body).digest('hex'),bytes:Buffer.byteLength(body)}]))});
  return files.has(key)?new Response(files.get(key)):new Response('missing',{status:404});
 };
 const options={baseUrl:base,indexedDB:factory,storage,fetcher};
 return {files,calls,values,factory,storage,fetcher,options,cache:(extra={})=>createAssetCache({...options,...extra})};
}
async function browser(f,run,href=base){
 const originals=Object.fromEntries(['indexedDB','localStorage','fetch','location','history'].map(key=>[key,Object.getOwnPropertyDescriptor(globalThis,key)]));
 Object.assign(globalThis,{indexedDB:f.factory,localStorage:f.storage,fetch:f.fetcher,location:{href},history:{state:{preserved:1},replaceState(state,_title,url){this.state=state;globalThis.location.href=url;}}});
 try{await run();}finally{for(const [key,descriptor] of Object.entries(originals)){if(descriptor)Object.defineProperty(globalThis,key,descriptor);else delete globalThis[key];}}
}

test('persistent model, texture and hierarchy hits survive a new viewer instance',async()=>{
 const f=fixture({'models/body.bin':'geometry','models/tissue.png':'texture','assets/cell-types.json':'{"nodes":[]}'});
 let cache=f.cache();await cache.initialize();
 for(const key of f.files.keys())assert.equal(await (await cache.load(base+key)).text(),f.files.get(key));
 cache=f.cache();await cache.initialize();
 for(const key of f.files.keys()){assert.equal(await (await cache.load(base+key)).text(),f.files.get(key));assert.equal(f.calls.get(key),1);}
 assert.equal(f.calls.get('asset-manifest.json'),2);
 assert.ok(f.values.has('human-atlas-assets-v1:/human-atlas/:index'));
});

test('changed files refresh, unchanged files keep their cache',async()=>{
 const f=fixture({'models/body.bin':'version1','assets/genes.json':'same'}),first=f.cache();await first.initialize();
 await first.load(base+'models/body.bin');await first.load(base+'assets/genes.json');
 f.files.set('models/body.bin','version2');const second=f.cache();await second.initialize();
 assert.equal(await (await second.load(base+'models/body.bin')).text(),'version2');await second.load(base+'assets/genes.json');
 assert.equal(f.calls.get('models/body.bin'),2);assert.equal(f.calls.get('assets/genes.json'),1);
});

test('rc=1 clears only atlas cache, consumes the flag and preserves view, invitation and fragment',async()=>{
 const f=fixture();f.values.set('atlas-saved-views','saved');f.values.set('atlas-slides','slides');
 await browser(f,async()=>{
  await initializeAssetCache('/human-atlas/');await fetchAsset(base+'models/body.bin');
  const href=base+'?model=cell&connect=invitation&rc=1&select=one&select=two#atlas-view=1.packed';location.href=href;
  await initializeAssetCache('/human-atlas/');
  assert.equal(location.href,href.replace('&rc=1',''));assert.deepEqual(history.state,{preserved:1});
  assert.equal(f.values.get('atlas-saved-views'),'saved');assert.equal(f.values.get('atlas-slides'),'slides');
  await fetchAsset(base+'models/body.bin');assert.equal(f.calls.get('models/body.bin'),2);
 });
});

test('reset does not let an older in-flight download repopulate the cache',async()=>{
 const f=fixture();let release,started;const begun=new Promise(resolve=>started=resolve),waiting=new Promise(resolve=>release=resolve);let delay=true;
 const cache=f.cache({fetcher:async(input,init)=>{if(String(input).endsWith('body.bin')&&delay){started();await waiting;}return f.fetcher(input,init);}});
 await cache.initialize();const pending=cache.load(base+'models/body.bin');await begun;
 await cache.reset();await cache.initialize();delay=false;release();await pending;
 await cache.load(base+'models/body.bin');assert.equal(f.calls.get('models/body.bin'),2);
});

test('LRU eviction uses metadata and respects its size budget',async()=>{
 const f=fixture({'models/a.bin':'a'.repeat(30),'models/b.bin':'b'.repeat(30),'models/c.bin':'c'.repeat(30)}),cache=f.cache({maxBytes:70});await cache.initialize();
 await cache.load(base+'models/a.bin');await cache.load(base+'models/b.bin');
 await new Promise(resolve=>setTimeout(resolve,2));await cache.load(base+'models/a.bin');await cache.load(base+'models/c.bin');
 await cache.load(base+'models/a.bin');await cache.load(base+'models/c.bin');assert.equal(f.calls.get('models/a.bin'),1);assert.equal(f.calls.get('models/c.bin'),1);
 await cache.load(base+'models/b.bin');assert.equal(f.calls.get('models/b.bin'),2);
});

test('blocked storage and quota failures fall back to downloads',async()=>{
 const f=fixture(),blocked=f.cache({indexedDB:{open(){throw new DOMException('Blocked','SecurityError');}},storage:{getItem(){throw Error('blocked');},setItem(){throw Error('blocked');},removeItem(){throw Error('blocked');}}});
 await blocked.initialize();assert.equal(await (await blocked.load(base+'models/body.bin')).text(),'body-v1');await blocked.reset();
 const cache=f.cache();await cache.initialize();const put=IDBObjectStore.prototype.put;
 IDBObjectStore.prototype.put=function(...args){if(this.name==='files')throw new DOMException('Full','QuotaExceededError');return put.apply(this,args);};
 try{for(let i=0;i<2;i++)assert.equal(await (await cache.load(base+'models/body.bin')).text(),'body-v1');}
 finally{IDBObjectStore.prototype.put=put;}
});

test('cached and network loads respect aborts',async()=>{
 const f=fixture(),cache=f.cache();await cache.initialize();await cache.load(base+'models/body.bin');
 const controller=new AbortController();const pending=cache.load(base+'models/body.bin',{signal:controller.signal});controller.abort();
 await assert.rejects(pending,{name:'AbortError'});await assert.rejects(cache.load(base+'models/body.bin',{signal:controller.signal}),{name:'AbortError'});
 assert.equal(f.calls.get('models/body.bin'),1);
 const networkController=new AbortController(),network=f.cache({fetcher:async(input,init)=>{const response=await f.fetcher(input,init);if(String(input).endsWith('body.bin'))networkController.abort();return response;}});await network.initialize();
 await assert.rejects(network.load(base+'models/body.bin',{signal:networkController.signal,cache:'reload'}),{name:'AbortError'});
});

test('cached compressed and HTTP-decoded downloads both decode once',async()=>{
 const data=Buffer.from('anatomy vertices and normals'),compressed=gzipSync(data),f=fixture({'models/body.bin.gz':compressed}),cache=f.cache();await cache.initialize();
 for(let i=0;i<2;i++)assert.deepEqual(Buffer.from(await decodeModelResponse(await cache.load(base+'models/body.bin.gz'),data.length,true)),data);
 assert.equal(f.calls.get('models/body.bin.gz'),1);
 const decoded=f.cache({fetcher:async(input,init)=>String(input).endsWith('.gz')?new Response(data,{headers:{'content-encoding':'gzip'}}):f.fetcher(input,init)});
 await decoded.initialize();await decoded.load(base+'models/body.bin.gz',{cache:'reload'});
 const hit=await decoded.load(base+'models/body.bin.gz');assert.equal(hit.headers.has('content-encoding'),false);assert.deepEqual(Buffer.from(await decodeModelResponse(hit,data.length,true)),data);
});

test('damaged cached binary and JSON recover with one fresh download',async()=>{
 const f=fixture({'models/body.bin':'valid-data','models/body.json':'{"valid":true}'});
 await browser(f,async()=>{
  await initializeAssetCache('/human-atlas/');await fetchAsset(base+'models/body.bin');await fetchAsset(base+'models/body.json');
  const db=await new Promise((resolve,reject)=>{const r=f.factory.open('human-atlas-assets-v1:/human-atlas/',1);r.onsuccess=()=>resolve(r.result);r.onerror=()=>reject(r.error);});
  for(const key of ['models/body.bin','models/body.json'])await new Promise((resolve,reject)=>{
   const tx=db.transaction('files','readwrite'),store=tx.objectStore('files'),get=store.get(base+key);
   get.onsuccess=()=>store.put({...get.result,body:new Blob(['!'])});tx.oncomplete=resolve;tx.onabort=()=>reject(tx.error);
  });
  assert.equal(Buffer.from(await loadModelBuffer(base+'models/body.bin',10,false,new AbortController().signal)).toString(),'valid-data');
  assert.deepEqual(await fetchAssetJson(base+'models/body.json'),{valid:true});
  assert.equal(f.calls.get('models/body.bin'),2);assert.equal(f.calls.get('models/body.json'),2);db.close();
 });
});

test('live APIs, external files, errors and query-bearing URLs bypass the cache',async()=>{
 const f=fixture(),cache=f.cache();await cache.initialize();
 for(let i=0;i<2;i++)for(const path of ['atlas-connect/identity','models/body.bin?live=1','missing.bin'])await cache.load(base+path);
 assert.equal(f.calls.get('atlas-connect/identity'),2);assert.equal(f.calls.get('models/body.bin'),2);assert.equal(f.calls.get('missing.bin'),2);
 let requests=0;const external=f.cache({fetcher:async(input,init)=>String(input).startsWith('http://external.test')?(requests++,new Response('external')):f.fetcher(input,init)});await external.initialize();
 await external.load('http://external.test/models/body.bin');await external.load('http://external.test/models/body.bin');assert.equal(requests,2);
});

test('cached index and files remain usable if the manifest cannot be fetched',async()=>{
 const f=fixture(),first=f.cache();await first.initialize();await first.load(base+'models/body.bin');
 const offline=f.cache({fetcher:async()=>{throw Error('offline');}});await offline.initialize();assert.equal(await (await offline.load(base+'models/body.bin')).text(),'body-v1');
});

test('build hashes change only with asset contents, including same-length changes',async()=>{
 const directory=await mkdtemp(join(tmpdir(),'atlas-cache-manifest-'));
 try{
  await mkdir(join(directory,'models'));await writeFile(join(directory,'models','body.bin'),'version1');
  const first=await assetManifest(directory),again=await assetManifest(directory);assert.deepEqual(first,again);
  await writeFile(join(directory,'models','body.bin'),'version2');const changed=await assetManifest(directory);
  assert.notEqual(first.files['models/body.bin'].version,changed.files['models/body.bin'].version);assert.equal(changed.files['models/body.bin'].bytes,8);
 }finally{await rm(directory,{recursive:true,force:true});}
});
