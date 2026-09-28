/** Large public files live in IndexedDB, never in the small localStorage quota. */
export interface AssetManifest {schema:1;files:Record<string,{version:string;bytes:number}>}
interface CachedFile {url:string;version:string;body:Blob;headers:[string,string][]}
interface Usage {url:string;version:string;size:number;used:number}
interface CacheOptions {
 baseUrl:string;
 indexedDB?:IDBFactory;
 storage?:Pick<Storage,'getItem'|'setItem'|'removeItem'>;
 fetcher?:typeof fetch;
 maxBytes?:number;
 manifestTimeout?:number;
}
const LIMIT=512*1024*1024;
const cachedResponses=new WeakSet<Response>();
const abort=(signal?:AbortSignal|null)=>{if(signal?.aborted)throw signal.reason??new DOMException('Aborted','AbortError');};
const request=<T>(operation:IDBRequest<T>)=>new Promise<T>((resolve,reject)=>{operation.onsuccess=()=>resolve(operation.result);operation.onerror=()=>reject(operation.error);});
const completion=(transaction:IDBTransaction)=>{const done=new Promise<void>((resolve,reject)=>{transaction.oncomplete=()=>resolve();transaction.onabort=transaction.onerror=()=>reject(transaction.error??new Error('Cache transaction failed.'));});void done.catch(()=>{});return done;};
function validManifest(value:unknown):value is AssetManifest {
 if(!value||typeof value!=='object')return false;
 const manifest=value as AssetManifest;
 return manifest.schema===1&&!!manifest.files&&typeof manifest.files==='object'&&Object.entries(manifest.files).every(([path,file])=>/^(models|assets|local-models)\//.test(path)&&file&&typeof file.version==='string'&&Number.isFinite(file.bytes)&&file.bytes>=0);
}
function responseFrom(file:CachedFile){
 const headers=new Headers(file.headers);
 // Fetch may already have decoded an HTTP Content-Encoding. The stored Blob is the delivered body.
 headers.delete('content-encoding');headers.set('content-length',String(file.body.size));
 return new Response(file.body,{headers});
}

export function createAssetCache(options:CacheOptions){
 const base=new URL(options.baseUrl),fetcher=options.fetcher??fetch,limit=options.maxBytes??LIMIT;
 const databaseName=`human-atlas-assets-v1:${base.pathname}`,indexKey=databaseName+':index';
 let database:Promise<IDBDatabase|null>|undefined,manifest:AssetManifest|undefined,generation=0,writes=true;
 const keyFor=(input:string)=>{const url=new URL(input,base);return url.origin===base.origin&&url.pathname.startsWith(base.pathname)&&!url.search?decodeURIComponent(url.pathname.slice(base.pathname.length)):'';};
 const indexRead=()=>{try{const value=JSON.parse(options.storage?.getItem(indexKey)??'null');return validManifest(value)?value:undefined;}catch{return undefined;}};
 const open=()=>database??=(async()=>{
  if(!options.indexedDB)return null;
  return new Promise<IDBDatabase|null>(resolve=>{
   let settled=false;
   const finish=(value:IDBDatabase|null)=>{if(settled){value?.close();return;}settled=true;clearTimeout(timer);resolve(value);};
   const timer=setTimeout(()=>finish(null),1500);
   try{
    const operation=options.indexedDB!.open(databaseName,1);
    operation.onupgradeneeded=()=>{const db=operation.result;db.createObjectStore('files',{keyPath:'url'});db.createObjectStore('usage',{keyPath:'url'});};
    operation.onsuccess=()=>{const db=operation.result;db.onversionchange=()=>{db.close();database=undefined;};finish(db);};
    operation.onerror=operation.onblocked=()=>finish(null);
   }catch{finish(null);}
  });
 })();
 const remove=async(url:string)=>{
  try{const db=await open();if(!db)return;const tx=db.transaction(['files','usage'],'readwrite'),done=completion(tx);tx.objectStore('files').delete(url);tx.objectStore('usage').delete(url);await done;}catch{/* Caching is optional. */}
 };
 const prune=async()=>{
  try{
   const db=await open();if(!db||!manifest)return;
   const tx=db.transaction(['files','usage'],'readwrite'),done=completion(tx);
   const entries=await request(tx.objectStore('usage').getAll()) as Usage[];
   for(const item of entries)if(manifest.files[keyFor(item.url)]?.version!==item.version){tx.objectStore('files').delete(item.url);tx.objectStore('usage').delete(item.url);}
   await done;
  }catch{/* Keep loading normally when storage is unavailable. */}
 };
 const read=async(url:string,version:string)=>{
  try{
   const db=await open();if(!db)return;
   const tx=db.transaction('files','readonly'),done=completion(tx);
   const file=await request(tx.objectStore('files').get(url)) as CachedFile|undefined;await done;
   if(!file)return;
   if(file.version!==version||!(file.body instanceof Blob)||!file.body.size){await remove(url);return;}
   // Only small metadata is read when choosing files to evict.
   const touch=db.transaction('usage','readwrite'),touched=completion(touch);
   touch.objectStore('usage').put({url,version,size:file.body.size,used:Date.now()});void touched.catch(()=>{});
   return file;
  }catch{return undefined;}
 };
 const write=async(file:CachedFile,started:number)=>{
  if(!writes||file.body.size>limit||started!==generation)return;
  try{
   const db=await open();if(!db||started!==generation)return;
   const save=async(reclaim=false)=>{
    const tx=db.transaction(['files','usage'],'readwrite'),done=completion(tx);
    // Attach a handler immediately: a quota abort can arrive before a request promise settles.
    void done.catch(()=>{});
    const usage=tx.objectStore('usage'),entries=(await request(usage.getAll()) as Usage[]).filter(item=>item.url!==file.url).sort((a,b)=>a.used-b.used);
    let total=entries.reduce((sum,item)=>sum+item.size,0)+file.body.size;
    const target=reclaim?Math.min(limit,total/2):limit;
    for(const entry of entries){if(total<=target)break;tx.objectStore('files').delete(entry.url);usage.delete(entry.url);total-=entry.size;}
    if(started!==generation){tx.abort();await done;return;}
    tx.objectStore('files').put(file);usage.put({url:file.url,version:file.version,size:file.body.size,used:Date.now()});await done;
   };
   try{await save();}catch(error){if(error instanceof DOMException&&error.name==='QuotaExceededError')await save(true);else throw error;}
  }catch{if(started===generation)writes=false;}
 };
 const initialize=async()=>{
  manifest=indexRead();
  const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),options.manifestTimeout??4000);
  try{
   const response=await fetcher(new URL('asset-manifest.json',base),{cache:'no-cache',signal:controller.signal});
   if(response.ok){const next:unknown=await response.json();if(validManifest(next)){manifest=next;try{options.storage?.setItem(indexKey,JSON.stringify(next));}catch{/* IndexedDB works even if localStorage is blocked. */}}}
  }catch{/* A previously cached index still permits reading cached models offline. */}
  finally{clearTimeout(timer);}
  // Old file versions are dropped before loading a new catalogue.
  await prune();
 };
 const reset=async()=>{
  generation++;writes=true;manifest=undefined;
  try{options.storage?.removeItem(indexKey);}catch{/* No other browser preferences are removed. */}
  try{const db=await open();if(!db)return;const tx=db.transaction(['files','usage'],'readwrite'),done=completion(tx);tx.objectStore('files').clear();tx.objectStore('usage').clear();await done;}catch{/* Storage may be disabled by the browser. */}
 };
 const load=async(input:string,init:RequestInit={})=>{
  abort(init.signal);
  const url=new URL(input,base).href,entry=manifest?.files[keyFor(url)],started=generation;
  // Live APIs, external hierarchies and query-bearing requests keep their normal network semantics.
  if(!entry||init.method&&init.method!=='GET')return fetcher(input,init);
  if(init.cache!=='reload'&&init.cache!=='no-store'){
   const file=await read(url,entry.version);abort(init.signal);
   if(file&&started===generation){const response=responseFrom(file);cachedResponses.add(response);return response;}
  }
  const response=await fetcher(input,init);abort(init.signal);
  if(!response.ok||response.status!==200||init.cache==='no-store')return response;
  const body=await response.blob();abort(init.signal);
  const file:CachedFile={url,version:entry.version,body,headers:[...response.headers]};
  // Persist compressed downloads rather than decoded geometry arrays, avoiding base64 and JSON copies.
  if(body.size)await write(file,started);abort(init.signal);
  return responseFrom(file);
 };
 return {initialize,reset,load,remove:(input:string)=>remove(new URL(input,base).href)};
}

let active:ReturnType<typeof createAssetCache>|undefined;
export async function initializeAssetCache(basePath:string){
 let storage:Storage|undefined,factory:IDBFactory|undefined;
 try{storage=localStorage;}catch{/* Blocked browser storage. */}
 try{factory=indexedDB;}catch{/* Blocked browser storage. */}
 active=createAssetCache({baseUrl:new URL(basePath,location.href).href,indexedDB:factory,storage});
 const url=new URL(location.href);
 if(url.searchParams.get('rc')==='1'){
  await active.reset();url.searchParams.delete('rc');history.replaceState(history.state,'',url.href);
 }
 await active.initialize();
}
export const fetchAsset=(url:string,init?:RequestInit)=>active?active.load(url,init):fetch(url,init);
export const forgetAsset=(url:string)=>active?.remove(url)??Promise.resolve();
export const isCachedAsset=(response:Response)=>cachedResponses.has(response);

export async function fetchAssetJson<T>(url:string,init:RequestInit={},message='The anatomy catalogue could not be loaded.'):Promise<T>{
 const response=await fetchAsset(url,init);
 if(!response.ok)throw new Error(message);
 try{const value=await response.json() as T;abort(init.signal);return value;}
 catch(error){
  if(init.signal?.aborted)throw error;
  await forgetAsset(url);if(!isCachedAsset(response))throw error;
  const fresh=await fetchAsset(url,{...init,cache:'reload'});if(!fresh.ok)throw new Error(message);
  try{const value=await fresh.json() as T;abort(init.signal);return value;}catch(error){if(!init.signal?.aborted)await forgetAsset(url);throw error;}
 }
}
