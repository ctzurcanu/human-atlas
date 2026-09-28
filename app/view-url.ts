import {gzipSync,gunzipSync,strToU8,strFromU8} from 'fflate';

const REQUEST_LIMIT=4096;
const STATE_LIMIT=2_000_000;
const PREFIX='#atlas-view=1.';
let cachedHash='',cachedQuery='';

function fragmentQuery(hash:string){
 if(!hash.startsWith(PREFIX))return '';
 if(hash===cachedHash)return cachedQuery;
 try{
  const encoded=hash.slice(PREFIX.length);
  if(encoded.length>STATE_LIMIT*2||!/^[\w-]+$/.test(encoded))return '';
  const binary=atob(encoded.replaceAll('-','+').replaceAll('_','/'));
  const bytes=Uint8Array.from(binary,char=>char.charCodeAt(0));
  if(bytes.length<18||bytes[0]!==31||bytes[1]!==139)return '';
  const size=new DataView(bytes.buffer).getUint32(bytes.length-4,true);
  if(size>STATE_LIMIT)return '';
  cachedQuery=strFromU8(gunzipSync(bytes,new Uint8Array(size)));cachedHash=hash;
  return cachedQuery;
 }catch{return '';}
}

// Accept both legacy query strings and complete URLs, including saved slides.
export function viewParameters(input:string){
 const value=/^https?:\/\//.test(input)?new URL(input):null;
 const separator=input.indexOf('#');
 const search=value?.search??(separator<0?input:input.slice(0,separator));
 const hash=value?.hash??(separator<0?'':input.slice(separator));
 const params=new URLSearchParams(fragmentQuery(hash));
 const query=new URLSearchParams(search);
 // Query controls can override a saved view, e.g. iframe controls or hierarchy.
 for(const key of new Set(query.keys())){
  const values=query.getAll(key),saved=params.getAll(key);
  if(values.length===saved.length&&values.every((value,index)=>value===saved[index]))continue;
  params.delete(key);for(const value of values)params.append(key,value);
 }
 return params;
}

export function compactViewUrl(input:string){
 const url=new URL(input),params=viewParameters(input);
 url.search=params.toString();url.hash='';
 if(url.href.length<=REQUEST_LIMIT)return url.href;
 const data=strToU8(params.toString());
 if(data.length>STATE_LIMIT)throw new Error('This view is too large to save as a link.');
 const compressed=gzipSync(data,{level:6,mtime:0});
 let binary='';for(const byte of compressed)binary+=String.fromCharCode(byte);
 url.hash=PREFIX+btoa(binary).replaceAll('+','-').replaceAll('/','_').replace(/=+$/,'');
 url.search='';
 // Keep routing controls readable to existing integrations. The full state is
 // self-contained in the fragment; it never reaches a static hosting server.
 for(const key of ['model','tree','embed','ui','connect','connectJoin','relay']){
  const value=params.get(key);if(value!==null&&value.length<=512)url.searchParams.set(key,value);
 }
 return url.href;
}
