import {useCallback,useEffect,useRef,useState} from 'react';
import {assetUrl} from './asset-url';
import {parseGuestHierarchy,type GuestHierarchy} from './guest-hierarchy';

const builtIn=assetUrl('/assets/chakras.json');
const canonical=(input:string)=>{const url=new URL(input,location.href);if(!['http:','https:'].includes(url.protocol))throw new Error('Use an HTTP or HTTPS hierarchy URL.');return url.href;};

export function useGuestHierarchies(){
 const [loaded,setLoaded]=useState<{url:string;hierarchy:GuestHierarchy}[]>([]),[sources,setSources]=useState<string[]>([]),[error,setError]=useState('');
 const loadedRef=useRef(loaded);loadedRef.current=loaded;
 const load=useCallback(async(input:string)=>{
  const url=canonical(input),existing=loadedRef.current.find(item=>item.url===url);if(existing)return existing.hierarchy;
  const response=await fetch(url);if(!response.ok)throw new Error(`Hierarchy URL returned ${response.status}.`);
  const length=Number(response.headers.get('content-length')??0);if(length>500_000)throw new Error('Hierarchy file is too large.');
  const content=await response.text();if(content.length>500_000)throw new Error('Hierarchy file is too large.');
  let raw:unknown;try{raw=JSON.parse(content);}catch{throw new Error('Hierarchy URL must provide JSON.');}
  const hierarchy=parseGuestHierarchy(raw);
  if(loadedRef.current.some(item=>item.hierarchy.id===hierarchy.id&&item.url!==url))throw new Error(`A hierarchy named ${hierarchy.id} is already loaded.`);
  setLoaded(current=>current.some(item=>item.url===url)?current:[...current,{url,hierarchy}]);
  if(url!==canonical(builtIn))setSources(current=>current.includes(url)?current:[...current,url]);
  setError('');return hierarchy;
 },[]);
 useEffect(()=>{const inputs=[builtIn,...new URLSearchParams(location.search).getAll('guest')];void Promise.allSettled(inputs.map(input=>load(input).catch(cause=>{setError(cause instanceof Error?cause.message:'Could not load a guest hierarchy.');throw cause;})));},[load]);
 const restore=useCallback((search:string)=>{for(const input of new URLSearchParams(search).getAll('guest'))void load(input).catch(cause=>setError(cause instanceof Error?cause.message:'Could not load a guest hierarchy.'));},[load]);
 return {hierarchies:loaded.map(item=>item.hierarchy),sources,error,load,restore};
}
