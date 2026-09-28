import {useCallback,useEffect,useRef,useState} from 'react';
import {assetUrl} from './asset-url';
import {parseGuestHierarchy,type GuestHierarchy} from './guest-hierarchy';
export const BUILT_IN_GUESTS=['genes','cell-types','physiology','dermatomes-myotomes','drugs','physical-exercise','chakras'] as const;
const names=['Genes','Cell Types','Physiology','Dermatomes and Myotomes','Drugs','Physical Exercise','Chakras'];
const builtIns=BUILT_IN_GUESTS.map(id=>assetUrl(`/assets/${id}.json`));
const canonical=(input:string)=>{const url=new URL(input,location.href);if(!['http:','https:'].includes(url.protocol))throw new Error('Use an HTTP or HTTPS hierarchy URL.');return url.href;};
export function useGuestHierarchies(){
 const [loaded,setLoaded]=useState<{url:string;hierarchy:GuestHierarchy}[]>([]),[sources,setSources]=useState<string[]>([]),[error,setError]=useState('');
 const loadedRef=useRef(loaded);loadedRef.current=loaded;
 const pending=useRef(new Map<string,Promise<GuestHierarchy>>());
 const load=useCallback((input:string):Promise<GuestHierarchy>=>{
  const url=canonical(input),existing=loadedRef.current.find(item=>item.url===url);
  if(existing)return Promise.resolve(existing.hierarchy);
  const inFlight=pending.current.get(url);if(inFlight)return inFlight;
  const builtIn=builtIns.some(input=>canonical(input)===url),limit=builtIn?20_000_000:500_000;
  const promise=(async()=>{
   const response=await fetch(url);if(!response.ok)throw new Error(`Hierarchy URL returned ${response.status}.`);
   const length=Number(response.headers.get('content-length')??0);if(length>limit)throw new Error('Hierarchy file is too large.');
   const content=await response.text();if(content.length>limit)throw new Error('Hierarchy file is too large.');
   let raw:unknown;try{raw=JSON.parse(content);}catch{throw new Error('Hierarchy URL must provide JSON.');}
   const hierarchy=parseGuestHierarchy(raw);
   if(!builtIn&&(BUILT_IN_GUESTS as readonly string[]).includes(hierarchy.id))throw new Error('This hierarchy ID is reserved for a built-in hierarchy.');
   if(loadedRef.current.some(item=>item.hierarchy.id===hierarchy.id&&item.url!==url))throw new Error(`A hierarchy named ${hierarchy.id} is already loaded.`);
   setLoaded(current=>current.some(item=>item.url===url)?current:[...current,{url,hierarchy}]);
   if(!builtIn)setSources(current=>current.includes(url)?current:[...current,url]);
   setError('');return hierarchy;
  })();pending.current.set(url,promise);void promise.finally(()=>pending.current.delete(url)).catch(()=>{});return promise;
 },[]);
 const extend=useCallback(async(id:string,nodeId:string)=>{
  const entry=loadedRef.current.find(item=>item.hierarchy.id===id),node=entry?.hierarchy.nodes.find(node=>node.id===nodeId);
  if(!node?.extension)return;
  const url=assetUrl(node.extension);const key=canonical(url);const active=pending.current.get(key);if(active){await active;return;}
  const promise=(async()=>{
   const response=await fetch(url);if(!response.ok)throw new Error('Could not load transcript branches.');
   const text=await response.text();if(text.length>5_000_000)throw new Error('Transcript file is too large.');
   const extension=parseGuestHierarchy(JSON.parse(text));if(extension.schema!=='human-atlas-hierarchy/v2'||extension.roots?.length!==1||extension.roots[0]!==nodeId)throw new Error('Invalid transcript extension.');
   // Add a chunk to the latest graph, preserving other concurrently opened genes.
   setLoaded(current=>current.map(item=>{if(item.hierarchy.id!==id)return item;const added=new Map(extension.nodes.map(node=>[node.id!,node]));const existing=new Set(item.hierarchy.nodes.map(node=>node.id));
    const nodes=item.hierarchy.nodes.map(node=>added.get(node.id!)??node).concat(extension.nodes.filter(node=>!existing.has(node.id)));
    return {...item,hierarchy:{...parseGuestHierarchy({...item.hierarchy,nodes}),revision:(item.hierarchy.revision??0)+1}};
   }));return extension;
  })();pending.current.set(key,promise);try{await promise;}catch(cause){setError(cause instanceof Error?cause.message:'Could not load transcripts.');throw cause;}finally{pending.current.delete(key);}
 },[]);
 const ensure=useCallback(async(id:string)=>{const index=(BUILT_IN_GUESTS as readonly string[]).indexOf(id);if(index<0)return;try{return await load(builtIns[index]);}catch(cause){setError(cause instanceof Error?cause.message:'Could not load hierarchy.');throw cause;}},[load]);
 useEffect(()=>{
  const params=new URLSearchParams(location.search),id=params.get('tree')?.replace(/^guest:/,'');
  const inputs=[builtIns[BUILT_IN_GUESTS.indexOf('chakras')],...(id&&(BUILT_IN_GUESTS as readonly string[]).includes(id)?[builtIns[(BUILT_IN_GUESTS as readonly string[]).indexOf(id)]]:[]),...params.getAll('guest')];
  void Promise.allSettled(inputs.map(input=>load(input).catch(cause=>{setError(cause instanceof Error?cause.message:'Could not load a guest hierarchy.');throw cause;})));
 },[load]);
 const restore=useCallback((search:string)=>{
  const id=new URLSearchParams(search).get('tree')?.replace(/^guest:/,'');if(id)void ensure(id).catch(()=>{});
  for(const input of new URLSearchParams(search).getAll('guest'))void load(input).catch(cause=>setError(cause instanceof Error?cause.message:'Could not load a guest hierarchy.'));
 },[load,ensure]);
 const rank=(id:string)=>{const index=(BUILT_IN_GUESTS as readonly string[]).indexOf(id);return index<0?BUILT_IN_GUESTS.length:index;};
 const hierarchies=loaded.map(item=>item.hierarchy).sort((a,b)=>rank(a.id)-rank(b.id));
 const menuHierarchies:GuestHierarchy[]=[...BUILT_IN_GUESTS.map((id,index)=>hierarchies.find(item=>item.id===id)??{schema:'human-atlas-hierarchy/v1' as const,id,name:names[index],nodes:[]}),...hierarchies.filter(item=>!(BUILT_IN_GUESTS as readonly string[]).includes(item.id))];
 return {hierarchies,menuHierarchies,sources,error,load,restore,ensure,extend};
}
