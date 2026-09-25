import {useEffect,useRef,useState} from 'react';
import {ChevronDown,Link2} from 'lucide-react';
import type {GuestHierarchy} from './guest-hierarchy';

type Props={hierarchies:GuestHierarchy[];active:string;onSelect:(id:string)=>void;onLoad:(url:string)=>Promise<void>};

export default function GuestHierarchyMenu({hierarchies,active,onSelect,onLoad}:Props){
 const [open,setOpen]=useState(false),[url,setUrl]=useState(''),[loading,setLoading]=useState(false),[error,setError]=useState('');
 const host=useRef<HTMLDivElement>(null);
 useEffect(()=>{if(!open)return;const outside=(event:PointerEvent)=>{if(!host.current?.contains(event.target as Node))setOpen(false);};const escape=(event:KeyboardEvent)=>{if(event.key==='Escape')setOpen(false);};document.addEventListener('pointerdown',outside);document.addEventListener('keydown',escape);return()=>{document.removeEventListener('pointerdown',outside);document.removeEventListener('keydown',escape);};},[open]);
 const submit=async(event:React.FormEvent)=>{event.preventDefault();if(!url.trim())return;setLoading(true);setError('');try{await onLoad(url.trim());setUrl('');setOpen(false);}catch(cause){setError(cause instanceof Error?cause.message:'Could not load this hierarchy.');}finally{setLoading(false);}};
 const selected=hierarchies.find(item=>active===`guest:${item.id}`);
 return <div className="guest-menu-wrap" ref={host}>
  <button type="button" className={selected?'active':''} aria-label="Guest hierarchies" title={selected?.name??'Guest hierarchies'} aria-haspopup="menu" aria-expanded={open} onClick={()=>setOpen(value=>!value)}>Guest <ChevronDown size={12}/></button>
  {open&&<div className="guest-menu glass" role="menu" aria-label="Guest hierarchies">
   <div className="guest-menu-title">Guest hierarchies</div>
   {hierarchies.map(item=><button type="button" role="menuitem" className="guest-menu-item" key={item.id} onClick={()=>{onSelect(item.id);setOpen(false);}}>{item.name}</button>)}
   {!hierarchies.length&&<span className="guest-menu-empty">Loading hierarchies…</span>}
   <form onSubmit={submit} className="guest-menu-form">
    <label htmlFor="guest-hierarchy-url"><Link2 size={13}/> Load from URL</label>
    <input id="guest-hierarchy-url" type="url" value={url} onChange={event=>setUrl(event.target.value)} placeholder="https://…/hierarchy.json" required/>
    <button type="submit" disabled={loading}>{loading?'Loading…':'Load hierarchy'}</button>
    {error&&<p role="alert">{error}</p>}
   </form>
  </div>}
 </div>;
}
