import {useState} from 'react';
import {BookmarkPlus,Copy,ExternalLink,Trash2,X} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {DEFAULT_EMBED_UI,EMBED_UI,embedUrl,iframeCode,type EmbedUi} from './embed';

interface Bookmark {id:string;name:string;url:string;createdAt:number}
const storageKey='human-atlas-bookmarks-v1';
function readBookmarks():Bookmark[]{
 try{
  const value=JSON.parse(localStorage.getItem(storageKey)??'[]');
  if(!Array.isArray(value))return [];
  return value.filter(item=>{
   if(typeof item?.id!=='string'||typeof item?.name!=='string'||typeof item?.url!=='string')return false;
   try{return new URL(item.url).origin===location.origin;}catch{return false;}
  });
 }catch{return [];}
}

interface Props {viewUrl:string;defaultName:string;close:()=>void}
export default function ShareView({viewUrl,defaultName,close}:Props){
 const [ui,setUi]=useState<EmbedUi[]>(DEFAULT_EMBED_UI),[height,setHeight]=useState(600);
 const [bookmarks,setBookmarks]=useState(readBookmarks),[name,setName]=useState(defaultName),[status,setStatus]=useState('');
 const code=iframeCode(viewUrl,ui,height),preview=embedUrl(viewUrl,ui);
 const copy=async(value:string,label:string)=>{try{await navigator.clipboard.writeText(value);setStatus(`${label} copied`);}catch{setStatus('Select the text above to copy it.');}};
 const toggle=(id:EmbedUi)=>setUi(current=>current.includes(id)?current.filter(item=>item!==id):EMBED_UI.map(item=>item.id).filter(item=>item===id||current.includes(item)));
 const persist=(next:Bookmark[])=>{try{localStorage.setItem(storageKey,JSON.stringify(next));setBookmarks(next);return true;}catch{setStatus('This browser could not save views locally.');return false;}};
 const save=()=>{
  const title=name.trim();if(!title){setStatus('Give this view a name first.');return;}
  const next=[{id:crypto.randomUUID(),name:title,url:viewUrl,createdAt:Date.now()},...bookmarks].slice(0,100);
  if(persist(next)){setName('');setStatus('View saved in this browser.');}
 };
 return <div className="share-view glass" role="dialog" aria-label="Save or open views">
  <div className="panel-heading"><span>Save or open views</span><Button variant="ghost" className="icon-button" onClick={close} aria-label="Close saved views"><X size={18}/></Button></div>
  <section className="bookmark-section" aria-label="Save current view"><h3>Save current view</h3><p>Your model, camera angle, zoom, selection, and layers are saved in this browser.</p>
   <div className="bookmark-add"><input aria-label="Name this view" placeholder="Name this view" value={name} onChange={e=>setName(e.target.value)} onKeyDown={e=>{if(e.key==='Enter')save();}}/><Button variant="outline" onClick={save}><BookmarkPlus size={15}/> Save view</Button></div>
   {bookmarks.length>0&&<><h4>Saved views</h4><ul>{bookmarks.map(item=><li key={item.id}><Button variant="ghost" onClick={()=>location.assign(item.url)} title={item.url}>{item.name}</Button><Button variant="ghost" className="icon-button" onClick={()=>{if(persist(bookmarks.filter(bookmark=>bookmark.id!==item.id)))setStatus('Saved view removed.');}} aria-label={`Remove ${item.name}`}><Trash2 size={15}/></Button></li>)}</ul></>}
  </section>
  <p role="status">{status}</p>
  <h3 className="share-subheading">Share this view</h3>
  <label>View URL<input aria-label="View URL" readOnly value={viewUrl} onFocus={e=>e.target.select()}/></label>
  <Button variant="outline" onClick={()=>void copy(viewUrl,'Link')}><Copy size={15}/> Copy link</Button>
  <fieldset className="embed-options"><legend>Controls in iframe</legend>
   <div className="embed-option-actions"><Button variant="ghost" onClick={()=>setUi([])}>Viewer only</Button><Button variant="ghost" onClick={()=>setUi(EMBED_UI.map(item=>item.id))}>Full controls</Button></div>
   <div className="embed-option-grid">{EMBED_UI.map(item=><label key={item.id}><input type="checkbox" checked={ui.includes(item.id)} onChange={()=>toggle(item.id)}/>{item.label}</label>)}</div>
  </fieldset>
  <label>Iframe height (px)<input aria-label="Iframe height" type="number" min={340} max={1200} step={10} value={height} onChange={e=>setHeight(Number(e.target.value))}/></label>
  <label>Embed on a website<textarea aria-label="Iframe embed code" readOnly rows={4} value={code} onFocus={e=>e.target.select()}/></label>
  <div className="embed-actions"><Button variant="outline" onClick={()=>void copy(code,'Embed code')}><Copy size={15}/> Copy embed code</Button><a href={preview} target="_blank" rel="noopener noreferrer"><ExternalLink size={15}/> Preview iframe</a></div>
 </div>;
}
