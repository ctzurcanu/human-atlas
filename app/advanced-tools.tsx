import {useEffect,useRef,useState,type PointerEvent as ReactPointerEvent} from 'react';
import {GripVertical,Play,Plus,SkipBack,SkipForward,Square,StepBack,StepForward,Trash2,X} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {Dialog,DialogContent,DialogDescription,DialogFooter,DialogHeader,DialogTitle} from '@/components/ui/dialog';
import type {SceneState} from './anatomy';
import ConnectTools from './connect-tools';
import type {AtlasConnection} from './use-atlas-connection';
import type {AdvancedPresentation} from './connect-state';

type ToolType='slides'|'quiz'|'layers'|'connect';
type CapturedScene={model:string;hierarchy:string;guestSources:string[];state:SceneState};
type SavedView={id:string;name:string;url:string;model:string;addedAt:string;scene?:CapturedScene};
export type ToolDocument={schemaVersion:1;type:ToolType;views:SavedView[]};
type ViewDrag={id:string;pointerId:number;startY:number;insertion:number};
const STORAGE_KEY='human-atlas-advanced-tools-v1';
const emptyDocument:ToolDocument={schemaVersion:1,type:'slides',views:[]};
const validType=(value:unknown):value is ToolType=>value==='slides'||value==='quiz'||value==='layers'||value==='connect';
function validateDocument(value:unknown):value is ToolDocument{
 if(!value||typeof value!=='object')return false;
 const document=value as Partial<ToolDocument>;
 return document.schemaVersion===1&&validType(document.type)&&Array.isArray(document.views)&&document.views.every(view=>{
  if(!view||typeof view!=='object'||typeof view.id!=='string'||typeof view.name!=='string'||typeof view.url!=='string'||typeof view.model!=='string'||typeof view.addedAt!=='string')return false;
  try{return new URL(view.url).origin===location.origin;}catch{return false;}
 });
}
function readDocument():ToolDocument{
 try{const saved=localStorage.getItem(STORAGE_KEY);if(saved){const parsed:unknown=JSON.parse(saved);if(validateDocument(parsed))return parsed;}}catch{/* Storage may be unavailable in an embedded viewer. */}
 return emptyDocument;
}

interface Props{close:()=>void;currentViewUrl:()=>string;captureScene:()=>CapturedScene;openView:(url:string,animate:boolean)=>void;model:string;viewName:string;ready:boolean;connection:AtlasConnection;onPresentationChange:(value:AdvancedPresentation)=>void;remotePresentation?:AdvancedPresentation|null}
export default function AdvancedTools({close,currentViewUrl,captureScene,openView:openSavedView,model,viewName,ready,connection,onPresentationChange,remotePresentation}:Props){
 const [document,setDocument]=useState<ToolDocument>(()=>new URLSearchParams(location.search).has('connect')?{...readDocument(),type:'connect'}:readDocument());
 const [editorOpen,setEditorOpen]=useState(false);
 const [editorText,setEditorText]=useState('');
 const [editorError,setEditorError]=useState('');
 const [newViewName,setNewViewName]=useState('');
 const [editingView,setEditingView]=useState<{id:string;name:string}|null>(null);
 const [slideIndex,setSlideIndex]=useState<number|null>(null);
 const slidePlaying=slideIndex!==null&&!!document.views[slideIndex];
 const presenting=slidePlaying||document.type==='connect'&&connection.active&&(connection.role==='guest'||!connection.controlsOpen);
 const presentationCallback=useRef(onPresentationChange);presentationCallback.current=onPresentationChange;
 useEffect(()=>presentationCallback.current({document,slideIndex}),[document,slideIndex]);
 useEffect(()=>{if(connection.role==='host'&&connection.controlsOpen){setDocument(current=>({...current,type:'connect'}));setSlideIndex(null);}},[connection.controlsOpen,connection.role]);
 useEffect(()=>{
  if(!remotePresentation)return;
  setDocument({...remotePresentation.document,views:remotePresentation.document.views.map(view=>{const url=new URL(view.url);return {...view,url:new URL(url.pathname+url.search,location.origin).href};})});
  setSlideIndex(remotePresentation.slideIndex);
 },[remotePresentation]);
 const [draggingId,setDraggingId]=useState<string|null>(null);
 const [dropIndex,setDropIndex]=useState<number|null>(null);
 const dragRef=useRef<ViewDrag|null>(null);
 const listRef=useRef<HTMLDivElement>(null);
 const viewClickRef=useRef<ReturnType<typeof setTimeout>|null>(null);
 const cancelViewClick=()=>{if(viewClickRef.current!==null){clearTimeout(viewClickRef.current);viewClickRef.current=null;}};
 useEffect(()=>()=>{if(viewClickRef.current!==null)clearTimeout(viewClickRef.current);},[]);
 const update=(next:ToolDocument)=>{setDocument(next);try{localStorage.setItem(STORAGE_KEY,JSON.stringify(next));}catch{/* The current session still retains the document. */}};
 const addView=()=>{
  const next:SavedView={id:crypto.randomUUID(),name:newViewName.trim()||viewName||`View ${document.views.length+1}`,url:currentViewUrl(),model,addedAt:new Date().toISOString(),scene:captureScene()};
  update({...document,views:[...document.views,next]});
  setNewViewName('');
 };
 const removeView=(id:string)=>update({...document,views:document.views.filter(view=>view.id!==id)});
 const openEditor=()=>{setEditorText(JSON.stringify(document,null,2));setEditorError('');setEditorOpen(true);};
 const applyEditor=()=>{
  try{
   const parsed:unknown=JSON.parse(editorText);
   if(!validateDocument(parsed))throw new Error('Use schemaVersion 1, a valid type, and views with same-site URLs.');
   update(parsed);setEditorOpen(false);
  }catch(error){setEditorError(error instanceof Error?error.message:'Invalid JSON.');}
 };
 const openView=(view:SavedView,animate=false)=>openSavedView(view.url,animate);
 const clickView=(view:SavedView)=>{
  cancelViewClick();
  viewClickRef.current=setTimeout(()=>{viewClickRef.current=null;openView(view);},300);
 };
 const editView=(view:SavedView)=>{cancelViewClick();setEditingView({id:view.id,name:view.name});};
 const saveViewName=()=>{
  if(!editingView)return;
  const name=editingView.name.trim();
  if(name)update({...document,views:document.views.map(view=>view.id===editingView.id?{...view,name}:view)});
  setEditingView(null);
 };
 const showSlide=(index:number)=>{
  const count=document.views.length;if(!count)return;
  const next=(index%count+count)%count,view=document.views[next];
  const sameModel=new URL(view.url).searchParams.get('model')===new URL(location.href).searchParams.get('model');
  setSlideIndex(next);openView(view,sameModel&&!window.matchMedia('(prefers-reduced-motion: reduce)').matches);
 };
 const stopSlides=()=>setSlideIndex(null);
 useEffect(()=>{
  if(slideIndex===null)return;
  const onKey=(event:KeyboardEvent)=>{
   if(event.defaultPrevented||event.altKey||event.ctrlKey||event.metaKey)return;
   const target=event.target;
   if(target instanceof HTMLElement&&(target.isContentEditable||target.closest('input,textarea,select,[role="textbox"]')))return;
   const next=event.code==='Space'||event.key==='ArrowRight'||event.key==='>';
   const previous=event.key==='ArrowLeft'||event.key==='<';
   const stop=event.key==='ArrowDown';
   if(!next&&!previous&&!stop)return;
   event.preventDefault();
   if(event.repeat)return;
   if(stop)stopSlides();else showSlide(slideIndex+(next?1:-1));
  };
  window.addEventListener('keydown',onKey);
  return()=>window.removeEventListener('keydown',onKey);
 },[slideIndex,document.views]);
 const moveView=(id:string,target:number)=>{
  const from=document.views.findIndex(view=>view.id===id);
  if(from<0||target<0||target>=document.views.length||from===target)return;
  const views=[...document.views],picked=views.splice(from,1)[0];views.splice(target,0,picked);
  update({...document,views});
 };
 const insertionAt=(y:number)=>{
  const rows=listRef.current?.querySelectorAll<HTMLLIElement>('li[data-view-id]');
  if(!rows)return 0;
  for(let index=0;index<rows.length;index++){const rect=rows[index].getBoundingClientRect();if(y<rect.top+rect.height/2)return index;}
  return rows.length;
 };
 const startDrag=(event:ReactPointerEvent<HTMLButtonElement>,view:SavedView,index:number)=>{
  if(event.button!==0)return;
  event.preventDefault();event.currentTarget.setPointerCapture(event.pointerId);
  dragRef.current={id:view.id,pointerId:event.pointerId,startY:event.clientY,insertion:index};
  setDraggingId(view.id);
 };
 const dragMove=(event:ReactPointerEvent<HTMLButtonElement>)=>{
  const drag=dragRef.current;if(!drag||drag.pointerId!==event.pointerId)return;
  event.preventDefault();
  const list=listRef.current;if(list){const rect=list.getBoundingClientRect();if(event.clientY<rect.top+12)list.scrollTop-=14;else if(event.clientY>rect.bottom-12)list.scrollTop+=14;}
  if(Math.abs(event.clientY-drag.startY)<4&&dropIndex===null)return;
  drag.insertion=insertionAt(event.clientY);setDropIndex(drag.insertion);
 };
 const endDrag=(event:ReactPointerEvent<HTMLButtonElement>,cancel=false)=>{
  const drag=dragRef.current;if(!drag||drag.pointerId!==event.pointerId)return;
  if(event.currentTarget.hasPointerCapture(event.pointerId))event.currentTarget.releasePointerCapture(event.pointerId);
  if(!cancel&&Math.abs(event.clientY-drag.startY)>=4){
   const from=document.views.findIndex(view=>view.id===drag.id),insertion=insertionAt(event.clientY);
   moveView(drag.id,insertion>from?insertion-1:insertion);
  }
  dragRef.current=null;setDraggingId(null);setDropIndex(null);
 };
 const typePicker=<select id="advanced-type" aria-label="Advanced tool type" value={document.type} onChange={event=>{const type=event.target.value as ToolType;if(type==='connect'&&connection.active&&connection.role==='host')connection.openControls();else connection.closeControls();setSlideIndex(null);update({...document,type});}}><option value="slides">Slides</option><option value="quiz">Quiz</option><option value="layers">Layers</option><option value="connect">Connect</option></select>;
 return <>
  <section className={`advanced-panel ${presenting?'is-presenting':'glass'}`} aria-label="Advanced tools">
   {document.type==='connect'?<div className="advanced-connect">{(!connection.active||connection.role==='host'&&connection.controlsOpen)&&<div className="connect-tool-type">{typePicker}</div>}<ConnectTools connection={connection}/></div>:slidePlaying&&slideIndex!==null?<div className="advanced-playback" role="toolbar" aria-label="Slide playback">
    <Button variant="ghost" onClick={()=>showSlide(0)} disabled={slideIndex===0} aria-label="First slide" title="First slide"><SkipBack size={18}/></Button>
    <Button variant="ghost" onClick={()=>showSlide(slideIndex-1)} aria-label="Previous slide" title="Previous slide (←)"><StepBack size={19}/></Button>
    <div className="advanced-slide-caption" aria-live="polite"><span className="advanced-slide-count">{slideIndex+1}/{document.views.length}</span><span className="advanced-slide-name" title={document.views[slideIndex].name}>{document.views[slideIndex].name}</span></div>
    <Button variant="ghost" onClick={stopSlides} aria-label="Stop slides" title="Stop slides (↓)"><Square size={18} fill="none"/></Button>
    <Button variant="ghost" onClick={()=>showSlide(slideIndex+1)} aria-label="Next slide" title="Next slide (Space or →)"><StepForward size={19}/></Button>
    <Button variant="ghost" onClick={()=>showSlide(document.views.length-1)} disabled={slideIndex>=document.views.length-1} aria-label="Last slide" title="Last slide"><SkipForward size={18}/></Button>
   </div>:<div className="advanced-grid">
    <div className="advanced-controls">
     {typePicker}
     <div className="advanced-add-row">
      <input type="text" aria-label="Name for new view" placeholder="View name" value={newViewName} onChange={event=>setNewViewName(event.target.value)}/>
      <Button variant="outline" onClick={addView} disabled={!ready} aria-label="Add present view" title="Add present view"><Plus size={18}/></Button>
     </div>
     <div className="advanced-actions" role="toolbar" aria-label="Advanced tool actions">
      <Button variant="outline" onClick={openEditor} aria-label="Open JSON editor" title="Open JSON editor"><span className="advanced-code-icon" aria-hidden="true">{'<>'}</span></Button>
      <Button variant="outline" onClick={()=>showSlide(0)} disabled={document.type!=='slides'||document.views.length===0} aria-label="Play slides" title="Play slides"><Play size={17}/></Button>
     </div>
    </div>
    <div ref={listRef} className="advanced-view-list" aria-label="Added views"><ol>{document.views.map((view,index)=><li key={view.id} data-view-id={view.id} data-dragging={draggingId===view.id} data-drop-before={dropIndex===index} data-drop-after={dropIndex===document.views.length&&index===document.views.length-1}>
     <button type="button" className="advanced-drag-handle" aria-label={`Reorder ${view.name}`} title="Drag to reorder; arrow keys also work" onPointerDown={event=>startDrag(event,view,index)} onPointerMove={dragMove} onPointerUp={event=>endDrag(event)} onPointerCancel={event=>endDrag(event,true)} onKeyDown={event=>{if(event.key==='ArrowUp'||event.key==='ArrowDown'){event.preventDefault();moveView(view.id,index+(event.key==='ArrowUp'?-1:1));}}}><GripVertical size={15}/></button>
     {editingView?.id===view.id?<div className="advanced-view-open advanced-view-edit"><span className="advanced-view-number">{index+1}</span><input autoFocus aria-label="Slide name" value={editingView.name} onFocus={event=>event.currentTarget.select()} onChange={event=>setEditingView({id:view.id,name:event.target.value})} onBlur={saveViewName} onKeyDown={event=>{if(event.nativeEvent.isComposing)return;if(event.key==='Enter'){event.preventDefault();event.stopPropagation();event.currentTarget.blur();}else if(event.key==='Escape'){event.preventDefault();event.stopPropagation();setEditingView(null);}}}/></div>:<button type="button" className="advanced-view-open" onClick={()=>clickView(view)} onDoubleClick={()=>editView(view)} onKeyDown={event=>{if(event.key==='F2'){event.preventDefault();editView(view);}}} title={`Open ${view.name}; double-click to rename`}><span className="advanced-view-number">{index+1}</span><span className="advanced-view-name">{view.name}</span></button>}
     <Button type="button" variant="ghost" className="advanced-view-delete" onClick={()=>removeView(view.id)} aria-label={`Delete slide ${view.name}`} title={`Remove ${view.name} from set`}><Trash2 size={15}/></Button>
    </li>)}</ol></div>
   </div>}
   {!presenting&&<Button variant="ghost" className="advanced-close" onClick={()=>{connection.closeControls();close();}} aria-label="Close advanced tools" title="Close advanced tools"><X size={17}/></Button>}
  </section>
  <Dialog open={editorOpen} onOpenChange={setEditorOpen}><DialogContent className="advanced-json-dialog" showCloseButton={false}><DialogHeader><DialogTitle>Advanced tools JSON</DialogTitle><DialogDescription>Edit the current type and saved views.</DialogDescription></DialogHeader><textarea aria-label="Advanced tools JSON" spellCheck={false} value={editorText} onChange={event=>{setEditorText(event.target.value);setEditorError('');}}/>{editorError&&<p className="advanced-json-error" role="alert">{editorError}</p>}<DialogFooter><Button variant="outline" onClick={()=>setEditorOpen(false)}>Cancel</Button><Button onClick={applyEditor}>Apply JSON</Button></DialogFooter></DialogContent></Dialog>
 </>;
}
