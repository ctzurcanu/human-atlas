import {useState} from 'react';
import {Plus,X} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {Dialog,DialogContent,DialogDescription,DialogFooter,DialogHeader,DialogTitle} from '@/components/ui/dialog';
import type {SceneState} from './anatomy';

type ToolType='slides'|'quiz'|'layers';
type CapturedScene={model:string;hierarchy:string;guestSources:string[];state:SceneState};
type SavedView={id:string;name:string;url:string;model:string;addedAt:string;scene?:CapturedScene};
type ToolDocument={schemaVersion:1;type:ToolType;views:SavedView[]};
const STORAGE_KEY='human-atlas-advanced-tools-v1';
const emptyDocument:ToolDocument={schemaVersion:1,type:'slides',views:[]};
const validType=(value:unknown):value is ToolType=>value==='slides'||value==='quiz'||value==='layers';
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

interface Props{close:()=>void;currentViewUrl:()=>string;captureScene:()=>CapturedScene;model:string;viewName:string}
export default function AdvancedTools({close,currentViewUrl,captureScene,model,viewName}:Props){
 const [document,setDocument]=useState<ToolDocument>(readDocument);
 const [editorOpen,setEditorOpen]=useState(false);
 const [editorText,setEditorText]=useState('');
 const [editorError,setEditorError]=useState('');
 const update=(next:ToolDocument)=>{setDocument(next);try{localStorage.setItem(STORAGE_KEY,JSON.stringify(next));}catch{/* The current session still retains the document. */}};
 const addView=()=>{
  const next:SavedView={id:crypto.randomUUID(),name:viewName||`View ${document.views.length+1}`,url:currentViewUrl(),model,addedAt:new Date().toISOString(),scene:captureScene()};
  update({...document,views:[...document.views,next]});
 };
 const openEditor=()=>{setEditorText(JSON.stringify(document,null,2));setEditorError('');setEditorOpen(true);};
 const applyEditor=()=>{
  try{
   const parsed:unknown=JSON.parse(editorText);
   if(!validateDocument(parsed))throw new Error('Use schemaVersion 1, a valid type, and views with same-site URLs.');
   update(parsed);setEditorOpen(false);
  }catch(error){setEditorError(error instanceof Error?error.message:'Invalid JSON.');}
 };
 const openView=(view:SavedView)=>{location.assign(view.url);};
 return <>
  <section className="advanced-panel glass" aria-label="Advanced tools">
   <div className="advanced-grid">
    <div className="advanced-controls">
     <select id="advanced-type" aria-label="Advanced tool type" value={document.type} onChange={event=>update({...document,type:event.target.value as ToolType})}><option value="slides">Slides</option><option value="quiz">Quiz</option><option value="layers">Layers</option></select>
     <div className="advanced-actions" role="toolbar" aria-label="Advanced tool actions">
      <Button variant="outline" onClick={openEditor} aria-label="Open JSON editor" title="Open JSON editor"><span className="advanced-code-icon" aria-hidden="true">{'<>'}</span></Button>
      <Button variant="outline" onClick={addView} aria-label="Add present view" title="Add present view"><Plus size={18}/><span className="sr-only">Add view</span></Button>
      <Button variant="ghost" onClick={close} aria-label="Close advanced tools" title="Close advanced tools"><X size={17}/></Button>
     </div>
    </div>
    <div className="advanced-view-list" aria-label="Added views"><ol>{document.views.map((view,index)=><li key={view.id}><button type="button" onClick={()=>openView(view)} title={`Open ${view.name}`}><span className="advanced-view-number">{index+1}</span><span className="advanced-view-name">{view.name}</span></button></li>)}</ol></div>
   </div>
  </section>
  <Dialog open={editorOpen} onOpenChange={setEditorOpen}><DialogContent className="advanced-json-dialog" showCloseButton={false}><DialogHeader><DialogTitle>Advanced tools JSON</DialogTitle><DialogDescription>Edit the current type and saved views.</DialogDescription></DialogHeader><textarea aria-label="Advanced tools JSON" spellCheck={false} value={editorText} onChange={event=>{setEditorText(event.target.value);setEditorError('');}}/>{editorError&&<p className="advanced-json-error" role="alert">{editorError}</p>}<DialogFooter><Button variant="outline" onClick={()=>setEditorOpen(false)}>Cancel</Button><Button onClick={applyEditor}>Apply JSON</Button></DialogFooter></DialogContent></Dialog>
 </>;
}
