import {useMemo,useState,type CSSProperties,type ReactNode} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {type Atlas,type Concept,type Part,type SceneState} from './anatomy';
import {DEPTH_LAYERS,depthLayerFor} from './depth-layers';
import {createDepthOrder} from './depth-sort';
import {buildAnatomyNodes,depthPathFor,entryLabel,hierarchyEntries,type AnatomyEntry,type AnatomyNode} from './anatomy-hierarchy';

type Update=SceneState|((state:SceneState)=>SceneState);
type Props={atlas:Atlas;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:Concept)=>void};

export default function DepthTree({atlas,state,setState,onChoose}:Props){
 const hierarchy=useMemo(()=>hierarchyEntries(atlas),[atlas]);
 const availableParts=useMemo(()=>atlas.parts.filter(part=>!part.suppressed),[atlas]);
 const compareEntries=useMemo(()=>createDepthOrder(atlas),[atlas]);
 const groups=useMemo(()=>{
  const byLayer=new Map<string,AnatomyEntry[]>(DEPTH_LAYERS.map(layer=>[layer.id,[]]));
  for(const entry of hierarchy){
   const partsByLayer=new Map<string,Part[]>();
   for(const part of entry.parts){const id=depthLayerFor(part),parts=partsByLayer.get(id)??[];parts.push(part);partsByLayer.set(id,parts);}
   for(const [id,parts] of partsByLayer)byLayer.get(id)!.push({...entry,parts});
  }
  return DEPTH_LAYERS.map((layer,index)=>{
   const entries=byLayer.get(layer.id)!;
   return {layer,index,parts:entries.flatMap(entry=>entry.parts),nodes:buildAnatomyNodes(entries,entry=>depthPathFor(entry,layer.name),compareEntries)};
  }).filter(group=>group.parts.length);
 },[hierarchy,compareEntries]);
 const [expanded,setExpanded]=useState<Set<string>>(()=>new Set(['all']));
 const hidden=new Set(state.hidden??[]),depthHidden=new Set(state.depthHidden??[]),visible=new Set(state.visible);
 const isOn=(part:Part)=>!part.suppressed&&!depthHidden.has(depthLayerFor(part))&&visible.has(part.system)&&!hidden.has(part.id)&&(depthLayerFor(part)!=='skin'||(state.skinOpacity??1)>0);
 const toggleOpen=(id:string)=>setExpanded(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});
 const toggleAll=()=>setState(current=>{
  const hiddenNow=new Set(current.hidden??[]),depthNow=new Set(current.depthHidden??[]),visibleNow=new Set(current.visible);
  const anyOn=availableParts.some(part=>!hiddenNow.has(part.id)&&!depthNow.has(depthLayerFor(part))&&visibleNow.has(part.system)&&(depthLayerFor(part)!=='skin'||(current.skinOpacity??1)>0));
  return {...current,visible:anyOn?current.visible:[...new Set(availableParts.map(part=>part.system))],hidden:anyOn?current.hidden:[],depthHidden:anyOn?DEPTH_LAYERS.map(layer=>layer.id):[],skinOpacity:anyOn?current.skinOpacity:(current.skinOpacity??0)>0?current.skinOpacity:.1,selected:[],isolate:false};
 });
 const toggleParts=(parts:Part[])=>setState(current=>{
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]),depthNow=new Set(current.depthHidden??[]);
  const anyOn=parts.some(part=>current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!depthNow.has(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(current.skinOpacity??1)>0));
  if(anyOn){for(const id of ids)hiddenNow.add(id);return {...current,hidden:[...hiddenNow],selected:current.selected.filter(id=>!ids.has(id)),isolate:false};}
  const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
  const restoredDepth=new Set(parts.map(depthLayerFor).filter(id=>depthNow.has(id)));
  for(const part of availableParts)if((newlyEnabled.has(part.system)||restoredDepth.has(depthLayerFor(part)))&&!ids.has(part.id))hiddenNow.add(part.id);
  for(const part of parts){hiddenNow.delete(part.id);depthNow.delete(depthLayerFor(part));}
  return {...current,visible:[...new Set([...current.visible,...parts.map(part=>part.system)])],hidden:[...hiddenNow],depthHidden:[...depthNow],skinOpacity:parts.some(part=>depthLayerFor(part)==='skin')&&!(current.skinOpacity??1) ? .1 : current.skinOpacity,isolate:false};
 });
 const toggleLayer=(id:string,parts:Part[])=>setState(current=>{
  const depthNow=new Set(current.depthHidden??[]);
  const anyOn=parts.some(part=>current.visible.includes(part.system)&&!(current.hidden??[]).includes(part.id)&&!depthNow.has(id)&&(id!=='skin'||(current.skinOpacity??1)>0));
  if(anyOn){depthNow.add(id);return {...current,depthHidden:[...depthNow],selected:current.selected.filter(selected=>!parts.some(part=>part.id===selected)),isolate:false};}
  depthNow.delete(id);
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
  for(const part of availableParts)if(newlyEnabled.has(part.system)&&!ids.has(part.id))hiddenNow.add(part.id);
  for(const id of ids)hiddenNow.delete(id);
  return {...current,visible:[...new Set([...current.visible,...parts.map(part=>part.system)])],hidden:[...hiddenNow],depthHidden:[...depthNow],skinOpacity:id==='skin'&&!(current.skinOpacity??1) ? .1 : current.skinOpacity,isolate:false};
 });
 const check=(label:string,parts:Part[],onClick:()=>void)=>{
  const enabled=parts.filter(isOn).length,status=enabled===0?'false':enabled===parts.length?'true':'mixed';
  return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${label}`} aria-checked={status} data-state={status} disabled={!parts.length} onClick={onClick} title={enabled>0?`Hide ${label}`:`Show ${label}`}>
   {status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}
  </button>;
 };
 const rowStyle=(depth:number)=>({'--tree-depth':depth} as CSSProperties);
 const entryRow=(entry:AnatomyEntry,depth:number,key:string)=>{
  const name=entryLabel(entry);
  return <div className="tree-row tree-leaf depth-leaf" style={rowStyle(depth)} key={`${key}:${entry.id}`}>
   <button type="button" className="tree-label" title={name} onClick={()=>onChoose({id:entry.id,name,elements:entry.parts.map(part=>part.id)})}><span className="tree-name">{name}</span>{entry.parts.length>1&&<span className="tree-count">{entry.parts.length}</span>}</button>
   {check(name,entry.parts,()=>toggleParts(entry.parts))}
  </div>;
 };
 const renderNode=(node:AnatomyNode,depth:number,parentKey:string):ReactNode=>{
  const key=`${parentKey}:${node.id}`;
  if(node.kind==='entry')return entryRow(node.entry,depth,key);
  const open=expanded.has(key);
  return <div key={key}>
   <div className={`tree-row ${node.kind==='bilateral'?'tree-bilateral':'tree-region'}`} style={rowStyle(depth)}>
    <button type="button" className="tree-label" aria-expanded={open} onClick={()=>toggleOpen(key)} title={node.name}><span className="tree-chevron">{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">{node.name}</span><span className="tree-count">{node.parts.length.toLocaleString()}</span></button>
    {check(node.name,node.parts,()=>toggleParts(node.parts))}
   </div>
   {open&&(node.kind==='group'?node.nodes.map(child=>renderNode(child,depth+1,key)):node.entries.map(entry=>entryRow(entry,depth+1,key)))}
  </div>;
 };
 return <nav className="system-tree depth-tree" aria-label="Anatomical depth layers and structures">
  <div className="tree-row tree-root" style={{'--tree-depth':0} as CSSProperties}>
   <button type="button" className="tree-label" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}><span className="tree-chevron">{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">All</span><span className="tree-count">{availableParts.length.toLocaleString()}</span></button>
   {check('all anatomy',availableParts,toggleAll)}
  </div>
  {expanded.has('all')&&groups.map(({layer,index,parts,nodes})=>{
   const open=expanded.has(layer.id);
   return <div key={layer.id}>
    <div className="tree-row tree-system depth-row" style={{'--tree-depth':0} as CSSProperties}>
     <button type="button" className="tree-label" aria-expanded={open} onClick={()=>toggleOpen(layer.id)} title={layer.name}><span className="depth-index">{String(index+1).padStart(2,'0')}</span><span className="tree-chevron">{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">{layer.name}</span><span className="tree-count">{parts.length.toLocaleString()}</span></button>
     {check(layer.name,parts,()=>toggleLayer(layer.id,parts))}
    </div>
    {open&&nodes.map(node=>renderNode(node,2,layer.id))}
   </div>;
  })}
 </nav>;
}
