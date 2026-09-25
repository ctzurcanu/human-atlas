import {useMemo,useState,type CSSProperties} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import type {Atlas,Concept,Part,SceneState} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {resolveGuestHierarchy,type GuestHierarchy,type ResolvedGuestNode} from './guest-hierarchy';

type Update=SceneState|((state:SceneState)=>SceneState);
type Props={atlas:Atlas;hierarchy:GuestHierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:Concept)=>void};
const style=(depth:number)=>({'--tree-depth':depth} as CSSProperties);

export default function GuestTree({atlas,hierarchy,state,setState,onChoose}:Props){
 const nodes=useMemo(()=>resolveGuestHierarchy(atlas,hierarchy),[atlas,hierarchy]);
 const all=useMemo(()=>[...new Map(nodes.flatMap(node=>node.parts).map(part=>[part.id,part])).values()],[nodes]);
 const [expanded,setExpanded]=useState<Set<string>>(()=>new Set(['all']));
 const hidden=new Set(state.hidden??[]),depthHidden=new Set(state.depthHidden??[]),visible=new Set(state.visible);
 const isOn=(part:Part)=>!part.suppressed&&visible.has(part.system)&&!hidden.has(part.id)&&!depthHidden.has(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(state.skinOpacity??1)>0);
 const toggleOpen=(id:string)=>setExpanded(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});
 const toggleParts=(parts:Part[])=>setState(current=>{
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const allOn=parts.every(part=>!part.suppressed&&current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(current.skinOpacity??1)>0));
  if(allOn){for(const id of ids)hiddenNow.add(id);return {...current,hidden:[...hiddenNow],selected:current.selected.filter(id=>!ids.has(id)),isolate:false};}
  const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
  const restoredDepth=new Set<string>(parts.map(depthLayerFor).filter(id=>current.depthHidden?.includes(id)));
  for(const part of atlas.parts)if((newlyEnabled.has(part.system)||restoredDepth.has(depthLayerFor(part)))&&!ids.has(part.id))hiddenNow.add(part.id);
  for(const id of ids)hiddenNow.delete(id);
  return {...current,visible:[...new Set([...current.visible,...parts.map(part=>part.system)])],hidden:[...hiddenNow],depthHidden:current.depthHidden?.filter(id=>!restoredDepth.has(id)),skinOpacity:parts.some(part=>depthLayerFor(part)==='skin')&&!(current.skinOpacity??1)?.1:current.skinOpacity,isolate:false};
 });
 const check=(name:string,parts:Part[])=>{
  if(!parts.length)return null;
  const enabled=parts.filter(isOn).length,status=enabled===0?'false':enabled===parts.length?'true':'mixed';
  return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${name}`} aria-checked={status} data-state={status} onClick={()=>toggleParts(parts)} title={enabled===parts.length?`Hide ${name}`:`Show ${name}`}>
   {status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}
  </button>;
 };
 const renderNode=(node:ResolvedGuestNode,depth:number)=>{
  const open=expanded.has(node.id),hasChildren=node.children.length>0;
  const chooseDirect=()=>onChoose({id:`guest:${hierarchy.id}:${node.id}`,name:node.name,elements:node.directParts.map(part=>part.id)});
  return <div key={node.id}>
   <div className={`tree-row ${hasChildren?'tree-region':'tree-leaf'} ${node.parts.length?'':'guest-unmodeled'}`} style={style(depth)}>
    {hasChildren&&<button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${node.name}`} aria-expanded={open} onClick={()=>toggleOpen(node.id)}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>}
    <button type="button" className="tree-label" title={node.directParts.length?`Select ${node.name}`:node.parts.length?node.name:`${node.name} · no matching structure in this model`} aria-expanded={hasChildren&&!node.directParts.length?open:undefined} onClick={()=>node.directParts.length?chooseDirect():hasChildren?toggleOpen(node.id):undefined}>
     <span className="tree-name">{node.name}</span>{node.parts.length>0&&<span className="tree-count">{node.parts.length.toLocaleString()}</span>}
    </button>
    {check(node.name,node.parts)}
   </div>
   {open&&node.children.map(child=>renderNode(child,depth+1))}
  </div>;
 };
 return <nav className="system-tree guest-tree" aria-label={`${hierarchy.name} hierarchy`}>
  <div className="tree-row tree-root" style={style(0)}>
   <button type="button" className="tree-label" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}><span className="tree-chevron">{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">All</span><span className="tree-count">{all.length.toLocaleString()}</span></button>
   {check('all anatomy',all)}
  </div>
  {expanded.has('all')&&nodes.map(node=>renderNode(node,1))}
  {hierarchy.source&&<a className="guest-source" href={hierarchy.source} target="_blank" rel="noopener noreferrer">Source: {hierarchy.name}</a>}
 </nav>;
}
