import {useMemo,useState,type CSSProperties} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {structureName,type Atlas,type Part,type SceneState} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {resolveGuestHierarchy,type GuestHierarchy,type ResolvedGuestNode} from './guest-hierarchy';
import {displayLaterality,lateralityClass} from './laterality';
import {hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';
import {terminologyForConcept,terminologyForGroup,terminologyTitle} from './anatomical-terminology';

type Update=SceneState|((state:SceneState)=>SceneState);
type Props={atlas:Atlas;hierarchy:GuestHierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:HierarchyChoice,toggle?:boolean)=>void};
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
 const guestChoice=(node:ResolvedGuestNode):HierarchyChoice=>{
  const nested=node.children.map(guestChoice),descendantIds=new Set(node.children.flatMap(child=>child.parts.map(part=>part.id)));
  const direct=[...new Map(node.directParts.filter(part=>!descendantIds.has(part.id)).map(part=>[part.conceptId,part])).values()].map(part=>({id:part.conceptId,name:structureName(part.name),elements:node.directParts.filter(member=>member.conceptId===part.conceptId&&!descendantIds.has(member.id)).map(member=>member.id),terminology:terminologyForConcept(part.conceptId)}));
  return hierarchyChoice(`guest:${hierarchy.id}:${node.id}`,node.name,node.parts,[...nested,...direct],terminologyForGroup(node.name));
 };
 const renderNode=(node:ResolvedGuestNode,depth:number)=>{
  const open=expanded.has(node.id),hasChildren=node.children.length>0,display=displayLaterality(node.name);
  const chooseDirect=(toggle=false)=>onChoose(guestChoice(node),toggle);
  return <div key={node.id}>
   <div className={`tree-row ${hasChildren?'tree-region':'tree-leaf'} ${node.parts.length?'':'guest-unmodeled'}`} style={style(depth)}>
    {hasChildren&&<button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${node.name}`} data-connect-key={node.id} aria-expanded={open} onClick={()=>toggleOpen(node.id)}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>}
    <button type="button" className="tree-label" title={node.parts.length?terminologyTitle(node.name,terminologyForGroup(node.name),`hierarchy:guest:${hierarchy.id}:${node.id}`):`${node.name} · no matching structure in this model`} aria-label={node.name} onClick={event=>node.parts.length?chooseDirect(event.ctrlKey||event.metaKey):undefined} onContextMenu={event=>{if(event.ctrlKey&&node.parts.length){event.preventDefault();chooseDirect(true);}}}>
     <span className={`tree-name ${lateralityClass(display.side)}`}>{display.label}</span>{node.parts.length>0&&<span className="tree-count">{node.parts.length.toLocaleString()}</span>}
    </button>
    {check(node.name,node.parts)}
   </div>
   {open&&node.children.map(child=>renderNode(child,depth+1))}
  </div>;
 };
 return <nav className="system-tree guest-tree" aria-label={`${hierarchy.name} hierarchy`}>
  <div className="tree-row tree-root" style={style(0)}>
   <button type="button" className="tree-expander" aria-label={`${expanded.has('all')?'Collapse':'Expand'} All`} data-connect-key="all" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}>{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
   <button type="button" className="tree-label" onClick={()=>onChoose(hierarchyChoice(`guest:${hierarchy.id}:all`,'All',all,nodes.map(guestChoice),terminologyForGroup('All')))} title={terminologyTitle('All',terminologyForGroup('All'),`hierarchy:guest:${hierarchy.id}:all`)}><span className="tree-name">All</span><span className="tree-count">{all.length.toLocaleString()}</span></button>
   {check('all anatomy',all)}
  </div>
  {expanded.has('all')&&nodes.map(node=>renderNode(node,1))}
  {hierarchy.source&&<a className="guest-source" href={hierarchy.source} target="_blank" rel="noopener noreferrer">Source: {hierarchy.name}</a>}
 </nav>;
}
