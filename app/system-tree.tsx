import {useMemo,useState,type CSSProperties,type ReactNode} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {SYSTEMS,type Atlas,type Concept,type Part,type SceneState,type SystemId} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {REGION_ORDER,buildAnatomyNodes,displaySystemFor,entryLabel,hierarchyEntries,regionPathFor,systemPathFor,type AnatomyEntry,type AnatomyNode} from './anatomy-hierarchy';

type Update=SceneState|((state:SceneState)=>SceneState);
export type Hierarchy='systems'|'regions';
type Branch={id:string;name:string;color?:string;nodes:AnatomyNode[];parts:Part[]};
const CELL_COMPARTMENTS=['Cell boundary','Nucleus','Cytoplasm'];
function buildTree(entries:AnatomyEntry[],mode:Hierarchy,scope?:string):Branch[]{
 const categories=mode==='systems'?SYSTEMS.map(system=>({...system,id:system.id})):(scope==='cell'?CELL_COMPARTMENTS:REGION_ORDER).map(name=>({id:name,name,color:undefined}));
 return categories.map(category=>{
  const children=entries.filter(entry=>mode==='systems'?displaySystemFor(entry.system)===category.id:entry.region===category.id);
  const path=(entry:AnatomyEntry)=>mode==='systems'?systemPathFor(entry,scope):regionPathFor(entry);
  return {id:category.id,name:category.name,color:category.color,nodes:buildAnatomyNodes(children,path),parts:children.flatMap(entry=>entry.parts)};
 }).filter(branch=>branch.parts.length);
}
function TreeCheck({name,enabled,total,onClick}: {name:string;enabled:number;total:number;onClick:()=>void}){
 const status=enabled===0?'false':enabled===total?'true':'mixed';
 return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${name}`} aria-checked={status} data-state={status} onClick={onClick} title={enabled>0?`Hide ${name}`:`Show ${name}`}>
  {status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}
 </button>;
}

export default function SystemTree({atlas,mode,state,setState,onChoose}: {atlas:Atlas;mode:Hierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:Concept)=>void}){
 const entries=useMemo(()=>hierarchyEntries(atlas),[atlas]);
 const availableParts=useMemo(()=>atlas.parts.filter(part=>!part.suppressed),[atlas]);
 const tree=useMemo(()=>buildTree(entries,mode,atlas.scope),[entries,mode,atlas.scope]);
 const [expanded,setExpanded]=useState<Set<string>>(()=>new Set(['all']));
 const hidden=new Set(state.hidden??[]);
 const isOn=(part:Part)=>!part.suppressed&&state.visible.includes(part.system)&&!hidden.has(part.id)&&!state.depthHidden?.includes(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(state.skinOpacity??1)>0);
 const countOn=(parts:Part[])=>parts.reduce((sum,part)=>sum+Number(isOn(part)),0);
 const toggleOpen=(id:string)=>setExpanded(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});
 const toggleAll=()=>setState(current=>{
  const active=[...new Set(availableParts.map(part=>part.system))];
  const anyOn=availableParts.some(part=>current.visible.includes(part.system)&&!(current.hidden??[]).includes(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(current.skinOpacity??1)>0));
  return {...current,visible:anyOn?[]:active,hidden:anyOn?current.hidden:[],depthHidden:anyOn?current.depthHidden:[],skinOpacity:anyOn?current.skinOpacity:(current.skinOpacity??0)>0?current.skinOpacity:.1,selected:[],isolate:false};
 });
 const toggleSystem=(branch:Branch)=>setState(current=>{
  const systems=new Set(branch.parts.map(part=>part.system)),ids=new Set(branch.parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const anyOn=branch.parts.some(part=>current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(current.skinOpacity??1)>0));
  if(anyOn)return {...current,visible:current.visible.filter(id=>!systems.has(id)),selected:current.selected.filter(id=>!ids.has(id)),isolate:false};
  const restored=new Set(branch.parts.map(depthLayerFor).filter(id=>current.depthHidden?.includes(id)));
  for(const part of availableParts)if(restored.has(depthLayerFor(part))&&!ids.has(part.id))hiddenNow.add(part.id);
  for(const part of branch.parts)hiddenNow.delete(part.id);
  return {...current,visible:[...new Set([...current.visible,...systems])],hidden:[...hiddenNow],depthHidden:current.depthHidden?.filter(id=>!restored.has(id as ReturnType<typeof depthLayerFor>)),skinOpacity:branch.parts.some(part=>depthLayerFor(part)==='skin')&&!(current.skinOpacity??0)?.1:current.skinOpacity,isolate:false};
 });
 const toggleParts=(parts:Part[])=>setState(current=>{
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const anyOn=parts.some(part=>current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&(depthLayerFor(part)!=='skin'||(current.skinOpacity??1)>0));
  const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
  const restored=new Set(parts.map(depthLayerFor).filter(id=>current.depthHidden?.includes(id)));
  if(!anyOn)for(const part of availableParts)if((newlyEnabled.has(part.system)||restored.has(depthLayerFor(part)))&&!ids.has(part.id))hiddenNow.add(part.id);
  if(anyOn)for(const id of ids)hiddenNow.add(id);else for(const id of ids)hiddenNow.delete(id);
  const systems=[...new Set(parts.map(part=>part.system))];
  return {...current,visible:anyOn?current.visible:[...new Set([...current.visible,...systems])],hidden:[...hiddenNow],depthHidden:anyOn?current.depthHidden:current.depthHidden?.filter(id=>!restored.has(id as ReturnType<typeof depthLayerFor>)),skinOpacity:!anyOn&&parts.some(part=>depthLayerFor(part)==='skin')&&!(current.skinOpacity??0)?.1:current.skinOpacity,selected:anyOn?current.selected.filter(id=>!ids.has(id)):current.selected,isolate:false};
 });
 const rowStyle=(depth:number)=>({'--tree-depth':depth} as CSSProperties);
 const entryRow=(entry:AnatomyEntry,depth:number,key:string)=>{
  const name=entryLabel(entry),parts=entry.parts;
  return <div className="tree-row tree-leaf" style={rowStyle(depth)} key={`${key}:${entry.id}`}>
   <button type="button" className="tree-label" title={name} onClick={()=>onChoose({id:entry.id,name,elements:parts.map(part=>part.id)})}><span className="tree-name">{name}</span>{parts.length>1&&<span className="tree-count">{parts.length}</span>}</button>
   <TreeCheck name={name} enabled={countOn(parts)} total={parts.length} onClick={()=>toggleParts(parts)}/>
  </div>;
 };
 const renderNode=(node:AnatomyNode,depth:number,parentKey:string):ReactNode=>{
  const key=`${parentKey}:${node.id}`;
  if(node.kind==='entry')return entryRow(node.entry,depth,key);
  const open=expanded.has(key);
  return <div key={key}>
   <div className={`tree-row ${node.kind==='bilateral'?'tree-bilateral':'tree-region'}`} style={rowStyle(depth)}>
    <button type="button" className="tree-label" aria-expanded={open} onClick={()=>toggleOpen(key)} title={node.name}><span className="tree-chevron">{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">{node.name}</span><span className="tree-count">{node.parts.length.toLocaleString()}</span></button>
    <TreeCheck name={node.name} enabled={countOn(node.parts)} total={node.parts.length} onClick={()=>toggleParts(node.parts)}/>
   </div>
   {open&&(node.kind==='group'?node.nodes.map(child=>renderNode(child,depth+1,key)):node.entries.map(entry=>entryRow(entry,depth+1,key)))}
  </div>;
 };
 return <nav className="system-tree" aria-label={atlas.scope==='cell'?(mode==='systems'?'Cell functions and components':'Cell compartments and components'):mode==='systems'?'Anatomical systems and structures':'Anatomical regions and structures'}>
  <div className="tree-row tree-root" style={rowStyle(0)}>
   <button type="button" className="tree-label" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}><span className="tree-chevron">{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">All</span><span className="tree-count">{availableParts.length.toLocaleString()}</span></button>
   <TreeCheck name={atlas.scope==='cell'?'all cell components':'all anatomy'} enabled={countOn(availableParts)} total={availableParts.length} onClick={toggleAll}/>
  </div>
  {expanded.has('all')&&tree.map(branch=>{
   const branchKey=`${mode}:${branch.id}`,open=expanded.has(branchKey);
   return <div key={branchKey}>
    <div className="tree-row tree-system" style={rowStyle(1)}>
     <button type="button" className="tree-label" aria-expanded={open} onClick={()=>toggleOpen(branchKey)} title={branch.name}><span className="tree-chevron">{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span>{branch.color&&<span className="tree-dot" style={{background:branch.color}}/>}<span className="tree-name">{branch.name}</span><span className="tree-count">{branch.parts.length}</span></button>
     <TreeCheck name={branch.name} enabled={countOn(branch.parts)} total={branch.parts.length} onClick={()=>mode==='systems'?toggleSystem(branch):toggleParts(branch.parts)}/>
    </div>
    {open&&branch.nodes.map(node=>renderNode(node,2,branchKey))}
   </div>;
  })}
 </nav>;
}
