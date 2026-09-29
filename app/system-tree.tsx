import {useMemo,useState,type CSSProperties,type ReactNode} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {SYSTEMS,type Atlas,type Part,type SceneState} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {partLayerOpacity} from './depth-control';
import {REGION_ORDER,MAJOR_SYSTEMS,buildAnatomyNodes,entryLabel,hierarchyEntries,majorSystemFor,regionPathFor,systemPathFor,type AnatomyEntry,type AnatomyNode} from './anatomy-hierarchy';
import {terminologyForGroup,terminologyTitle,type Terminology} from './anatomical-terminology';
import {anatomyNodeChoice,hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';
import {displayLaterality,lateralityClass} from './laterality';

type Update=SceneState|((state:SceneState)=>SceneState);
export type Hierarchy='systems'|'regions';
type Branch={id:string;name:string;color?:string;nodes:AnatomyNode[];parts:Part[];terminology:Terminology};
const CELL_COMPARTMENTS=['Cell boundary','Nucleus','Cytoplasm'];
function buildTree(entries:AnatomyEntry[],mode:Hierarchy,scope?:string):Branch[]{
 const categories=mode==='systems'?(scope==='cell'?SYSTEMS:MAJOR_SYSTEMS.map(system=>({...system,color:SYSTEMS.find(source=>source.id===system.source)?.color})))
  :(scope==='cell'?CELL_COMPARTMENTS:REGION_ORDER).map(name=>({id:name,name,color:undefined}));
 return categories.map(category=>{
  const children=entries.filter(entry=>mode==='systems'?(scope==='cell'?entry.system===category.id:majorSystemFor(entry)===category.id):entry.region===category.id);
  const path=(entry:AnatomyEntry)=>mode==='systems'?systemPathFor(entry,scope):regionPathFor(entry);
  return {id:category.id,name:category.name,color:category.color,nodes:buildAnatomyNodes(children,path),parts:children.flatMap(entry=>entry.parts),terminology:terminologyForGroup(category.name)};
 }).filter(branch=>branch.parts.length);
}
function TreeCheck({name,enabled,total,onClick}: {name:string;enabled:number;total:number;onClick:()=>void}){
 const status=enabled===0?'false':enabled===total?'true':'mixed';
 return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${name}`} aria-checked={status} data-state={status} onClick={onClick} title={enabled>0?`Hide ${name}`:`Show ${name}`}>
  {status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}
 </button>;
}

export default function SystemTree({atlas,mode,state,setState,onChoose}: {atlas:Atlas;mode:Hierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:HierarchyChoice,toggle?:boolean)=>void}){
 const entries=useMemo(()=>hierarchyEntries(atlas),[atlas]);
 const availableParts=useMemo(()=>atlas.parts.filter(part=>!part.suppressed),[atlas]);
 const tree=useMemo(()=>buildTree(entries,mode,atlas.scope),[entries,mode,atlas.scope]);
 const [expanded,setExpanded]=useState<Set<string>>(()=>new Set(['all']));
 const hidden=new Set(state.hidden??[]);
 const isOn=(part:Part)=>!part.suppressed&&state.visible.includes(part.system)&&!hidden.has(part.id)&&!state.depthHidden?.includes(depthLayerFor(part))&&partLayerOpacity(part,state)>0;
 const countOn=(parts:Part[])=>parts.reduce((sum,part)=>sum+Number(isOn(part)),0);
 const toggleOpen=(id:string)=>setExpanded(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});
 const toggleAll=()=>setState(current=>{
  const active=[...new Set(availableParts.map(part=>part.system))];
  const anyOn=availableParts.some(part=>current.visible.includes(part.system)&&!(current.hidden??[]).includes(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&partLayerOpacity(part,current)>0);
  return {...current,visible:anyOn?[]:active,hidden:anyOn?current.hidden:[],depth:anyOn?current.depth:undefined,depthHidden:anyOn?current.depthHidden:[],skinOpacity:anyOn?current.skinOpacity:(current.skinOpacity??0)>0?current.skinOpacity:.1,selected:[],isolate:false};
 });
 const toggleParts=(parts:Part[])=>setState(current=>{
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const anyOn=parts.some(part=>current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&partLayerOpacity(part,current)>0);
  const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
  const restored=new Set(parts.map(depthLayerFor).filter(id=>current.depthHidden?.includes(id)));
  if(!anyOn)for(const part of availableParts)if((newlyEnabled.has(part.system)||restored.has(depthLayerFor(part)))&&!ids.has(part.id))hiddenNow.add(part.id);
  if(anyOn)for(const id of ids)hiddenNow.add(id);else for(const id of ids)hiddenNow.delete(id);
  const systems=[...new Set(parts.map(part=>part.system))];
  return {...current,visible:anyOn?current.visible:[...new Set([...current.visible,...systems])],hidden:[...hiddenNow],depthHidden:anyOn?current.depthHidden:current.depthHidden?.filter(id=>!restored.has(id as ReturnType<typeof depthLayerFor>)),skinOpacity:!anyOn&&parts.some(part=>depthLayerFor(part)==='skin')&&!(current.skinOpacity??0)?.1:current.skinOpacity,selected:anyOn?current.selected.filter(id=>!ids.has(id)):current.selected,isolate:false};
 });
 const rowStyle=(depth:number)=>({'--tree-depth':depth} as CSSProperties);
 const entryRow=(entry:AnatomyEntry,depth:number,key:string)=>{
  const name=entryLabel(entry),parts=entry.parts,display=displayLaterality(name);
  return <div className="tree-row tree-leaf" style={rowStyle(depth)} key={`${key}:${entry.id}`}>
   <button type="button" className="tree-label" title={terminologyTitle(name,entry.terminology,entry.id)} aria-label={name} onClick={event=>onChoose({id:entry.id,name,elements:parts.map(part=>part.id)},event.ctrlKey||event.metaKey)} onContextMenu={event=>{if(event.ctrlKey){event.preventDefault();onChoose({id:entry.id,name,elements:parts.map(part=>part.id)},true);}}}><span className={`tree-name ${lateralityClass(display.side)}`}>{display.label}</span>{parts.length>1&&<span className="tree-count">{parts.length}</span>}</button>
   <TreeCheck name={name} enabled={countOn(parts)} total={parts.length} onClick={()=>toggleParts(parts)}/>
  </div>;
 };
 const renderNode=(node:AnatomyNode,depth:number,parentKey:string):ReactNode=>{
  const key=`${parentKey}:${node.id}`;
  if(node.kind==='entry')return entryRow(node.entry,depth,key);
  const open=expanded.has(key),display=displayLaterality(node.name);
  return <div key={key}>
   <div className={`tree-row ${node.kind==='bilateral'?'tree-bilateral':'tree-region'}`} style={rowStyle(depth)}>
    <button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${node.name}`} data-connect-key={key} aria-expanded={open} onClick={()=>toggleOpen(key)}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
    <button type="button" className="tree-label" aria-label={node.name} onClick={event=>onChoose(anatomyNodeChoice(node,`${mode}:${key}`),event.ctrlKey||event.metaKey)} title={terminologyTitle(node.name,node.terminology,`hierarchy:${mode}:${key}`)}><span className={`tree-name ${lateralityClass(display.side)}`}>{display.label}</span><span className="tree-count">{node.parts.length.toLocaleString()}</span></button>
    <TreeCheck name={node.name} enabled={countOn(node.parts)} total={node.parts.length} onClick={()=>toggleParts(node.parts)}/>
   </div>
   {open&&(node.kind==='group'?node.nodes.map(child=>renderNode(child,depth+1,key)):node.entries.map(entry=>entryRow(entry,depth+1,key)))}
  </div>;
 };
 return <nav className="system-tree" aria-label={atlas.scope==='cell'?(mode==='systems'?'Cell functions and components':'Cell compartments and components'):mode==='systems'?'Anatomical systems and structures':'Anatomical regions and structures'}>
  <div className="tree-row tree-root" style={rowStyle(0)}>
   <button type="button" className="tree-expander" aria-label={`${expanded.has('all')?'Collapse':'Expand'} All`} data-connect-key="all" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}>{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
   <button type="button" className="tree-label" onClick={()=>onChoose(hierarchyChoice(`${mode}:all`,'All',availableParts,tree.map(branch=>hierarchyChoice(`${mode}:${branch.id}`,branch.name,branch.parts,branch.nodes.map((node,index)=>anatomyNodeChoice(node,`${mode}:${branch.id}:${index}`)),branch.terminology)),terminologyForGroup('All')))} title={terminologyTitle(atlas.scope==='cell'?'All cell components':'All · human body',terminologyForGroup('All'),`hierarchy:${mode}:all`)}><span className="tree-name">All</span><span className="tree-count">{availableParts.length.toLocaleString()}</span></button>
   <TreeCheck name={atlas.scope==='cell'?'all cell components':'all anatomy'} enabled={countOn(availableParts)} total={availableParts.length} onClick={toggleAll}/>
  </div>
  {expanded.has('all')&&tree.map(branch=>{
   const branchKey=`${mode}:${branch.id}`,open=expanded.has(branchKey);
   return <div key={branchKey}>
    <div className="tree-row tree-system" style={rowStyle(1)}>
     <button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${branch.name}`} data-connect-key={branchKey} aria-expanded={open} onClick={()=>toggleOpen(branchKey)}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
     <button type="button" className="tree-label" onClick={event=>onChoose(hierarchyChoice(`${mode}:${branchKey}`,branch.name,branch.parts,branch.nodes.map((node,index)=>anatomyNodeChoice(node,`${branchKey}:${index}`)),branch.terminology),event.ctrlKey||event.metaKey)} title={terminologyTitle(branch.name,branch.terminology,`hierarchy:${mode}:${branchKey}`)}>{branch.color&&<span className="tree-dot" style={{background:branch.color}}/>}<span className="tree-name">{branch.name}</span><span className="tree-count">{branch.parts.length}</span></button>
     <TreeCheck name={branch.name} enabled={countOn(branch.parts)} total={branch.parts.length} onClick={()=>toggleParts(branch.parts)}/>
    </div>
    {open&&branch.nodes.map(node=>renderNode(node,2,branchKey))}
   </div>;
  })}
 </nav>;
}
