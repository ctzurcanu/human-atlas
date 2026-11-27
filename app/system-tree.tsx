import {meshCount,separateOrganModels} from './organ-models';
import SeparateOrganModelLinks from './separate-organ-model-links';
import HierarchyCount from './hierarchy-count';
import {toggleSystemParts} from './hierarchy-visibility';
import {useMemo,useState,type CSSProperties,type ReactNode} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {SYSTEMS,type Atlas,type Part,type SceneState} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {partLayerOpacity} from './depth-control';
import {REGION_ORDER,MAJOR_SYSTEMS,buildAnatomyNodes,entryLabel,hierarchyEntries,majorSystemFor,regionPathFor,systemPathFor,type AnatomyEntry,type AnatomyNode} from './anatomy-hierarchy';
import {terminologyForGroup,terminologyTitle,type Terminology} from './anatomical-terminology';
import {anatomyNodeChoice,hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';
import {displayLaterality,lateralityClass} from './laterality';
import {withPulmonarySystem} from './pulmonary-system';
import {withLiverSystem} from './liver-system';
import {withTongueSystem} from './tongue-system';
import {withHeartSystem} from './heart-system';
import {withSpinalCordSystem} from './spinal-cord-system';
import {withEyeSystem} from './eye-system';
import {withLingualArterialSystem} from './lingual-arterial-system';
import {withNativeForearmVenousEntries} from './native-forearm-venous-controls';
import {withSourceVesselFamilies} from './source-vessel-families';
import {anatomicalReferenceChoices as liverReferenceChoices,anatomicalReferenceActive as liverReferenceActive,toggleAnatomicalReferenceChoice as toggleLiverReferenceChoice} from './anatomical-reference-selection';

type Update=SceneState|((state:SceneState)=>SceneState);
export type Hierarchy='systems'|'regions';
type Branch={id:string;name:string;color?:string;nodes:AnatomyNode[];parts:Part[];terminology:Terminology};
const CELL_COMPARTMENTS=['Cell boundary','Nucleus','Cytoplasm'];
export function buildSystemTree(entries:AnatomyEntry[],mode:Hierarchy,atlas:Atlas):Branch[]{
 const scope=atlas.scope;
 const scopedEntries=mode==='systems'?withNativeForearmVenousEntries(atlas,entries):entries;
 const categories=mode==='systems'?(scope==='cell'?SYSTEMS:MAJOR_SYSTEMS.map(system=>({...system,color:SYSTEMS.find(source=>source.id===system.source)?.color})))
  :(scope==='cell'?CELL_COMPARTMENTS:REGION_ORDER).map(name=>({id:name,name,color:undefined}));
 return categories.map(category=>{
  const children=scopedEntries.filter(entry=>mode==='systems'?(scope==='cell'?entry.system===category.id:majorSystemFor(entry)===category.id):entry.region===category.id);
  const path=(entry:AnatomyEntry)=>mode==='systems'?systemPathFor(entry,scope):regionPathFor(entry);
  const initialNodes=buildAnatomyNodes(children,path),eyeAssembly=mode==='systems'&&category.id==='nervous';
  const nodes=eyeAssembly?withEyeSystem(atlas,withSpinalCordSystem(atlas,initialNodes)):initialNodes;
  const branchParts=eyeAssembly?[...new Map(nodes.flatMap(node=>node.parts).map(part=>[part.id,part])).values()]:children.flatMap(entry=>entry.parts);
   return {id:category.id,name:category.name,color:category.color,nodes:mode==='systems'&&category.id==='respiratory'?withPulmonarySystem(atlas,nodes):mode==='systems'&&category.id==='digestive'?withTongueSystem(atlas,withLiverSystem(atlas,nodes)):mode==='systems'&&category.id==='circulatory'?withSourceVesselFamilies(atlas,withHeartSystem(atlas,withLingualArterialSystem(atlas,nodes))):nodes,parts:branchParts,terminology:terminologyForGroup(category.name)};
 }).filter(branch=>branch.parts.length);
}
function TreeCheck({name,enabled,total,onClick}: {name:string;enabled:number;total:number;onClick:()=>void}){
 const status=enabled===0?'false':enabled===total?'true':'mixed';
 return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${name}`} aria-checked={status} data-state={status} disabled={!total} onClick={onClick} title={!total?'No geometry installed':enabled>0?`Hide ${name}`:`Show ${name}`}>
  {status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}
 </button>;
}

export default function SystemTree({atlas,mode,state,setState,onChoose,choice}: {atlas:Atlas;mode:Hierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:HierarchyChoice,toggle?:boolean)=>void;choice?:HierarchyChoice}){
 const entries=useMemo(()=>hierarchyEntries(atlas),[atlas]);
 const availableParts=useMemo(()=>atlas.parts.filter(part=>!part.suppressed),[atlas]);
 const tree=useMemo(()=>buildSystemTree(entries,mode,atlas),[entries,mode,atlas]);
 const references=useMemo(()=>liverReferenceChoices(atlas),[atlas]);
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
 const toggleParts=(parts:Part[])=>setState(current=>toggleSystemParts(current,parts,availableParts));
 const rowStyle=(depth:number)=>({'--tree-depth':depth} as CSSProperties);
 const entryRow=(entry:AnatomyEntry,depth:number,key:string)=>{
  const name=entryLabel(entry),parts=entry.parts,display=displayLaterality(name);
  const reference=references.get(entry.ta98Term??entry.terminology?.ta98??'');
  return <div key={`${key}:${entry.id}`}><div className={`tree-row tree-leaf ${parts.length?'':'guest-unmodeled'}`} style={rowStyle(depth)}>
   <button type="button" className="tree-label" disabled={!parts.length} title={terminologyTitle(name,entry.terminology,entry.id)+(!parts.length?'\nNo geometry installed':'')} aria-label={name} onClick={event=>onChoose({id:entry.id,name,elements:parts.map(part=>part.id),terminology:entry.terminology,ta98Term:entry.ta98Term,ta98Kind:entry.ta98Kind},event.ctrlKey||event.metaKey)} onContextMenu={event=>{if(event.ctrlKey){event.preventDefault();onChoose({id:entry.id,name,elements:parts.map(part=>part.id),terminology:entry.terminology,ta98Term:entry.ta98Term,ta98Kind:entry.ta98Kind},true);}}}><span className={`tree-name ${lateralityClass(display.side)}`}>{display.label}</span>{(meshCount(parts)>1||separateOrganModels(atlas,{id:entry.id,name,term:entry.ta98Term??entry.terminology?.ta98??undefined}).length>0)&&<span className="tree-count">{meshCount(parts)}</span>}</button>
   <TreeCheck name={name} enabled={reference?Number(liverReferenceActive(reference,choice)):countOn(parts)} total={reference?1:parts.length} onClick={()=>reference?onChoose(toggleLiverReferenceChoice(reference,choice)):toggleParts(parts)}/>
  </div><SeparateOrganModelLinks atlas={atlas} id={entry.id} name={name} term={entry.ta98Term??entry.terminology?.ta98??undefined} depth={depth+1}/></div>;
 };
 const renderNode=(node:AnatomyNode,depth:number,parentKey:string):ReactNode=>{
  const key=`${parentKey}:${node.id}`;
  if(node.kind==='entry')return entryRow(node.entry,depth,key);
  const open=expanded.has(key),display=displayLaterality(node.name);
  const reference=references.get(node.terminology?.ta98??'');
  return <div key={key}>
   <div className={`tree-row ${node.kind==='bilateral'?'tree-bilateral':'tree-region'} ${node.parts.length?'':'guest-unmodeled'}`} style={rowStyle(depth)}>
    <button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${node.name}`} data-connect-key={key} aria-expanded={open} onClick={()=>toggleOpen(key)}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
    <button type="button" className="tree-label" disabled={!node.parts.length} aria-label={node.name} onClick={event=>onChoose(anatomyNodeChoice(node,`${mode}:${key}`),event.ctrlKey||event.metaKey)} title={terminologyTitle(node.name,node.terminology,`hierarchy:${mode}:${key}`)+(!node.parts.length?'\nNo geometry installed':'')}><span className={`tree-name ${lateralityClass(display.side)}`}>{display.label}</span><HierarchyCount parts={node.parts} itemCount={node.kind==='group'?node.itemCount:undefined}/></button>
    <TreeCheck name={node.name} enabled={reference?Number(liverReferenceActive(reference,choice)):countOn(node.parts)} total={reference?1:node.parts.length} onClick={()=>reference?onChoose(toggleLiverReferenceChoice(reference,choice)):toggleParts(node.parts)}/>
   </div>
   {open&&<><SeparateOrganModelLinks atlas={atlas} id={node.id} name={node.name} term={node.terminology?.ta98??undefined} depth={depth+1}/>{node.kind==='group'?node.nodes.map(child=>renderNode(child,depth+1,key)):node.entries.map(entry=>entryRow(entry,depth+1,key))}</>}
  </div>;
 };
 return <nav className="system-tree" aria-label={atlas.scope==='cell'?(mode==='systems'?'Cell functions and components':'Cell compartments and components'):mode==='systems'?'Anatomical systems and structures':'Anatomical regions and structures'}>
  <div className="tree-row tree-root" style={rowStyle(0)}>
   <button type="button" className="tree-expander" aria-label={`${expanded.has('all')?'Collapse':'Expand'} All`} data-connect-key="all" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}>{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
   <button type="button" className="tree-label" onClick={()=>onChoose(hierarchyChoice(`${mode}:all`,'All',availableParts,tree.map(branch=>hierarchyChoice(`${mode}:${branch.id}`,branch.name,branch.parts,branch.nodes.map((node,index)=>anatomyNodeChoice(node,`${mode}:${branch.id}:${index}`)),branch.terminology)),terminologyForGroup('All')))} title={terminologyTitle(atlas.scope==='cell'?'All cell components':'All · human body',terminologyForGroup('All'),`hierarchy:${mode}:all`)}><span className="tree-name">All</span><span className="tree-count">{meshCount(availableParts).toLocaleString()}</span></button>
   <TreeCheck name={atlas.scope==='cell'?'all cell components':'all anatomy'} enabled={countOn(availableParts)} total={availableParts.length} onClick={toggleAll}/>
  </div>
  {expanded.has('all')&&tree.map(branch=>{
   const branchKey=`${mode}:${branch.id}`,open=expanded.has(branchKey);
   return <div key={branchKey}>
    <div className="tree-row tree-system" style={rowStyle(1)}>
     <button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${branch.name}`} data-connect-key={branchKey} aria-expanded={open} onClick={()=>toggleOpen(branchKey)}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
     <button type="button" className="tree-label" onClick={event=>onChoose(hierarchyChoice(`${mode}:${branchKey}`,branch.name,branch.parts,branch.nodes.map((node,index)=>anatomyNodeChoice(node,`${branchKey}:${index}`)),branch.terminology),event.ctrlKey||event.metaKey)} title={terminologyTitle(branch.name,branch.terminology,`hierarchy:${mode}:${branchKey}`)}>{branch.color&&<span className="tree-dot" style={{background:branch.color}}/>}<span className="tree-name">{branch.name}</span><span className="tree-count">{meshCount(branch.parts).toLocaleString()}</span></button>
     <TreeCheck name={branch.name} enabled={countOn(branch.parts)} total={branch.parts.length} onClick={()=>toggleParts(branch.parts)}/>
    </div>
    {open&&branch.nodes.map(node=>renderNode(node,2,branchKey))}
   </div>;
  })}
 </nav>;
}
