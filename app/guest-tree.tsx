import {useEffect,useMemo,useState,type CSSProperties} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {structureName,type Atlas,type Part,type SceneState} from './anatomy';
import {depthLayerFor} from './depth-layers';
import {partLayerOpacity} from './depth-control';
import {resolveGuestHierarchy,type GuestHierarchy,type GuestLink,type ResolvedGuestNode} from './guest-hierarchy';
import {displayLaterality,lateralityClass} from './laterality';
import {guestNodeChoice,guestVisibleChildren} from './guest-choice';
import {hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';
import {terminologyForConcept,terminologyForGroup,terminologyTitle} from './anatomical-terminology';
import {dermatomeColors,dermatomeSelection,dermatomeTerritoriesFor} from './dermatome-colors';

type Update=SceneState|((state:SceneState)=>SceneState);
type Props={dermatomeReference?:boolean;atlas:Atlas;hierarchy:GuestHierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:HierarchyChoice,toggle?:boolean)=>void;onNavigate?:(link:GuestLink)=>void;onLoadExtension?:(id:string)=>Promise<void>};
const style=(depth:number)=>({'--tree-depth':depth} as CSSProperties);
const hierarchyNames:Record<string,string>={genes:'Genes','cell-types':'Cell Types',physiology:'Physiology','dermatomes-myotomes':'Dermatomes and Myotomes',drugs:'Drugs','physical-exercise':'Physical Exercise',chakras:'Chakras'};

export default function GuestTree({atlas,hierarchy,state,setState,onChoose,onNavigate,onLoadExtension,dermatomeReference}:Props){
 const nodes=useMemo(()=>resolveGuestHierarchy(atlas,hierarchy),[atlas,hierarchy]);
 const anatomyNames=useMemo(()=>new Map(atlas.concepts.map(concept=>[concept.id,concept.name])),[atlas]);
 const all=useMemo(()=>[...new Map(nodes.flatMap(node=>node.parts).map(part=>[part.id,part])).values()],[nodes]);
 const [extensionBusy,setExtensionBusy]=useState<Set<string>>(()=>new Set()),[extensionError,setExtensionError]=useState('');
 const loadExtension=async(id:string)=>{setExtensionBusy(current=>new Set([...current,id]));setExtensionError('');try{await onLoadExtension?.(id);}catch(cause){setExtensionError(cause instanceof Error?cause.message:'Could not load transcripts.');}finally{setExtensionBusy(current=>{const next=new Set(current);next.delete(id);return next;});}};
 const [expanded,setExpanded]=useState<Set<string>>(()=>new Set(['all'])),[searchClosed,setSearchClosed]=useState<Set<string>>(()=>new Set());
 useEffect(()=>setSearchClosed(new Set()),[state.guestQuery,hierarchy.roots]);
 const hidden=new Set(state.hidden??[]),depthHidden=new Set(state.depthHidden??[]),visible=new Set(state.visible);
 const isOn=(part:Part)=>!part.suppressed&&visible.has(part.system)&&!hidden.has(part.id)&&!depthHidden.has(depthLayerFor(part))&&partLayerOpacity(part,state)>0;
 const toggleOpen=(id:string)=>{if(state.guestQuery&&id!=='all'&&search?.included.has(id)){setSearchClosed(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});return;}setExpanded(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});};
 const toggleParts=(parts:Part[])=>setState(current=>{
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const allOn=parts.every(part=>!part.suppressed&&current.visible.includes(part.system)&&!hiddenNow.has(part.id)&&!current.depthHidden?.includes(depthLayerFor(part))&&partLayerOpacity(part,current)>0);
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
 const territoryCheck=(name:string,id:string)=>{
  const keys=dermatomeTerritoriesFor(id),hidden=new Set(state.dermatomeHidden??[]),enabled=keys.filter(key=>!hidden.has(key)).length;
  const status=enabled===0?'false':enabled===keys.length?'true':'mixed';
  return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${name}`} aria-checked={status} data-state={status} disabled={!keys.length} title={!keys.length?'No separate color territory in the source map':enabled===keys.length?`Hide ${name}`:`Show ${name}`} onClick={()=>setState(current=>{
   const next=new Set(current.dermatomeHidden??[]),allOn=keys.every(key=>!next.has(key));for(const key of keys){if(allOn)next.add(key);else next.delete(key);}return {...current,dermatomeHidden:[...next]};
  })}>{status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}</button>;
 };
 const usesTerritories=(id:string)=>import.meta.env.DEV&&hierarchy.id==='dermatomes-myotomes'&&(/^(DERMATOME|TRIGEMINAL)(:|$)/.test(id)||!!dermatomeReference&&/^(ROOT|REGION):/.test(id));
 const isCatalogueLeaf=(node:ResolvedGuestNode)=>hierarchy.id==='dermatomes-myotomes'&&/^(ANATOMY|INNERVATION|CRANIAL|AUTONOMIC):/.test(node.id)&&!node.children.length&&!!node.directParts.length&&node.directParts.every(part=>part.conceptId===node.directParts[0].conceptId);
 const guestChoice=(node:ResolvedGuestNode)=>isCatalogueLeaf(node)?{id:node.directParts[0].conceptId,name:node.name,elements:node.parts.map(part=>part.id),terminology:terminologyForConcept(node.directParts[0].conceptId)}:guestNodeChoice(hierarchy.id,node);
 const [limits,setLimits]=useState<Record<string,number>>({});
 const query=hierarchy.schema==='human-atlas-hierarchy/v2'?(state.guestQuery??'').trim().toLowerCase():'';
 const index=useMemo(()=>{
  const byId=new Map<string,ResolvedGuestNode>(),byQueryId=new Map<string,ResolvedGuestNode>(),parents=new Map<string,string[]>();
  const visit=(node:ResolvedGuestNode)=>{if(byId.has(node.id))return;byId.set(node.id,node);byQueryId.set(node.id.toLowerCase(),node);for(const child of node.children){parents.set(child.id,[...(parents.get(child.id)??[]),node.id]);visit(child);}};
  nodes.forEach(visit);return {byId,byQueryId,parents};
 },[nodes]);
 const search=useMemo(()=>{
  if(!query)return null;
  const exact=index.byQueryId.get(query);
  const matches=exact?[exact]:[...index.byId.values()].filter(node=>node.name.toLowerCase().includes(query)||node.id.toLowerCase().includes(query)||node.synonyms?.some(term=>term.toLowerCase().includes(query)));
  const included=new Set<string>();
  const include=(id:string)=>{if(included.has(id))return;included.add(id);for(const parent of index.parents.get(id)??[])include(parent);};
  matches.slice(0,200).forEach(node=>include(node.id));return {included,matches:new Set(matches.slice(0,200).map(node=>node.id)),total:matches.length};
 },[index,query]);
 const anatomyLeaves=(node:ResolvedGuestNode)=>{
  const groups=new Map<string,Part[]>();for(const part of node.directParts){const parts=groups.get(part.conceptId)??[];parts.push(part);groups.set(part.conceptId,parts);}
  return [...groups].map(([id,parts])=>({id,name:structureName(anatomyNames.get(id)??parts[0].name),parts}));
 };
 const renderedSearch=new Set<string>();
 const renderNode=(node:ResolvedGuestNode,depth:number,subtree=false)=>{
  if(search&&renderedSearch.has(node.id))return null;if(search)renderedSearch.add(node.id);
  const visibleChildren=guestVisibleChildren(hierarchy.id,node);
  const open=search?(search.included.has(node.id)&&!searchClosed.has(node.id))||expanded.has(node.id):expanded.has(node.id),hasChildren=visibleChildren.length>0||!!node.extension||(hierarchy.schema==='human-atlas-hierarchy/v2'&&(!isCatalogueLeaf(node)&&node.directParts.length>0||!!node.links?.length||hierarchy.id==='ta98'&&!!node.description)),display=displayLaterality(node.name);
  const chooseDirect=(toggle=false)=>onChoose(guestChoice(node),toggle);
  const includeSubtree=subtree||!!search?.matches.has(node.id)||(visibleChildren!==node.children&&node.children.some(child=>search?.matches.has(child.id))),showDetails=!search||includeSubtree;
  const leaves=open&&!isCatalogueLeaf(node)&&hierarchy.schema==='human-atlas-hierarchy/v2'&&showDetails?anatomyLeaves(node):[],leafLimit=limits[node.id+':anatomy']??150;
  const children=search&&!includeSubtree?visibleChildren.filter(child=>search.included.has(child.id)):visibleChildren,limit=limits[node.id]??150;
  return <div key={node.id}>
   <div className={`tree-row ${hasChildren?'tree-region':'tree-leaf'} ${node.parts.length?'':'guest-unmodeled'}`} style={style(depth)}>
    {hasChildren&&<button type="button" className="tree-expander" aria-label={`${open?'Collapse':'Expand'} ${node.name}`} data-connect-key={node.id} aria-expanded={open} onClick={()=>{if(!open&&node.extension)void loadExtension(node.id);toggleOpen(node.id);}}>{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>}
    <button type="button" className="tree-label" title={hierarchy.id==='ta98'||hierarchy.id==='dermatomes-myotomes'?node.description??node.name:node.name} aria-label={node.name} onClick={event=>node.parts.length||hierarchy.id==='ta98'||hierarchy.id==='dermatomes-myotomes'?chooseDirect(event.ctrlKey||event.metaKey):hasChildren?toggleOpen(node.id):undefined} onContextMenu={event=>{if(event.ctrlKey&&node.parts.length){event.preventDefault();chooseDirect(true);}}}>
     {hierarchy.id==='dermatomes-myotomes'&&dermatomeColors[dermatomeSelection(node.id).level??'']&&<span aria-hidden="true" style={{width:7,height:7,borderRadius:'50%',flexShrink:0,background:`rgb(${dermatomeColors[dermatomeSelection(node.id).level!].join(',')})`}}/>}<span className={`tree-name ${lateralityClass(display.side)}`}>{display.label}</span>{(usesTerritories(node.id)?dermatomeTerritoriesFor(node.id).length:node.parts.length)>0&&<span className="tree-count">{(usesTerritories(node.id)?dermatomeTerritoriesFor(node.id).length:node.parts.length).toLocaleString()}</span>}
    </button>
    {usesTerritories(node.id)?territoryCheck(node.name,node.id):check(node.name,node.parts)}
   </div>
   {open&&<>
    {hierarchy.id==='ta98'&&showDetails&&node.description&&<div className="guest-node-info" style={style(depth+1)}><p style={{whiteSpace:'pre-line'}}>{node.description}</p></div>}
    {hierarchy.schema==='human-atlas-hierarchy/v2'&&showDetails&&(node.extension||!!node.links?.length)&&<div className="guest-node-info" style={style(depth+1)}>
     {node.extension&&<button type="button" className="guest-more" disabled={extensionBusy.has(node.id)} onClick={()=>void loadExtension(node.id)}>{extensionBusy.has(node.id)?'Loading transcripts…':'Load transcripts and protein isoforms'}</button>}
     {!!node.links?.length&&<div className="guest-biology-links" aria-label={hierarchy.id==='genes'?'Associated cell types':hierarchy.id==='cell-types'?'Marker genes':'Related hierarchies'}>{node.links.map(link=><button type="button" key={`${link.hierarchy}:${link.id}`} onClick={()=>onNavigate?.(link)} title={`Open ${link.name} in ${hierarchyNames[link.hierarchy]??link.hierarchy}`}>{link.name}</button>)}</div>}
    </div>}
    {children.slice(0,limit).map(child=>renderNode(child,depth+1,includeSubtree))}
    {children.length>limit&&<button type="button" className="guest-more" onClick={()=>setLimits(current=>({...current,[node.id]:limit+150}))}>Show more ({(children.length-limit).toLocaleString()} remaining)</button>}
    {leaves.slice(0,leafLimit).map(leaf=><div className="tree-row tree-leaf guest-anatomy-leaf" style={style(depth+1)} key={`anatomy:${leaf.id}`}>
     <button type="button" className="tree-label" onClick={event=>onChoose({id:leaf.id,name:leaf.name,elements:leaf.parts.map(part=>part.id),terminology:terminologyForConcept(leaf.id)},event.ctrlKey||event.metaKey)}><span className="tree-name">{leaf.name}</span><span className="tree-count">{leaf.parts.length}</span></button>{check(leaf.name,leaf.parts)}
    </div>)}
    {leaves.length>leafLimit&&<button type="button" className="guest-more" onClick={()=>setLimits(current=>({...current,[node.id+':anatomy']:leafLimit+150}))}>Show more anatomy ({(leaves.length-leafLimit).toLocaleString()} remaining)</button>}
   </>}
  </div>;
 };
 return <nav className={`system-tree guest-tree ${hierarchy.schema==='human-atlas-hierarchy/v2'?'guest-biology':''}`} aria-label={`${hierarchy.name} hierarchy`}>
  {hierarchy.schema==='human-atlas-hierarchy/v2'&&<div className="guest-biology-controls">
   {!!hierarchy.views?.length&&<div className="guest-view-tabs" role="group" aria-label={`${hierarchy.name} view`}>{hierarchy.views.map(view=><button type="button" key={view.id} className={(state.guestView??hierarchy.views![0].id)===view.id?'active':''} onClick={()=>setState(current=>({...current,guestView:view.id,guestQuery:undefined}))}>{view.name}</button>)}</div>}
   <input type="search" aria-label={`Search ${hierarchy.name}`} placeholder={`Search ${hierarchy.name.toLowerCase()}…`} value={state.guestQuery??''} onChange={event=>setState(current=>({...current,guestQuery:event.target.value}))}/>
   {search&&<span role="status">{search.total.toLocaleString()} matches{search.total>200?' · first 200 shown':''}</span>}
  </div>}
  {import.meta.env.DEV&&hierarchy.id==='dermatomes-myotomes'&&<p className="search-note">Dermatome overlay · 80% · approximate territories</p>}
  <div className="tree-row tree-root" style={style(0)}>
   <button type="button" className="tree-expander" aria-label={`${expanded.has('all')?'Collapse':'Expand'} All`} data-connect-key="all" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}>{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</button>
   <button type="button" className="tree-label" onClick={()=>onChoose(hierarchyChoice(`guest:${hierarchy.id}:all`,'All',all,nodes.filter(node=>hierarchy.id==='ta98'||node.parts.length).map(guestChoice),terminologyForGroup('All')))} title={terminologyTitle('All',terminologyForGroup('All'),`hierarchy:guest:${hierarchy.id}:all`)}><span className="tree-name">All</span><span className="tree-count">{(dermatomeReference?dermatomeTerritoriesFor('all').length:all.length).toLocaleString()}</span></button>
   {dermatomeReference?territoryCheck('all dermatomes','all'):check('all anatomy',all)}
  </div>
  {expanded.has('all')&&nodes.filter(node=>!search||search.included.has(node.id)).map(node=>renderNode(node,1))}
  {extensionError&&<p role="alert">{extensionError}</p>}
 </nav>;
}
