import {useMemo,useState,type CSSProperties} from 'react';
import {Check,ChevronDown,ChevronRight,Minus} from 'lucide-react';
import {SYSTEMS,type Atlas,type Concept,type Part,type SceneState,type SystemId} from './anatomy';

type Update=SceneState|((state:SceneState)=>SceneState);
export type Hierarchy='systems'|'regions';
type Entry={id:string;name:string;system:SystemId;region:string;parts:Part[]};
type EntryNode={kind:'entry';name:string;entry:Entry;parts:Part[]};
type BilateralNode={kind:'bilateral';name:string;entries:Entry[];parts:Part[]};
type Branch={id:string;name:string;color?:string;nodes:(EntryNode|BilateralNode)[];parts:Part[]};
const REGION_ORDER=['Head & neck','Torso & pelvis','Upper limbs','Lower limbs','Whole body & spanning'];
const CELL_COMPARTMENTS=['Cell boundary','Nucleus','Cytoplasm'];
const sortNames=(a:{name:string},b:{name:string})=>a.name.localeCompare(b.name,undefined,{numeric:true,sensitivity:'base'});
function labelFor(entry:Entry){
 const side=/^(left|right) deferent duct$/i.exec(entry.name);
 const name=side?`Ductus deferens (${side[1].toLowerCase()})`:entry.name;
 return name.replace(/^./,letter=>letter.toUpperCase());
}
function bilateral(name:string){
 const suffix=/^(.*) \((left|right)\)$/i.exec(name);
 if(suffix)return {base:suffix[1],side:suffix[2].toLowerCase()};
 const prefix=/^(left|right) (.+)$/i.exec(name);
 return prefix?{base:prefix[2].replace(/^./,letter=>letter.toUpperCase()),side:prefix[1].toLowerCase()}:null;
}
function groupEntries(entries:Entry[]):(EntryNode|BilateralNode)[]{
 const bilateralEntries=new Map<string,Entry[]>();
 for(const entry of entries){const parsed=bilateral(labelFor(entry));if(parsed){const list=bilateralEntries.get(parsed.base)??[];list.push(entry);bilateralEntries.set(parsed.base,list);}}
 const grouped=new Set<string>(),nodes:(EntryNode|BilateralNode)[]=[];
 for(const entry of entries){
  if(grouped.has(entry.id))continue;
  const parsed=bilateral(labelFor(entry)),siblings=parsed?bilateralEntries.get(parsed.base)??[]:[];
  if(parsed&&siblings.some(sibling=>bilateral(labelFor(sibling))?.side==='left')&&siblings.some(sibling=>bilateral(labelFor(sibling))?.side==='right')){
   siblings.forEach(sibling=>grouped.add(sibling.id));
   nodes.push({kind:'bilateral',name:parsed.base,entries:siblings.sort((a,b)=>labelFor(a).localeCompare(labelFor(b))),parts:siblings.flatMap(sibling=>sibling.parts)});
  }else nodes.push({kind:'entry',name:labelFor(entry),entry,parts:entry.parts});
 }
 return nodes.sort(sortNames);
}
function regionFor(parts:Part[],scope?:string){
 if(scope==='cell')return parts[0].groups?.find(group=>CELL_COMPARTMENTS.includes(group))??'Cytoplasm';
 const minY=Math.min(...parts.map(part=>part.bounds[0][1]));
 const maxY=Math.max(...parts.map(part=>part.bounds[1][1]));
 if(maxY-minY>1.05)return 'Whole body & spanning';
 const groups=new Set(parts.flatMap(part=>part.groups??[]).map(group=>group.toLowerCase()));
 const has=(pattern:RegExp)=>[...groups].some(group=>pattern.test(group));
 if(has(/^(head|neck|brain|face|skull|cranium|cerebrum|cerebellum|left head|right head|head and neck)$/))return 'Head & neck';
 if(has(/^(left |right )?(upper limb|arm|forearm|hand|wrist|shoulder)( region)?$/))return 'Upper limbs';
 if(has(/^(left |right )?(lower limb|leg|thigh|foot|ankle|talocrural region)( region)?$/))return 'Lower limbs';
 if(has(/^(trunk|torso|thorax|chest|abdomen|pelvis|pelvic region)$/))return 'Torso & pelvis';
 const name=parts[0].name.toLowerCase();
 if(/^(skin of body|whole body|body surface)$/.test(name))return 'Whole body & spanning';
 const x=parts.reduce((sum,part)=>sum+(part.bounds[0][0]+part.bounds[1][0])/2,0)/parts.length;
 const y=(minY+maxY)/2;
 return y<.73?'Lower limbs':Math.abs(x)>.19&&y>.7&&y<1.55?'Upper limbs':y>1.4?'Head & neck':'Torso & pelvis';
}
export function hierarchyEntries(atlas:Atlas):Entry[]{
 const names=new Map(atlas.concepts.map(concept=>[concept.id,concept.name]));
 const byConcept=new Map<string,Part[]>();
 for(const part of atlas.parts){const group=byConcept.get(part.conceptId)??[];group.push(part);byConcept.set(part.conceptId,group);}
 return [...byConcept].map(([id,parts])=>({id,name:names.get(id)??parts[0].name,system:parts[0].system,region:regionFor(parts,atlas.scope),parts}));
}
function buildTree(entries:Entry[],mode:Hierarchy,scope?:string):Branch[]{
 const categories=mode==='systems'?SYSTEMS.map(system=>({...system,id:system.id})):(scope==='cell'?CELL_COMPARTMENTS:REGION_ORDER).map(name=>({id:name,name,color:undefined}));
 return categories.map(category=>{
  const children=entries.filter(entry=>mode==='systems'?entry.system===category.id:entry.region===category.id);
  return {id:category.id,name:category.name,color:category.color,nodes:groupEntries(children),parts:children.flatMap(entry=>entry.parts)};
 }).filter(branch=>branch.parts.length);
}
function TreeCheck({name,enabled,total,onClick}: {name:string;enabled:number;total:number;onClick:()=>void}){
 const status=enabled===0?'false':enabled===total?'true':'mixed';
 return <button type="button" className="tree-check" role="checkbox" aria-label={`Show ${name}`} aria-checked={status} data-state={status} onClick={onClick} title={enabled===total?`Hide ${name}`:`Show ${name}`}>
  {status==='true'?<Check size={12}/>:status==='mixed'?<Minus size={12}/>:null}
 </button>;
}

export default function SystemTree({atlas,mode,state,setState,onChoose}: {atlas:Atlas;mode:Hierarchy;state:SceneState;setState:(update:Update)=>void;onChoose:(concept:Concept)=>void}){
 const entries=useMemo(()=>hierarchyEntries(atlas),[atlas]);
 const tree=useMemo(()=>buildTree(entries,mode,atlas.scope),[entries,mode,atlas.scope]);
 const [expanded,setExpanded]=useState<Set<string>>(()=>new Set(['all']));
 const hidden=new Set(state.hidden??[]);
 const isOn=(part:Part)=>state.visible.includes(part.system)&&!hidden.has(part.id);
 const countOn=(parts:Part[])=>parts.reduce((sum,part)=>sum+Number(isOn(part)),0);
 const toggleOpen=(id:string)=>setExpanded(current=>{const next=new Set(current);if(next.has(id))next.delete(id);else next.add(id);return next;});
 const toggleAll=()=>setState(current=>{
  const active=[...new Set(atlas.parts.map(part=>part.system))];
  const allOn=atlas.parts.every(part=>current.visible.includes(part.system)&&!(current.hidden??[]).includes(part.id));
  return {...current,visible:allOn?[]:active,hidden:allOn?current.hidden:[],selected:[],isolate:false};
 });
 const toggleSystem=(branch:Branch)=>setState(current=>{
  const system=branch.id as SystemId,hiddenNow=new Set(current.hidden??[]);
  const allOn=current.visible.includes(system)&&branch.parts.every(part=>!hiddenNow.has(part.id));
  if(allOn)return {...current,visible:current.visible.filter(id=>id!==system),selected:current.selected.filter(id=>!branch.parts.some(part=>part.id===id)),isolate:false};
  for(const part of branch.parts)hiddenNow.delete(part.id);
  return {...current,visible:[...new Set([...current.visible,system])],hidden:[...hiddenNow],isolate:false};
 });
 const toggleParts=(parts:Part[])=>setState(current=>{
  const ids=new Set(parts.map(part=>part.id)),hiddenNow=new Set(current.hidden??[]);
  const allOn=parts.every(part=>current.visible.includes(part.system)&&!hiddenNow.has(part.id));
  const newlyEnabled=new Set(parts.filter(part=>!current.visible.includes(part.system)).map(part=>part.system));
  if(!allOn)for(const part of atlas.parts)if(newlyEnabled.has(part.system)&&!ids.has(part.id))hiddenNow.add(part.id);
  if(allOn)for(const id of ids)hiddenNow.add(id);else for(const id of ids)hiddenNow.delete(id);
  const systems=[...new Set(parts.map(part=>part.system))];
  return {...current,visible:allOn?current.visible:[...new Set([...current.visible,...systems])],hidden:[...hiddenNow],selected:allOn?current.selected.filter(id=>!ids.has(id)):current.selected,isolate:false};
 });
 const rowStyle=(depth:number)=>({'--tree-depth':depth} as CSSProperties);
 const entryRow=(entry:Entry,depth:number)=>{
  const name=labelFor(entry),parts=entry.parts;
  return <div className="tree-row tree-leaf" style={rowStyle(depth)} key={entry.id}>
   <button type="button" className="tree-label" title={name} onClick={()=>onChoose({id:entry.id,name,elements:parts.map(part=>part.id)})}><span className="tree-name">{name}</span>{parts.length>1&&<span className="tree-count">{parts.length}</span>}</button>
   <TreeCheck name={name} enabled={countOn(parts)} total={parts.length} onClick={()=>toggleParts(parts)}/>
  </div>;
 };
 return <nav className="system-tree" aria-label={atlas.scope==='cell'?(mode==='systems'?'Cell functions and components':'Cell compartments and components'):mode==='systems'?'Anatomical systems and structures':'Anatomical regions and structures'}>
  <div className="tree-row tree-root" style={rowStyle(0)}>
   <button type="button" className="tree-label" aria-expanded={expanded.has('all')} onClick={()=>toggleOpen('all')}><span className="tree-chevron">{expanded.has('all')?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">All</span><span className="tree-count">{atlas.parts.length.toLocaleString()}</span></button>
   <TreeCheck name={atlas.scope==='cell'?'all cell components':'all anatomy'} enabled={countOn(atlas.parts)} total={atlas.parts.length} onClick={toggleAll}/>
  </div>
  {expanded.has('all')&&tree.map(branch=>{
   const branchKey=`${mode}:${branch.id}`,open=expanded.has(branchKey);
   return <div key={branchKey}>
    <div className="tree-row tree-system" style={rowStyle(1)}>
     <button type="button" className="tree-label" aria-expanded={open} onClick={()=>toggleOpen(branchKey)} title={branch.name}><span className="tree-chevron">{open?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span>{branch.color&&<span className="tree-dot" style={{background:branch.color}}/>}<span className="tree-name">{branch.name}</span><span className="tree-count">{branch.parts.length}</span></button>
     <TreeCheck name={branch.name} enabled={countOn(branch.parts)} total={branch.parts.length} onClick={()=>mode==='systems'?toggleSystem(branch):toggleParts(branch.parts)}/>
    </div>
    {open&&branch.nodes.map(node=>{
       if(node.kind==='entry')return entryRow(node.entry,2);
       const groupKey=`${branchKey}:${node.name}`,groupOpen=expanded.has(groupKey);
       return <div key={groupKey}>
        <div className="tree-row tree-bilateral" style={rowStyle(2)}>
         <button type="button" className="tree-label" aria-expanded={groupOpen} onClick={()=>toggleOpen(groupKey)} title={node.name}><span className="tree-chevron">{groupOpen?<ChevronDown size={13}/>:<ChevronRight size={13}/>}</span><span className="tree-name">{node.name}</span><span className="tree-count">{node.parts.length}</span></button>
         <TreeCheck name={node.name} enabled={countOn(node.parts)} total={node.parts.length} onClick={()=>toggleParts(node.parts)}/>
        </div>
        {groupOpen&&node.entries.map(entry=>entryRow(entry,3))}
       </div>;
      })}
   </div>;
  })}
 </nav>;
}
