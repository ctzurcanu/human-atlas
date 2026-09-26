import {useEffect,useState,type RefObject} from 'react';
import {Crosshair,ExternalLink,EyeOff,Focus} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {Slider} from '@/components/ui/slider';
import {enabledSections} from './section-stack';
import {SheetDescription,SheetTitle} from '@/components/ui/sheet';
import {SYSTEMS,explanation,structureName,surfaceRole,type Part,type SceneState} from './anatomy';
import {anatomyPath} from './anatomy-path';
import {localDescription,wikipediaDescription,type StructureDescription} from './structure-description';

interface Props {
 titleRef:RefObject<HTMLHeadingElement|null>;
 choice:{id:string;name:string};
 selectedParts:Part[];
 partById:Map<string,Part>;
 sex:'male'|'female';
 scope?:string;
 state:SceneState;
 covering:string[];
 depth:number;
 onOpacity:(value:number)=>void;
 onDepth:(value:number)=>void;
 onCenter:()=>void;
 onIsolate:()=>void;
 onHide:()=>void;
 onChoosePart:(id:string)=>void;
}

const wikipediaSearch=(term:string)=>`https://en.wikipedia.org/wiki/Special:Search?search=${encodeURIComponent(term.replace(/^\d+:\s*/,''))}`;
const wikipediaStructure=(name:string)=>{
 const term=structureName(name).replace(/\s*\((?:left|right)\)/gi,'').trim();
 if(term.toLowerCase()==='body of sternum')return 'https://en.wikipedia.org/wiki/Sternum#Body';
 if(term.toLowerCase()==='internal abdominal oblique muscle')return 'https://en.wikipedia.org/wiki/Abdominal_internal_oblique_muscle';
 return wikipediaSearch(term);
};
const sliderValue=(value:number|readonly number[])=>typeof value==='number'?value:value[0];

export default function SelectionInspector({titleRef,choice,selectedParts,partById,sex,scope,state,covering,depth,onOpacity,onDepth,onCenter,onIsolate,onHide,onChoosePart}:Props){
 const title=structureName(choice.name);
 const cell=scope==='cell';
 const selected=selectedParts[0],systems=[...new Set(selectedParts.map(part=>part.system))].map(id=>SYSTEMS.find(item=>item.id===id)).filter(item=>!!item),system=systems[0];
 const tags=[...new Set(selectedParts.flatMap(part=>[part.tissue,surfaceRole(part.name),...(part.regions??[])]).filter((tag):tag is string=>!!tag))];
 const path=anatomyPath(title,selectedParts,partById.values(),sex,scope);
 const next=partById.get(covering[depth]),previous=partById.get(covering[depth-1]);
 const context=Math.round((state.contextOpacity??1)*100);
 const [remote,setRemote]=useState<{title:string;value:StructureDescription}|null>(null);
 useEffect(()=>{
  if(explanation(title))return;
  let active=true;
  wikipediaDescription(title).then(value=>{if(active&&value)setRemote({title,value});});
  return()=>{active=false;};
 },[title]);
 const local=localDescription(title,selected,path);
 const resolved=!explanation(title)&&remote?.title===title?remote.value:local;
 const description=selected?resolved.text:'';

 return <>
  <div className="detail-header">
   <div className="detail-accent" style={{background:system?.color}}/>
   <div className="structure-title-row"><SheetTitle ref={titleRef} tabIndex={-1} className="structure-title">{title}</SheetTitle><a className="structure-wikipedia" href={wikipediaStructure(title)} target="_blank" rel="noopener noreferrer" aria-label={`Read about ${title} on Wikipedia`} title="Open Wikipedia in a new tab"><ExternalLink size={15}/></a></div>
  </div>
  <div className="detail-scroll" key={choice.id}>
   <nav className="structure-breadcrumbs" aria-label="Anatomical path">{path.map((node,index)=><span className="breadcrumb-step" key={`${index}:${node.label}`}>{index>0&&<span className="breadcrumb-separator" aria-hidden="true">›</span>}<a href={index===0?`https://en.wikipedia.org/wiki/${node.search}`:wikipediaSearch(node.search)} target="_blank" rel="noopener noreferrer" title={`Search Wikipedia for ${node.label}`}>{node.label}</a></span>)}</nav>
   {tags.length>0&&<div className="structure-tags">{tags.map(tag=><span className="structure-tag" key={tag}>{tag}</span>)}</div>}
   {!enabledSections(state).length&&<section className="inspection-section" aria-label={cell?'Surrounding components':'Surrounding anatomy'}>
    <div className="inspection-heading"><strong>{cell?'Surrounding components':'Surrounding anatomy'}</strong><output>{context}%</output></div>
    <Slider aria-label={cell?'Surrounding components opacity':'Surrounding anatomy opacity'} min={0} max={100} step={1} value={[context]} onValueChange={value=>onOpacity(sliderValue(value))}/>
    <div className="inspection-range"><span>Selection only</span><span>{cell?'Full cell':'Full anatomy'}</span></div>
   </section>}
   <section className="inspection-section" aria-label={cell?'Covering components':'Covering tissue'}>
    <div className="inspection-heading"><strong>{cell?'Covering components':'Covering tissue'}</strong></div>
    <div className="inspection-subheading"><span>Peel depth</span><output>{depth} / {covering.length}</output></div>
    <Slider aria-label="Covering tissue peel depth" min={0} max={covering.length} step={1} value={[depth]} disabled={!covering.length} onValueChange={value=>onDepth(sliderValue(value))}/>
    <div className="inspection-range"><span>Restore</span><span>Remove</span></div>
    <div className="peel-buttons">
     <Button variant="outline" disabled={depth===0} onClick={()=>onDepth(depth-1)}>{previous?`Restore ${structureName(previous.name)}`:'Nothing to restore'}</Button>
     <Button variant="outline" disabled={depth>=covering.length} onClick={()=>onDepth(depth+1)}>{next?`Remove ${structureName(next.name)}`:'Nothing to remove'}</Button>
    </div>
   </section>
   <SheetDescription className="structure-description">{description}</SheetDescription>
   {resolved.url&&<a className="context-note" href={resolved.url} target="_blank" rel="noopener noreferrer">Description from Wikipedia</a>}
   <div className="structure-meta"><span>Atlas reference<strong>{choice.id}</strong></span><span>Selected pieces<strong>{selectedParts.length.toLocaleString()}</strong></span></div>
   {new Set(selectedParts.map(part=>part.conceptId)).size>1&&<div className="member-list"><h3>Included structures</h3>{selectedParts.map(part=><Button variant="ghost" key={part.id} onClick={()=>onChoosePart(part.id)}><span>{structureName(part.name)}</span></Button>)}</div>}
  </div>
  <div className="detail-actions">
   <Button variant="outline" onClick={onCenter}><Crosshair size={16}/>Center selected structure</Button>
   <Button className={`primary-action ${state.isolate?'active':''}`} onClick={onIsolate}><Focus size={17}/>{state.isolate?`Show surrounding ${cell?'components':'anatomy'}`:`Show selected ${cell?'component':'structure'} only`}</Button>
   <Button variant="ghost" className="secondary-action" onClick={onHide}><EyeOff size={16}/>Hide selected {cell?'component':'structure'}</Button>
  </div>
 </>;
}
