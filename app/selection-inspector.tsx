import {useEffect,useMemo,useState,type RefObject} from 'react';
import {Crosshair,ExternalLink,EyeOff,Focus,Network} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {Slider} from '@/components/ui/slider';
import {Tooltip,TooltipContent,TooltipProvider,TooltipTrigger} from '@/components/ui/tooltip';
import {enabledSections} from './section-stack';
import {SheetDescription,SheetTitle} from '@/components/ui/sheet';
import {SYSTEMS,structureName,surfaceRole,type Part,type SceneState} from './anatomy';
import {anatomyPath} from './anatomy-path';
import {anatomicalRelations,isAnatomicalBone,type RelationKind} from './anatomical-relations';
import {localDescription,wikipediaDescription,type StructureDescription} from './structure-description';

interface Props {
 titleRef:RefObject<HTMLHeadingElement|null>;
 choice:{id:string;name:string};
 selectedParts:Part[];
 anchorPart?:Part;
 relationshipDepth:number;
 relationshipExhausted:boolean;
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
 onExpandRelationships:()=>void;
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
const relationLabels:Record<RelationKind,string>={before:'Before',after:'After',innervation:'Innervated by',arterial:'Arterial supply',venous:'Venous drainage',innervates:'Innervates',supplies:'Supplies',drains:'Drains',articulates:'Articulates with',connects:'Connects',connectedBy:'Ligaments / tendons',joint:'At this joint',continuous:'Continuous with',covers:'Covers',coveredBy:'Covered by',cartilages:'Cartilages',bones:'Bones',tendons:'Tendons',fascia:'Fascia',muscles:'Muscles',origin:'Origin on',insertion:'Inserts on'};
const relationName=(part:Part)=>`${structureName(part.name)}${surfaceRole(part.name)?` · ${surfaceRole(part.name)}`:''}`;

export default function SelectionInspector({titleRef,choice,selectedParts,anchorPart,relationshipDepth,relationshipExhausted,partById,sex,scope,state,covering,depth,onOpacity,onDepth,onCenter,onIsolate,onExpandRelationships,onHide,onChoosePart}:Props){
 const title=structureName(choice.name);
 const cell=scope==='cell';
 const inspectedParts=anchorPart?[anchorPart]:selectedParts;
 const selected=inspectedParts[0],systems=[...new Set(inspectedParts.map(part=>part.system))].map(id=>SYSTEMS.find(item=>item.id===id)).filter(item=>!!item),system=systems[0];
 const tags=[...new Set(inspectedParts.flatMap(part=>[part.tissue,surfaceRole(part.name),...(part.regions??[])]).filter((tag):tag is string=>!!tag))];
 const path=anatomyPath(title,inspectedParts,partById.values(),sex,scope);
 const next=partById.get(covering[depth]),previous=partById.get(covering[depth-1]);
 const context=Math.round((state.contextOpacity??1)*100);
 const [remote,setRemote]=useState<{title:string;value:StructureDescription}|null>(null);
 useEffect(()=>{
  let active=true;
  wikipediaDescription(title).then(value=>{if(active&&value)setRemote({title,value});});
  return()=>{active=false;};
 },[title]);
 const local=localDescription(title,selected,path);
 const resolved=remote?.title===title?remote.value:local;
 const description=selected?resolved.text:'';
 const relations=useMemo(()=>selected?anatomicalRelations(selected,partById.values()):[],[selected,partById]);
 const skeletal=selected?isAnatomicalBone(selected):false;
 const boneCategories:RelationKind[]=['cartilages','tendons','fascia','muscles','innervation','arterial','venous'];
 const availableKinds=(Object.keys(relationLabels) as RelationKind[]).filter(kind=>relations.some(relation=>relation.kind===kind)||skeletal&&boneCategories.includes(kind));
 const relationKinds=skeletal?[...boneCategories,...availableKinds.filter(kind=>!boneCategories.includes(kind))]:availableKinds;

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
   {selected&&<section className="anatomical-relations" aria-label="Anatomical relationships">
    <h3>Related structures</h3>
    {relationKinds.length?relationKinds.map(kind=><div className="anatomical-relation-row" key={kind}>
     <span className="anatomical-relation-label">{relationLabels[kind]}</span>
     <div className="anatomical-relation-links">{relations.filter(relation=>relation.kind===kind).map(({target,via,viaModeled})=><span className="anatomical-relation-item" key={target.id}><button type="button" onClick={()=>onChoosePart(target.id)} title={`Select ${relationName(target)}`}>{relationName(target)}</button>{via?.length&&<small>via {via.join(' → ')}{viaModeled?'':' (not modeled)'}</small>}</span>)}{skeletal&&!relations.some(relation=>relation.kind===kind)&&<span className="anatomical-relation-unmapped">No named atlas link</span>}</div>
    </div>):<p className="anatomical-relations-empty">No mapped relationships in this model.</p>}
   </section>}
   {anchorPart&&relationshipDepth>0&&<p className="relationship-selection-status">Related structures added · {relationshipDepth} {relationshipDepth===1?'step':'steps'}</p>}
   {!anchorPart&&new Set(selectedParts.map(part=>part.conceptId)).size>1&&<div className="member-list"><h3>Included structures</h3>{selectedParts.map(part=><Button variant="ghost" key={part.id} onClick={()=>onChoosePart(part.id)}><span>{structureName(part.name)}</span></Button>)}</div>}
  </div>
  <div className="detail-actions">
   <TooltipProvider delay={150}><div className="detail-action-pill" role="toolbar" aria-label="Selected structure tools">
    <Tooltip><TooltipTrigger render={<button type="button" onClick={onCenter} aria-label="Center selected structure"><Crosshair size={19}/></button>}/><TooltipContent>Center selected structure</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" className={state.isolate?'active':''} onClick={onIsolate} aria-pressed={state.isolate} aria-label={state.isolate?`Show surrounding ${cell?'components':'anatomy'}`:`Show selected ${cell?'component':'structure'} only`}><Focus size={20}/></button>}/><TooltipContent>{state.isolate?`Show surrounding ${cell?'components':'anatomy'}`:`Show selected ${cell?'component':'structure'} only`}</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" className={relationshipDepth>0?'active':''} onClick={onExpandRelationships} aria-disabled={relationshipExhausted} aria-label={relationshipExhausted?'No further related structures':relationshipDepth>0?'Add the next related structures':'Show related structures together'}><Network size={20}/></button>}/><TooltipContent>{relationshipExhausted?'No further related structures':relationshipDepth>0?'Add the next related structures':'Show related structures together'}</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" onClick={onHide} aria-label={`Hide selected ${cell?'component':'structure'}`}><EyeOff size={20}/></button>}/><TooltipContent>Hide selected {cell?'component':'structure'}</TooltipContent></Tooltip>
   </div></TooltipProvider>
  </div>
 </>;
}
