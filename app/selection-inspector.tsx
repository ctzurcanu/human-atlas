import {useEffect,useMemo,useRef,useState,type RefObject} from 'react';
import {Crosshair,EyeOff,Focus,Network} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {Accordion,AccordionContent,AccordionItem,AccordionTrigger} from '@/components/ui/accordion';
import {Slider} from '@/components/ui/slider';
import {Tooltip,TooltipContent,TooltipProvider,TooltipTrigger} from '@/components/ui/tooltip';
import {SheetDescription,SheetTitle} from '@/components/ui/sheet';
import {SYSTEMS,structureName,surfaceRole,type Part,type Concept,type SceneState} from './anatomy';
import {anatomicalRelations,anatomicalConceptFunctionalRelations,type RelationKind,type ResolvedRelation} from './anatomical-relations';
import {atlasIdentifier,identifierReference,terminologyForConcept,type IdentifierKind} from './anatomical-terminology';
import type {HierarchyChoice} from './hierarchy-choice';
import {ta98PartOfChoice} from './hierarchy-navigation';
import {anatomicalComponentIds,hasExpandableAnatomicalRelations} from './anatomical-components';
import {ta98ModeledParent} from './ta98-modeled-parent';
import {displayLaterality,lateralityClass} from './laterality';
import {localDescription,wikipediaDescription,type StructureDescription} from './structure-description';

interface Props {
 titleRef:RefObject<HTMLHeadingElement|null>;
 choice:HierarchyChoice;
 ancestors:HierarchyChoice[];
 selectedParts:Part[];
 anchorParts?:Part[];
 relationshipDepth:number;
 relationshipExhausted:boolean;
 partById:Map<string,Part>;
 concepts:Concept[];
 scope?:string;
 state:SceneState;
 viewUrl:string;
 covering:string[];
 depth:number;
 onDepth:(value:number)=>void;
 onCenter:()=>void;
 onIsolate:()=>void;
 onExpandRelationships:()=>void;
 onShowComponents:(ids:string[])=>void;
 onHide:()=>void;
 onHidePart:(id:string)=>void;
 onChoosePart:(id:string,toggle?:boolean)=>void;
 onChooseChild:(choice:HierarchyChoice,toggle?:boolean)=>void;
}

const wikipediaSearch=(term:string)=>`https://en.wikipedia.org/wiki/Special:Search?search=${encodeURIComponent(term.replace(/^\d+:\s*/,''))}`;
const wikipediaStructure=(name:string)=>{
 const term=structureName(name).replace(/\s*\((?:left|right)\)/gi,'').replace(/\.[eo]\d*[lr]$/i,'').replace(/\.[lr]$/i,'').trim();
 if(term.toLowerCase()==='body of sternum')return 'https://en.wikipedia.org/wiki/Sternum#Body';
 if(term.toLowerCase()==='internal abdominal oblique muscle')return 'https://en.wikipedia.org/wiki/Abdominal_internal_oblique_muscle';
 return wikipediaSearch(term);
};
const sliderValue=(value:number|readonly number[])=>typeof value==='number'?value:value[0];
const relationLabels:Record<RelationKind,string>={before:'Before',after:'After',innervation:'Innervated by',arterial:'Arterial supply',venous:'Venous drainage',innervates:'Innervates',supplies:'Supplies',drains:'Drains',articulates:'Articulates with',connects:'Connects',connectedBy:'Ligaments / tendons',joint:'At this joint',continuous:'Continuous with',covers:'Covers',coveredBy:'Covered by',cartilages:'Cartilages',bones:'Bones',tendons:'Tendons',fascia:'Fascia',muscles:'Muscles',origin:'Origin on',insertion:'Inserts on',partOf:'Part of',contains:'Contains',originFor:'Origin marker for',insertionFor:'Insertion marker for',originSites:'Origin markers',insertionSites:'Insertion markers',counterpart:'Opposite side'};
const relationName=(part:Part)=>`${structureName(part.name)}${surfaceRole(part.name)?` · ${surfaceRole(part.name)}`:''}`;
const relationCache=new WeakMap<Map<string,Part>,Map<string,ResolvedRelation[]>>();
function IdentifierPill({kind,value,href,language}:{kind:IdentifierKind;value:string;href:string;language?:string}){
 return <span role="listitem" className="identifier-item"><a className={`identifier-pill identifier-${kind.toLowerCase()}`} href={href} target="_blank" rel="noopener noreferrer" title={`Open ${kind} ${value}`}><strong>{kind}</strong><span lang={language}>{value}</span></a></span>;
}

export default function SelectionInspector({titleRef,choice,ancestors,selectedParts,anchorParts,relationshipDepth,relationshipExhausted,partById,concepts,scope,state,viewUrl,covering,depth,onDepth,onCenter,onIsolate,onExpandRelationships,onShowComponents,onHide,onHidePart,onChoosePart,onChooseChild}:Props){
 const title=structureName(choice.name);
 const group=!!choice.children;
 const hasGeometry=selectedParts.length>0;
 const titleSide=displayLaterality(title).side;
 const cell=scope==='cell';
 const inspectedParts=anchorParts?.length?anchorParts:selectedParts;
 const inspectedKey=inspectedParts.map(part=>part.id).join('|');
 const relationKey=`${choice.id}|${choice.name}|${!!anchorParts?.length}|${inspectedKey}`;
 const coverageLimitation=choice.coverageLimitation??[...new Set(inspectedParts.map(part=>part.coverageLimitation).filter(Boolean))].join(' ');
 const allParts=useMemo(()=>[...partById.values()],[partById]);
 const selected=inspectedParts[0],systems=[...new Set(inspectedParts.map(part=>part.system))].map(id=>SYSTEMS.find(item=>item.id===id)).filter(item=>!!item),system=systems[0];
 const choiceTerm=terminologyForConcept(choice.id,choice.ta98Term);
 const terminology=choice.terminology??(choiceTerm.ta98||choiceTerm.fma||choiceTerm.ontology||choice.elements.length!==1?choiceTerm:terminologyForConcept(selected?.conceptId??choice.id));
 const atlasId=atlasIdentifier(choice.id);
 const localIdentifier=atlasId?{kind:atlasId.slice(0,atlasId.indexOf(':')) as 'Atlas'|'HA-G',value:atlasId.slice(atlasId.indexOf(':')+1)}:null;
 const tags=group?[]:[...new Set(inspectedParts.flatMap(part=>[part.tissue,surfaceRole(part.name),...(part.regions??[])]).filter((tag):tag is string=>!!tag))];
 const path=ancestors.map(parent=>({label:parent.name,search:parent.name}));
 const partOfChoice=ta98PartOfChoice(choice,ancestors)??ta98ModeledParent(choice,concepts,partById);
 const next=partById.get(covering[depth]),previous=partById.get(covering[depth-1]);
 const reference=(kind:IdentifierKind)=>identifierReference(kind,terminology,viewUrl);
 const [remote,setRemote]=useState<{title:string;value:StructureDescription}|null>(null);
 const lastControlToggle=useRef<{id:string;at:number}|null>(null);
 const controlChoose=(id:string)=>{const now=performance.now(),last=lastControlToggle.current;if(last?.id===id&&now-last.at<150)return;lastControlToggle.current={id,at:now};onChoosePart(id,true);};
 useEffect(()=>{
  if(group)return;
  let active=true;
  wikipediaDescription(title).then(value=>{if(active&&value)setRemote({title,value});});
  return()=>{active=false;};
 },[title,group]);
 const local=group?null:localDescription(title,selected,path);
 const resolved=remote?.title===title?remote.value:local;
 const description=!hasGeometry?'No selectable mesh in this model.':group?`${selectedParts.length.toLocaleString()} modeled ${selectedParts.length===1?'piece':'pieces'} in this group.`:selected?resolved?.text??'':'';
 const [relationState,setRelationState]=useState<{key:string;model:Map<string,Part>;items:ResolvedRelation[]}|null>(null);
 useEffect(()=>{
  // Mapped child choices establish Contains links without expanding every
  // leaf of a large hierarchy group. Otherwise resolve the available organ
  // links, including single-mesh organs selected through a TA98 group node.
  if(group&&anatomicalComponentIds(choice,[],partById).length)return;
  let cancelled=false,timer:ReturnType<typeof setTimeout>|undefined;
  // Give the sheet a paint before building the relationship graph. The first
  // skeletal selection can otherwise hold up the entire details panel.
  const frame=requestAnimationFrame(()=>{timer=setTimeout(()=>{
   if(cancelled)return;
   let cache=relationCache.get(partById);
   if(!cache){cache=new Map();relationCache.set(partById,cache);}
   const seen=new Set<string>();
   const items=[...anatomicalConceptFunctionalRelations(choice,inspectedParts,allParts),...inspectedParts.flatMap(part=>{
    let relations=cache.get(part.id);
    if(!relations){relations=anatomicalRelations(part,allParts);cache.set(part.id,relations);}
    return relations;
   })].filter(relation=>{const key=`${relation.kind}:${relation.target.id}`;if(seen.has(key))return false;seen.add(key);return true;});
   if(!cancelled)setRelationState({key:relationKey,model:partById,items});
  },0);});
  return()=>{cancelled=true;cancelAnimationFrame(frame);if(timer)clearTimeout(timer);};
 },[relationKey,partById,allParts,group]);
 const relations=relationState?.key===relationKey&&relationState.model===partById?relationState.items:[];
 const componentIds=anatomicalComponentIds(choice,relations,partById);
 const relationsPending=!!inspectedParts.length&&!componentIds.length&&(relationState?.key!==relationKey||relationState?.model!==partById);
 const relatedDisabled=relationshipExhausted||relationsPending||(!componentIds.length&&!hasExpandableAnatomicalRelations(relations,partById));
 const relationKinds=(Object.keys(relationLabels) as RelationKind[]).filter(kind=>relations.some(relation=>relation.kind===kind));
 const relationGroups=relationKinds.map(kind=>({kind,items:relations.filter(relation=>relation.kind===kind)}));
 const initialRelationGroup=useMemo(()=>relationKinds[0]?[relationKinds[0]]:[],[relationKinds[0]]);

 return <>
  <div className="detail-header">
   <div className="detail-accent" style={{background:system?.color}}/>
   <div className="structure-title-row"><SheetTitle ref={titleRef} tabIndex={-1} className={`structure-title ${lateralityClass(titleSide)}`}>{group?title:<a className="structure-title-link" href={wikipediaStructure(title)} target="_blank" rel="noopener noreferrer" title={`Read about ${title} on Wikipedia`}>{title}</a>}</SheetTitle></div>
  </div>
  <div className="detail-scroll" key={choice.id}>
   {!!ancestors.length&&<nav className="structure-breadcrumbs" aria-label="Anatomical path">{ancestors.map((parent,index)=>{const display=displayLaterality(parent.name);return <span className="breadcrumb-step" key={`${index}:${parent.id}`}>{index>0&&<span className="breadcrumb-separator" aria-hidden="true">›</span>}<button type="button" className={lateralityClass(display.side)} onClick={()=>onChooseChild(parent)} aria-label={`Select ${parent.name}`} title={`Open ${parent.name}`}>{display.label}</button></span>;})}</nav>}
   {tags.length>0&&<div className="structure-tags">{tags.map(tag=><span className="structure-tag" key={tag}>{tag}</span>)}</div>}
   <SheetDescription className="structure-description">{description}</SheetDescription>
   {coverageLimitation&&<p className="context-note" aria-label="Model coverage">{coverageLimitation}</p>}
   {partOfChoice&&<section className="anatomical-relations" aria-label="Anatomical parent">
    <Accordion multiple defaultValue={['partOf']} className="anatomical-relation-accordion">
     <AccordionItem value="partOf" className="anatomical-relation-row">
      <AccordionTrigger className="anatomical-relation-trigger"><span>Part of</span><span className="anatomical-relation-count">1</span></AccordionTrigger>
      <AccordionContent className="anatomical-relation-panel"><div className="anatomical-relation-links"><span className="anatomical-relation-item"><button type="button" onClick={()=>onChooseChild(partOfChoice)} aria-label={`Select ${partOfChoice.name}`} title={`Open ${partOfChoice.name}`}>{structureName(partOfChoice.name)}</button><small>{partOfChoice.elements.length?`${partOfChoice.elements.length.toLocaleString()} ${partOfChoice.elements.length===1?'piece':'pieces'}`:'No mesh'}</small></span></div></AccordionContent>
     </AccordionItem>
    </Accordion>
   </section>}
   {group&&!!choice.children?.length&&<section className="anatomical-relations" aria-label="Anatomical relationships">
    <h3>Related structures <span className="anatomical-relations-count">{choice.children.length}</span></h3>
    <Accordion key={choice.id} multiple defaultValue={['contains']} className="anatomical-relation-accordion">
     <AccordionItem value="contains" className="anatomical-relation-row">
      <AccordionTrigger className="anatomical-relation-trigger"><span>Contains</span><span className="anatomical-relation-count">{choice.children.length}</span></AccordionTrigger>
      <AccordionContent className="anatomical-relation-panel"><div className="anatomical-relation-links">{choice.children.map(child=>{const display=displayLaterality(child.name);return <span className="anatomical-relation-item" key={child.id}><button type="button" className={lateralityClass(display.side)} onClick={event=>onChooseChild(child,event.ctrlKey||event.metaKey)} aria-label={`Select ${child.name}`} title={`Select ${child.name}; Control-click to add or remove`}>{display.label}</button><small>{child.elements.length?`${child.elements.length.toLocaleString()} ${child.elements.length===1?'piece':'pieces'}`:'No mesh'}</small></span>;})}</div></AccordionContent>
     </AccordionItem>
    </Accordion>
   </section>}
   {selected&&(relationsPending||relationKinds.length>0)&&<section className="anatomical-relations" aria-label="Anatomical relationships" aria-busy={relationsPending}>
    <h3>Related structures <span className="anatomical-relations-count">{relationsPending?'…':relations.length}</span></h3>
    {!relationsPending&&<Accordion key={inspectedKey} multiple defaultValue={initialRelationGroup} className="anatomical-relation-accordion">
     {relationGroups.map(({kind,items})=><AccordionItem value={kind} className="anatomical-relation-row" key={kind}>
      <AccordionTrigger className="anatomical-relation-trigger"><span>{relationLabels[kind]}</span><span className="anatomical-relation-count">{items.length}</span></AccordionTrigger>
      <AccordionContent className="anatomical-relation-panel"><div className="anatomical-relation-links">{items.map(({target,via,viaModeled,note})=>{const fullName=relationName(target),display=displayLaterality(fullName);return <span className="anatomical-relation-item" key={target.id}><button type="button" className={lateralityClass(display.side)} onContextMenu={event=>{if(event.ctrlKey){event.preventDefault();event.shiftKey?onHidePart(target.id):controlChoose(target.id);}}} onClick={event=>event.shiftKey?onHidePart(target.id):event.ctrlKey?controlChoose(target.id):onChoosePart(target.id,event.metaKey)} aria-label={`Select ${fullName}`} title={`Select ${fullName}; Control-click to add or remove; Shift-click to hide`}>{display.label}</button>{note&&<small>{note}</small>}{via?.length&&<small>via {via.join(' → ')}{viaModeled?'':' (not modeled)'}</small>}</span>;})}</div></AccordionContent>
     </AccordionItem>)}
    </Accordion>}
   </section>}
   {anchorParts?.length&&relationshipDepth>0&&<p className="relationship-selection-status">Related structures added · {relationshipDepth} {relationshipDepth===1?'step':'steps'}</p>}
   {!group&&!anchorParts?.length&&new Set(selectedParts.map(part=>part.conceptId)).size>1&&<div className="member-list"><h3>Included structures</h3>{selectedParts.map(part=>{const fullName=structureName(part.name),display=displayLaterality(fullName);return <Button variant="ghost" key={part.id} aria-label={fullName} title={fullName} onContextMenu={event=>{if(event.ctrlKey){event.preventDefault();event.shiftKey?onHidePart(part.id):controlChoose(part.id);}}} onClick={event=>event.shiftKey?onHidePart(part.id):event.ctrlKey?controlChoose(part.id):onChoosePart(part.id,event.metaKey)}><span className={lateralityClass(display.side)}>{display.label}</span></Button>;})}</div>}
   {!group&&hasGeometry&&<section className="inspection-section" aria-label="Peel depth">
    <div className="inspection-subheading"><span>Peel depth</span><output>{depth} / {covering.length}</output></div>
    <Slider aria-label="Peel depth" min={0} max={Math.max(1,covering.length)} step={1} value={[depth]} disabled={!covering.length} onValueChange={value=>onDepth(sliderValue(value))}/>
    <div className="peel-buttons">
     <Button variant="outline" disabled={depth===0} onClick={()=>onDepth(depth-1)}>{previous?`Restore ${structureName(previous.name)}`:'Nothing to restore'}</Button>
     <Button variant="outline" disabled={depth>=covering.length} onClick={()=>onDepth(depth+1)}>{next?`Remove ${structureName(next.name)}`:'Nothing to remove'}</Button>
    </div>
   </section>}
   {(terminology.latin||terminology.ta98||terminology.tha||terminology.fma||terminology.ontology||localIdentifier)&&<section className="identifiers-section" aria-label="Identificators"><h3 className="identifiers-heading"><span>Identificators</span></h3><div className="structure-terminology" role="list" aria-label="Anatomical identifiers">{terminology.ta98&&<IdentifierPill kind="TA98" value={terminology.ta98} href={reference('TA98')}/>}{terminology.tha&&<IdentifierPill kind="THA" value={terminology.tha.replace(/^THA:/,'')} href={reference('THA')}/>}{terminology.fma&&<IdentifierPill kind="FMA" value={terminology.fma.replace(/^FMA:/,'')} href={reference('FMA')}/>}{terminology.ontology&&<IdentifierPill kind="UBERON" value={terminology.ontology.replace(/^UBERON:/,'')} href={reference('UBERON')}/>}{localIdentifier&&<IdentifierPill kind={localIdentifier.kind} value={localIdentifier.value} href={reference(localIdentifier.kind)}/>}{terminology.latin&&<IdentifierPill kind="La" value={terminology.latin} language="la" href={reference('La')}/>}</div></section>}
  </div>
  {hasGeometry&&<div className="detail-actions">
   <TooltipProvider delay={150}><div className="detail-action-pill" role="toolbar" aria-label="Selected structure tools">
    <Tooltip><TooltipTrigger render={<button type="button" onClick={onCenter} aria-label="Center selected structure"><Crosshair size={19}/></button>}/><TooltipContent>Center selected structure</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" className={state.isolate?'active':''} onClick={onIsolate} aria-pressed={state.isolate} aria-label={state.isolate?`Show surrounding ${cell?'components':'anatomy'}`:`Show selected ${cell?'component':'structure'} only`}><Focus size={20}/></button>}/><TooltipContent>{state.isolate?`Show surrounding ${cell?'components':'anatomy'}`:`Show selected ${cell?'component':'structure'} only`}</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" className={relationshipDepth>0?'active':''} onClick={onExpandRelationships} disabled={relatedDisabled} aria-disabled={relatedDisabled} aria-label={relationshipExhausted?'No further related structures':relationshipDepth>0?'Add the next related structures':'Show related structures together'}><Network size={20}/></button>}/><TooltipContent>{relationshipExhausted?'No further related structures':relationshipDepth>0?'Add the next related structures':'Show related structures together'}</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" onClick={()=>onShowComponents(componentIds)} disabled={!componentIds.length} aria-label="Show components"><Network size={20} style={{transform:'scaleY(-1)'}}/></button>}/><TooltipContent>{componentIds.length?'Show components':'No mapped components'}</TooltipContent></Tooltip>
    <Tooltip><TooltipTrigger render={<button type="button" onClick={onHide} aria-label={`Hide selected ${cell?'component':'structure'}`}><EyeOff size={20}/></button>}/><TooltipContent>Hide selected {cell?'component':'structure'}</TooltipContent></Tooltip>
   </div></TooltipProvider>
  </div>}
 </>;
}
