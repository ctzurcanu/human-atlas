import {liverMicroscopicDetailUrl} from './liver-completion-controls';
import {cardiacSiteConceptCoverage} from './cardiac-label-sites';
import {carotidSurfaceChoiceDescription} from './carotid-surface-features';
import {spinalCordDisplayCount,spinalCordRelations} from './spinal-cord-presentation';
import {namedAnatomyMembers,namedAnatomyRelations} from './anatomical-item-presentation';
import {ileocaecalReferenceLimitation,ileocaecalReferenceDescription} from './ileocaecal-reference-relations';
import {colonicTaeniaSourceLimitation} from './colonic-taenia';
import {endocrineSourceLimitation,endocrineDescription} from './endocrine-relations';
import {digestiveDuctSourceLimitation,digestiveDuctDescription} from './digestive-duct-relations';
import {laryngealDescription,laryngealSelectionDescription} from './laryngeal-relations';
import {smallBowelSourceLimitation} from './small-bowel-relations';
import {hepaticSourceLimitation,hepaticDescription} from './hepatic-relations';
import {useEffect,useMemo,useRef,useState,type RefObject} from 'react';
import {Crosshair,EyeOff,Focus,Network} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {Accordion,AccordionContent,AccordionItem,AccordionTrigger} from '@/components/ui/accordion';
import {Slider} from '@/components/ui/slider';
import {Tooltip,TooltipContent,TooltipProvider,TooltipTrigger} from '@/components/ui/tooltip';
import {SheetDescription,SheetTitle} from '@/components/ui/sheet';
import {SYSTEMS,structureName,surfaceRole,type Part,type Concept,type SceneState} from './anatomy';
import {anatomicalRelations,anatomicalConceptFunctionalRelations,type RelationKind,type ResolvedRelation} from './anatomical-relations';
import {thyroidWholeFunctionalRelations} from './thyroid-functional-relations';
import {atlasIdentifier,identifierReference,terminologyForConcept,type IdentifierKind} from './anatomical-terminology';
import type {HierarchyChoice} from './hierarchy-choice';
import {ta98PartOfChoice} from './hierarchy-navigation';
import {anatomicalComponentRelations,anatomicalComponentIds,hasExpandableAnatomicalRelations,storedAssemblyParent} from './anatomical-components';
import {ta98ModeledParent} from './ta98-modeled-parent';
import {largeBowelWholeSource} from './large-bowel-frame';
import {gallbladderWholeSource} from './gallbladder-wall';
import {displayLaterality,lateralityClass,structureLabelName} from './laterality';
import {mergeRelationshipNotes} from './relationship-summary';
import {nerveSourceCoverageLimitation} from './nerve-source-review';
import {handVesselCoverageLimitation} from './hand-vascular-relations';
import {urinarySourceCoverageLimitation} from './urinary-relations';
import {coronarySourceLimitation} from './coronary-relations';
import {cardiacSourceLimitation} from './cardiac-relations';
import {pulmonarySourceLimitation} from './pulmonary-relations';
import {cranialSourceLimitation} from './cranial-relations';
import {localDescription,wikipediaDescription,type StructureDescription} from './structure-description';
import {itemDescription} from '../shared/item-description.mjs';

import type {WallPeelStep} from './bowel-wall-peeling';
import {intestinalVillusDetailUrl,intestinalCryptDetailUrl,intestinalLymphoidDetailUrl} from './intestine-detail-navigation';

interface Props {
 titleRef:RefObject<HTMLHeadingElement|null>;
 choice:HierarchyChoice;
 componentChoices?:HierarchyChoice[];
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
 wallSteps?:WallPeelStep[];
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
const relationLabels:Record<RelationKind,string>={before:'Before',after:'After',innervation:'Innervated by',arterial:'Arterial supply',venous:'Venous drainage',innervates:'Innervates',supplies:'Supplies',drains:'Drains',articulates:'Articulates with',connects:'Connects',connectedBy:'Ligaments / tendons',joint:'At this joint',continuous:'Continuous with',covers:'Covers',coveredBy:'Covered by',cartilages:'Cartilages',bones:'Bones',tendons:'Tendons',fascia:'Fascia',muscles:'Muscles',origin:'Origin on',insertion:'Inserts on',partOf:'Part of',contains:'Contains',originFor:'Origin marker for',insertionFor:'Insertion marker for',originSites:'Origin markers',insertionSites:'Insertion markers',counterpart:'Opposite side',adjacent:'Adjacent structures',lymphaticDrainage:'Lymph drains into',lymphaticTributaries:'Receives lymph from'};
const relationName=(part:Part)=>`${structureName(part.name)}${surfaceRole(part.name)?` · ${surfaceRole(part.name)}`:''}`;
const relationCache=new WeakMap<Map<string,Part>,Map<string,ResolvedRelation[]>>();
function IdentifierPill({kind,value,href,language}:{kind:IdentifierKind;value:string;href:string;language?:string}){
 return <span role="listitem" className="identifier-item"><a className={`identifier-pill identifier-${kind.toLowerCase()}`} href={href} target="_blank" rel="noopener noreferrer" title={`Open ${kind} ${value}`}><strong>{kind}</strong><span lang={language}>{value}</span></a></span>;
}

export default function SelectionInspector({titleRef,choice,componentChoices,ancestors,selectedParts,anchorParts,relationshipDepth,relationshipExhausted,partById,concepts,scope,state,viewUrl,covering,wallSteps,depth,onDepth,onCenter,onIsolate,onExpandRelationships,onShowComponents,onHide,onHidePart,onChoosePart,onChooseChild}:Props){
 const title=structureLabelName(choice.name);
 const liverMicroDetail=import.meta.env.DEV?liverMicroscopicDetailUrl(choice):null;
 const microscopicDetail=import.meta.env.DEV?intestinalVillusDetailUrl(choice):null;
 const cryptDetail=import.meta.env.DEV?intestinalCryptDetailUrl(choice):null;
 const lymphoidDetail=import.meta.env.DEV?intestinalLymphoidDetailUrl(choice):null;
 const group=!!choice.children;
 const displayedComponents=componentChoices??choice.children;
 const spinalComponents=componentChoices!==undefined&&spinalCordDisplayCount(selectedParts)!==undefined;
 const selectionSet=choice.id==='selection-set';
 const hasGeometry=selectedParts.length>0;
 const titleSide=displayLaterality(title).side;
 const cell=scope==='cell';
 const inspectedParts=anchorParts?.length?anchorParts:selectedParts;
 const inspectedKey=inspectedParts.map(part=>part.id).join('|');
 const relationKey=`${choice.id}|${choice.name}|${!!anchorParts?.length}|${inspectedKey}`;

 const allParts=useMemo(()=>[...partById.values()],[partById]);
 const selected=inspectedParts[0],systems=[...new Set(inspectedParts.map(part=>part.system))].map(id=>SYSTEMS.find(item=>item.id===id)).filter(item=>!!item),system=systems[0];
 const conceptById=useMemo(()=>new Map(concepts.map(concept=>[concept.id,concept])),[concepts]);
 const includedChoices=useMemo(()=>namedAnatomyMembers(selectedParts,concepts),[selectedParts,concepts]);
 const memberCount=(ids:string[])=>new Set(ids.map(id=>partById.get(id)?.conceptId).filter(Boolean)).size;
 const liverSurfaceReference=concepts.find(c=>/^TA98LIVER:male:(surface|landmark):/.test(c.id)&&(c.id===choice.id||c.ta98Term===(choice.ta98Term??choice.terminology?.ta98))&&c.elements.length===choice.elements.length&&c.elements.every(id=>choice.elements.includes(id)));
 const laryngealSelection=laryngealSelectionDescription(inspectedParts,allParts);
 const coverageLimitation=itemDescription(cardiacSiteConceptCoverage(concepts,choice)??choice.coverageLimitation??liverSurfaceReference?.description??[...new Set(inspectedParts.map(part=>part.coverageLimitation??conceptById.get(part.conceptId)?.coverageLimitation??nerveSourceCoverageLimitation(part)??handVesselCoverageLimitation(part)??urinarySourceCoverageLimitation(part)??coronarySourceLimitation(part)??cardiacSourceLimitation(part)??pulmonarySourceLimitation(part)??cranialSourceLimitation(part)??hepaticSourceLimitation(part)??digestiveDuctSourceLimitation(part)??ileocaecalReferenceLimitation(part)??colonicTaeniaSourceLimitation(part)??smallBowelSourceLimitation(part)??endocrineSourceLimitation(part)).map(note=>itemDescription(note)).map(note=>inspectedParts.length>1&&note?.startsWith('Source identity established for the vessel member')?note.replace(/completeness of "[^"]+" is unverified/, 'completeness of the corresponding TA98 collection is unverified'):note).filter(Boolean))].join(' '));
 const identity=choice.children?choice:conceptById.get(choice.id)??choice;
 const choiceTerm=terminologyForConcept(identity.id,identity.ta98Term,identity.ta98Kind);
 const suppliedTerm=choice.terminology;
 const terminology=identity.ta98Kind==='parent'?choiceTerm:suppliedTerm&&(suppliedTerm.ta98||suppliedTerm.fma||suppliedTerm.ontology)?suppliedTerm:(choiceTerm.ta98||choiceTerm.fma||choiceTerm.ontology||choice.elements.length!==1?choiceTerm:terminologyForConcept(selected?.conceptId??choice.id));
 const atlasId=atlasIdentifier(choice.id);
 const localIdentifier=atlasId?{kind:atlasId.slice(0,atlasId.indexOf(':')) as 'Atlas'|'HA-G',value:atlasId.slice(atlasId.indexOf(':')+1)}:null;
 const tags=group||laryngealSelection?[]:[...new Set(inspectedParts.flatMap(part=>[part.tissue,surfaceRole(part.name),...(part.regions??[])]).filter((tag):tag is string=>!!tag))];
 const path=ancestors.map(parent=>({label:parent.name,search:parent.name}));
 const partOfChoice=storedAssemblyParent(choice,concepts,partById)??ta98PartOfChoice(choice,ancestors)??ta98ModeledParent(choice,concepts,partById);
 const partOfDisplay=displayLaterality(structureLabelName(partOfChoice?.name??''));
 const peelCount=wallSteps?.length??covering.length;
 const next=partById.get(covering[depth]),previous=partById.get(covering[depth-1]);
 const reference=(kind:IdentifierKind)=>identifierReference(kind,terminology,viewUrl);
 const [remote,setRemote]=useState<{title:string;value:StructureDescription}|null>(null);
 const lastControlToggle=useRef<{id:string;at:number}|null>(null);
 const controlChoose=(id:string)=>{const now=performance.now(),last=lastControlToggle.current;if(last?.id===id&&now-last.at<150)return;lastControlToggle.current={id,at:now};onChoosePart(id,true);};
 const reviewedDescription=selected&&(laryngealSelection??laryngealDescription(selected)??ileocaecalReferenceDescription(selected)??digestiveDuctDescription(selected)??hepaticDescription(selected)??endocrineDescription(selected));
 const describedAssembly=conceptById.get(choice.id);
 const assemblyDescription=!group&&!selectionSet&&inspectedParts.length>1&&describedAssembly?.description&&describedAssembly.elements.length===inspectedParts.length&&describedAssembly.elements.every(id=>inspectedParts.some(part=>part.id===id))?describedAssembly.description:undefined;
 const referenceDescription=carotidSurfaceChoiceDescription(choice)??liverSurfaceReference?.description;
 const descriptionPart=selected&&referenceDescription?{...selected,description:referenceDescription}:selected&&assemblyDescription?{...selected,description:assemblyDescription}:selected;
 const hasItemDescription=!!itemDescription(descriptionPart?.description);
 useEffect(()=>{
  if(group||selectionSet||hasItemDescription||reviewedDescription)return;
  let active=true;
  wikipediaDescription(title).then(value=>{if(active&&value)setRemote({title,value});});
  return()=>{active=false;};
 },[title,group,selectionSet,hasItemDescription,reviewedDescription?.text]);
 const local=group||selectionSet?null:localDescription(title,descriptionPart,path);
 const resolved=hasItemDescription?local:reviewedDescription??(remote?.title===title?remote.value:local);
 const colorTerritory=import.meta.env.DEV&&/^guest:dermatomes-myotomes:(?:DERMATOME(?::|-SELECTION$)|TRIGEMINAL:)/.test(choice.id);
 const description=colorTerritory?choice.id.includes(':C1:')?'C1 has no cutaneous dermatome.':choice.id.includes(':Co1:')?'No separate Co1 color territory is present in the source map.':'Approximate sensory territory, highlighted on the surface reference.':!hasGeometry?'No selectable mesh in this model.':selectionSet?`${includedChoices.length.toLocaleString()} selected structures. Choose a structure to inspect its description.`:group&&spinalComponents?`${componentChoices!.length} named anatomical components.`:group?`${displayedComponents?.length??0} named components in this group.`:selected?itemDescription(resolved?.text):'';
 const [relationState,setRelationState]=useState<{key:string;model:Map<string,Part>;items:ResolvedRelation[]}|null>(null);
 useEffect(()=>{
  // Mapped child choices establish Contains links without expanding every
  // leaf of a large hierarchy group. Otherwise resolve the available organ
  // links, including single-mesh organs selected through a TA98 group node.
  if(group&&anatomicalComponentIds(choice,[],partById,concepts).length&&!largeBowelWholeSource(inspectedParts,allParts)&&!gallbladderWholeSource(inspectedParts,allParts)&&!thyroidWholeFunctionalRelations(choice,inspectedParts,allParts)?.length)return;
  let cancelled=false,timer:ReturnType<typeof setTimeout>|undefined;
  // Give the sheet a paint before building the relationship graph. The first
  // skeletal selection can otherwise hold up the entire details panel.
  const frame=requestAnimationFrame(()=>{timer=setTimeout(()=>{
   if(cancelled)return;
   let cache=relationCache.get(partById);
   if(!cache){cache=new Map();relationCache.set(partById,cache);}
   const items=mergeRelationshipNotes([...anatomicalConceptFunctionalRelations(choice,inspectedParts,allParts),...inspectedParts.flatMap(part=>{
    let relations=cache.get(part.id);
    if(!relations){relations=anatomicalRelations(part,allParts);cache.set(part.id,relations);}
    return relations;
   })]);
   if(!cancelled)setRelationState({key:relationKey,model:partById,items});
  },0);});
  return()=>{cancelled=true;cancelAnimationFrame(frame);if(timer)clearTimeout(timer);};
 },[relationKey,partById,allParts,group,concepts]);
 const rawRelations=relationState?.key===relationKey&&relationState.model===partById?relationState.items:[];
 const relations=namedAnatomyRelations(componentChoices?(spinalComponents?spinalCordRelations(rawRelations.filter(r=>r.kind!=='contains'&&(!['continuous','adjacent'].includes(r.kind)||!choice.elements.includes(r.target.id))),concepts):rawRelations.filter(r=>r.kind!=='contains')):anatomicalComponentRelations(choice,rawRelations,partById,concepts),concepts);
 const componentIds=componentChoices?[...new Set(componentChoices.flatMap(c=>c.elements))]:anatomicalComponentIds(choice,relations,partById,concepts);
 const relationsPending=!!inspectedParts.length&&!componentIds.length&&(relationState?.key!==relationKey||relationState?.model!==partById);
 const relatedDisabled=relationshipExhausted||relationsPending||!hasExpandableAnatomicalRelations(relations,partById);
 const relationKinds=(Object.keys(relationLabels) as RelationKind[]).filter(kind=>(kind!=='contains'||!componentChoices)&&relations.some(relation=>relation.kind===kind));
 const relationGroups=relationKinds.map(kind=>({kind,items:relations.filter(relation=>relation.kind===kind)}));
 // The first group initializes the accordion once per semantic selection.
 // Asynchronous supply links must not mutate its uncontrolled default value.
 const relationDefault=useRef<{key:string;value:RelationKind[]}>({key:'',value:[]});
 if(relationDefault.current.key!==relationKey)relationDefault.current={key:relationKey,value:relationKinds[0]?[relationKinds[0]]:[]};
 if(!relationDefault.current.value.length&&relationKinds[0])relationDefault.current.value=[relationKinds[0]];
 const initialRelationGroup=relationDefault.current.value;

 return <>
  <div className="detail-header">
   <div className="detail-accent" style={{background:system?.color}}/>
   <div className="structure-title-row"><SheetTitle ref={titleRef} tabIndex={-1} className={`structure-title ${lateralityClass(titleSide)}`}>{group||selectionSet?title:<a className="structure-title-link" href={wikipediaStructure(title)} target="_blank" rel="noopener noreferrer" title={`Read about ${title} on Wikipedia`}>{title}</a>}</SheetTitle></div>
  </div>
  <div className="detail-scroll" key={choice.id}>
   {!!ancestors.length&&<nav className="structure-breadcrumbs" aria-label="Anatomical path">{ancestors.map((parent,index)=>{const display=displayLaterality(parent.name);return <span className="breadcrumb-step" key={`${index}:${parent.id}`}>{index>0&&<span className="breadcrumb-separator" aria-hidden="true">›</span>}<button type="button" className={lateralityClass(display.side)} onClick={()=>onChooseChild(parent)} aria-label={`Select ${parent.name}`} title={`Open ${parent.name}`}>{display.label}</button></span>;})}</nav>}
   {tags.length>0&&<div className="structure-tags">{tags.map(tag=>{const display=displayLaterality(tag);return <span className={`structure-tag${display.side?` tree-name ${lateralityClass(display.side)}`:''}`} key={tag} title={tag}>{display.label}</span>;})}</div>}
   <SheetDescription className="structure-description">{description}</SheetDescription>
   {coverageLimitation&&coverageLimitation!==description&&(selectionSet
    ?<details className="context-note"><summary>Model coverage</summary><p>{coverageLimitation}</p></details>
    :<p className="context-note" aria-label="Model coverage">{coverageLimitation}</p>)}
   {liverMicroDetail&&<p className="context-note"><a href={liverMicroDetail} target="_blank" rel="noopener noreferrer">Open magnified liver unit</a></p>}
   {microscopicDetail&&<p className="context-note"><a href={microscopicDetail} target="_blank" rel="noopener noreferrer">Open magnified villus detail</a><br/>Separate illustrative fragment; not registered to this organ.</p>}
   {cryptDetail&&<p className="context-note"><a href={cryptDetail} target="_blank" rel="noopener noreferrer">Open magnified intestinal gland detail</a><br/>Separate illustrative crypt; not registered to this organ.</p>}
   {lymphoidDetail&&<p className="context-note"><a href={lymphoidDetail} target="_blank" rel="noopener noreferrer">Open magnified lymphoid follicle detail</a><br/>Separate secondary-follicle illustration; organ placement and surrounding patch tissues remain missing.</p>}
   {partOfChoice&&<section className="anatomical-relations" aria-label="Anatomical parent">
    <Accordion multiple defaultValue={['partOf']} className="anatomical-relation-accordion">
     <AccordionItem value="partOf" className="anatomical-relation-row">
      <AccordionTrigger className="anatomical-relation-trigger"><span>Part of</span><span className="anatomical-relation-count">1</span></AccordionTrigger>
      <AccordionContent className="anatomical-relation-panel"><div className="anatomical-relation-links"><span className="anatomical-relation-item"><button type="button" className={lateralityClass(partOfDisplay.side)} onClick={()=>onChooseChild(partOfChoice)} aria-label={`Select ${partOfChoice.name}`} title={`Open ${partOfChoice.name}`}>{partOfDisplay.label}</button><small>{spinalComponents?'':partOfChoice.elements.length?`${memberCount(partOfChoice.elements)} named structures`:'No mesh'}</small></span></div></AccordionContent>
     </AccordionItem>
    </Accordion>
   </section>}
   {(group||componentChoices)&&!!displayedComponents?.length&&<section className="anatomical-relations" aria-label="Anatomical relationships">
    <h3>{colorTerritory?'Selected territories':'Related structures'} <span className="anatomical-relations-count">{displayedComponents.length}</span></h3>
    <Accordion key={choice.id} multiple defaultValue={['contains']} className="anatomical-relation-accordion">
     <AccordionItem value="contains" className="anatomical-relation-row">
      <AccordionTrigger className="anatomical-relation-trigger"><span>{colorTerritory?'Selection':'Contains'}</span><span className="anatomical-relation-count">{displayedComponents.length}</span></AccordionTrigger>
      <AccordionContent className="anatomical-relation-panel"><div className="anatomical-relation-links">{displayedComponents.map(child=>{const display=displayLaterality(child.name);return <span className="anatomical-relation-item" key={child.id}><button type="button" className={lateralityClass(display.side)} onClick={event=>onChooseChild(child,event.ctrlKey||event.metaKey)} aria-label={`Select ${child.name}`} title={`Select ${child.name}; Control-click to add or remove`}>{display.label}</button><small>{spinalComponents?'':colorTerritory?'Surface territory':cardiacSiteConceptCoverage(concepts,child)?'Provisional reference':child.elements.length?`${memberCount(child.elements)} named structures`:'No mesh'}</small></span>;})}</div></AccordionContent>
     </AccordionItem>
    </Accordion>
   </section>}
   {selected&&(relationsPending||relationKinds.length>0)&&<section className="anatomical-relations" aria-label="Anatomical relationships" aria-busy={relationsPending}>
    <h3>Related structures <span className="anatomical-relations-count">{relationsPending?'…':relationGroups.reduce((n,g)=>n+g.items.length,0)}</span></h3>
    {!relationsPending&&<Accordion key={relationKey} multiple defaultValue={initialRelationGroup} className="anatomical-relation-accordion">
     {relationGroups.map(({kind,items})=><AccordionItem value={kind} className="anatomical-relation-row" key={kind}>
      <AccordionTrigger className="anatomical-relation-trigger"><span>{relationLabels[kind]}</span><span className="anatomical-relation-count">{items.length}</span></AccordionTrigger>
      <AccordionContent className="anatomical-relation-panel"><div className="anatomical-relation-links">{items.map(({target,targetConcept,via,viaModeled,viaRoutes,note})=>{const fullName=targetConcept?.name??relationName(target),display=displayLaterality(fullName);return <span className="anatomical-relation-item" key={target.id}><button type="button" className={lateralityClass(display.side)} onContextMenu={event=>{if(event.ctrlKey){event.preventDefault();event.shiftKey?(targetConcept?.elements??[target.id]).forEach(onHidePart):targetConcept?onChooseChild(targetConcept,true):controlChoose(target.id);}}} onClick={event=>event.shiftKey?(targetConcept?.elements??[target.id]).forEach(onHidePart):targetConcept?onChooseChild(targetConcept,event.ctrlKey||event.metaKey):event.ctrlKey?controlChoose(target.id):onChoosePart(target.id,event.metaKey)} aria-label={`Select ${fullName}`} title={`Select ${fullName}; Control-click to add or remove; Shift-click to hide`}>{display.label}</button>{note&&<small>{note}</small>}{(viaRoutes??(via?.length?[{via,modeled:!!viaModeled}]:[])).map((route,index)=><small key={index}>via {route.via.join(' → ')}{route.modeled?'':' (not modeled)'}</small>)}</span>;})}</div></AccordionContent>
     </AccordionItem>)}
    </Accordion>}
   </section>}
   {anchorParts?.length&&relationshipDepth>0&&<p className="relationship-selection-status">Related structures added · {relationshipDepth} {relationshipDepth===1?'step':'steps'}</p>}
   {!group&&!componentChoices&&!anchorParts?.length&&includedChoices.length>1&&<div className="member-list"><h3>Included structures</h3>{includedChoices.map(child=>{const fullName=structureLabelName(child.name),display=displayLaterality(fullName);return <Button variant="ghost" key={child.id} aria-label={fullName} title={fullName} onContextMenu={event=>{if(event.ctrlKey){event.preventDefault();event.shiftKey?child.elements.forEach(id=>onHidePart(id)):onChooseChild(child,true);}}} onClick={event=>event.shiftKey?child.elements.forEach(id=>onHidePart(id)):onChooseChild(child,event.ctrlKey||event.metaKey)}><span className={lateralityClass(display.side)}>{display.label}</span></Button>;})}</div>}
   {!liverSurfaceReference&&(!group||wallSteps)&&hasGeometry&&<section className="inspection-section" aria-label="Peel depth">
    <div className="inspection-subheading"><span>Peel depth</span><output>{depth} / {peelCount}</output></div>
    <Slider aria-label="Peel depth" min={0} max={Math.max(1,peelCount)} step={1} value={[depth]} disabled={!peelCount} onValueChange={value=>onDepth(sliderValue(value))}/>
    <div className="peel-buttons">
     <Button variant="outline" disabled={depth===0} onClick={()=>onDepth(depth-1)}>{wallSteps?.[depth-1]?`Restore ${structureName(wallSteps[depth-1].name)}`:previous?`Restore ${structureName(previous.name)}`:'Nothing to restore'}</Button>
     <Button variant="outline" disabled={depth>=peelCount} onClick={()=>onDepth(depth+1)}>{wallSteps?(wallSteps[depth]?`Remove ${structureName(wallSteps[depth].name)}`:'Nothing to remove'):next?`Remove ${structureName(next.name)}`:'Nothing to remove'}</Button>
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
