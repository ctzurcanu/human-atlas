import {defaultModelCamera} from './model-defaults';
import {withReviewedCatalogueIdentities} from './reviewed-catalogue-identities';
import {assetUrl} from './asset-url';
import {fetchAssetJson} from './asset-cache';
import AdvancedTools from './advanced-tools';
import {useAtlasConnection} from './use-atlas-connection';
import {ConnectSessionButton} from './connect-tools';
import {expandedControls,applyExpandedControls} from './connect-dom';
import {validCamera} from '../shared/camera-frame.mjs';
import type {AdvancedPresentation,ConnectSnapshot,ScenePose} from './connect-state';
import SectionTools from './section-tools';
import {sectionStack} from './section-stack';
import {useViewHistory} from './use-view-history';
import {readViewUrl,viewUrl} from './viewer-state';
import {compactViewUrl,viewParameters} from './view-url';
import {explosionMarks,explosionToSlider,sliderToExplosion} from './transition-motion';
import {readEmbedUi,type EmbedUi} from './embed';
import ShareView from './share-view';
import ExportTools,{type SceneCapture} from './export-tools';
import type {SceneFrameCapture} from './page-capture';
import SelectionInspector from './selection-inspector';
import {requestedSelectionChoice} from './requested-selection-choice';
import {anatomicalRelations,anatomicalConceptFunctionalRelations} from './anatomical-relations';
import {expandRelationshipSelection,type RelationshipExpansion} from './relationship-selection';
import {displayLaterality,lateralityClass} from './laterality';
import SystemTree,{type Hierarchy} from './system-tree';
import DepthTree from './depth-tree';
import GuestTree from './guest-tree';
import type {Bounds3} from './section-plane';
import {anatomicalComponentRelations} from './anatomical-components';
import DermatomeSurface from './dermatome-surface';
import {dermatomeTerritoriesFor,updateDermatomeSelection,dermatomeChoice} from './dermatome-colors';
import GuestHierarchyMenu from './guest-hierarchy-menu';
import {useGuestHierarchies} from './use-guest-hierarchies';
import {DEPTH_LAYERS} from './depth-layers';
import {depthPosition,setDepthPosition} from './depth-control';
import {depthLayerFor} from './depth-layers';
import type {HierarchyChoice} from './hierarchy-choice';
import {hierarchyAncestors,hierarchyNavigation} from './hierarchy-navigation';
import {terminologyForConcept} from './anatomical-terminology';
import {flushSync} from 'react-dom';
import {registerAtlasTools,type AgentActions,type AgentViewChange} from './agent-tools';
import {useEffect,useMemo,useRef,useState} from 'react';
import {ArrowUpRight,BookmarkPlus,Layers3,Pause,RotateCcw,RotateCw,ScanLine,Search,SlidersHorizontal,X} from 'lucide-react';
import {Select,SelectTrigger,SelectValue,SelectContent,SelectItem} from '@/components/ui/select';
import {Button} from '@/components/ui/button';
import {Slider} from '@/components/ui/slider';
import {Sheet,SheetContent,SheetTitle,SheetDescription} from '@/components/ui/sheet';
import {Combobox,ComboboxInput,ComboboxContent,ComboboxList,ComboboxItem,ComboboxEmpty} from '@/components/ui/combobox';
import AnatomyScene from './scene';
import {CELL_VISIBLE,DEFAULT_VISIBLE,SYSTEMS,structureName,type Atlas,type Part,type Concept,type SceneState,type View} from './anatomy';
const urlChoice=(atlas:Atlas,state:SceneState,terms:string[])=>{
 if(state.dermatomeSelected?.length)return dermatomeChoice(state.dermatomeSelected);
 if(terms.length===1)return atlas.concepts.find(concept=>concept.id===terms[0])??null;
 if(terms.length||!state.selected.length)return null;
 const selected=new Set(state.selected);
 return atlas.concepts.find(concept=>concept.elements.length===selected.size&&concept.elements.every(id=>selected.has(id)))??null;
};
const initial:SceneState={focus:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:true,explode:0,visible:DEFAULT_VISIBLE,selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0};
type Model='male-full'|'male-detail'|'male'|'female'|'embryo'|'cell'|'local-male'|'local-female'|'local-reference'|'local-ta98'|'local-female-ta98';
type LayerHierarchy=Hierarchy|'depth'|`guest:${string}`;
const layerHierarchy=(tree:string|null,model:Model):LayerHierarchy=>tree==='regions'||tree==='depth'&&model!=='cell'?tree:tree?.startsWith('guest:')?tree as `guest:${string}`:'systems';
const modelChoices:{value:Model;label:string}[]=[
 {value:'male-full',label:'Male · full resolution + skin'},
 {value:'male-detail',label:'Male · detailed'},
 {value:'male',label:'Male · standard'},
 {value:'female',label:'Female · partial'},
 {value:'embryo',label:'Embryo 3month'},
 {value:'cell',label:'Human eukaryotic cell'},
 ...(import.meta.env.DEV?[{value:'local-ta98' as const,label:'Male · TA98 combined review'},{value:'local-female-ta98' as const,label:'Female · TA98 review'},{value:'local-reference' as const,label:'Male · reference GLBs'},{value:'local-male' as const,label:'Male · local Zygote'},{value:'local-female' as const,label:'Female · local Zygote'}]:[]),
];
const parseModel=(value:string|null):Model=>modelChoices.find(choice=>choice.value===value)?.value??'male-detail';
const cataloguePath=(model:Model)=>model==='local-female-ta98'?(new URLSearchParams(window.location.search).get('delivery')==='0'?'/local-models/ta98-female.json':'/local-models/ta98-female-runtime.json'):model==='local-ta98'?(new URLSearchParams(window.location.search).get('delivery')==='1'?'/local-models/ta98-runtime.json':'/local-models/ta98-male.json'):model==='local-reference'?'/local-models/reference.json':model==='local-male'?'/local-models/male.json':model==='local-female'?'/local-models/female.json':model==='male-full'||model==='male-detail'?'/models/atlas-male-complete.json':model==='female'?'/models/atlas-hra-female.json':model==='embryo'?(import.meta.env.VITE_EMBRYO_CATALOGUE||'/models/embryo-3month/atlas.json'):model==='cell'?'/models/atlas-cell.json':'/models/atlas.json';
const fetchCatalogue=(url:string,init:RequestInit,message:string)=>fetchAssetJson<Atlas>(url,init,message).catch(error=>{
 if(import.meta.env.DEV&&url.endsWith('/local-models/embryo-cs23/atlas.json')&&!init.signal?.aborted)return fetchAssetJson<Atlas>(assetUrl('/models/atlas-embryo.json'),init,message);
 throw error;
}).then(withReviewedCatalogueIdentities);
const defaultLayers=(model:Model,sex:'male'|'female')=>model==='embryo'?[...DEFAULT_VISIBLE,'integumentary' as const,'pregnancy' as const]:model==='cell'?CELL_VISIBLE:sex==='female'||model==='male-full'||model.startsWith('local-')?[...DEFAULT_VISIBLE,'integumentary' as const]:DEFAULT_VISIBLE;
export default function Home(){
 const detailTitle=useRef<HTMLHeadingElement>(null);
 const studioRef=useRef<HTMLElement>(null);
 const framed=window.self!==window.top;
 const [compactFrame,setCompactFrame]=useState(()=>framed&&(window.innerWidth<768||window.innerHeight<600));
 useEffect(()=>{if(!framed)return;const update=()=>setCompactFrame(window.innerWidth<768||window.innerHeight<600);window.addEventListener('resize',update);return()=>window.removeEventListener('resize',update);},[framed]);
 const embedParams=viewParameters(location.href);
 const embedded=embedParams.get('embed')==='1'||compactFrame;
 const embedUi=embedded?readEmbedUi(embedParams.toString()):null;
 const requestedUi=embedParams.get('ui')?.split(',')??[];
 const show=(control:EmbedUi)=>{
  if(control==='camera')return requestedUi.includes('camera')&&(!embedded||!!embedUi?.has('camera'));
  if(compactFrame&&control==='explode')return requestedUi.includes(control);
  return !embedded||!!embedUi?.has(control);
 };
 const [model,setModel]=useState<Model>(()=>parseModel(viewParameters(location.href).get('model')));
 const [hierarchy,setHierarchy]=useState<LayerHierarchy>(()=>layerHierarchy(viewParameters(location.href).get('tree'),model));
 const guest=useGuestHierarchies();
 const guestDefinition=guest.hierarchies.find(item=>hierarchy===`guest:${item.id}`);
 const [explodeSteps,setExplodeSteps]=useState(7),[layersVisible,setLayersVisible]=useState(true),[mobileLayersOpen,setMobileLayersOpen]=useState(false);
 const sex=model==='female'||model==='embryo'||(model==='local-female'||model==='local-female-ta98')?'female':'male';
 const [atlas,setAtlas]=useState<Atlas|null>(null),[progress,setProgress]=useState(0),[error,setError]=useState(''),[panel,setPanel]=useState<'search'|'advanced'|'sections'|null>(()=>embedParams.has('connect')?'advanced':null),[details,setDetails]=useState(false),[about,setAbout]=useState(false),[query,setQuery]=useState(''),[chosen,setChosen]=useState<HierarchyChoice|null>(null);
 const [dermatomeBounds,setDermatomeBounds]=useState<Bounds3|undefined>();
 const [frontPanel,setFrontPanel]=useState<'layers'|'search'|'advanced'|'sections'|'details'>('layers');
 const cameraRef=useRef<number[]|undefined>(undefined),captureRef=useRef<SceneCapture|null>(null),frameCaptureRef=useRef<SceneFrameCapture|null>(null),pendingUrl=useRef(location.href);
 const slideTransitionSerial=useRef(0);
 const slideNavigationPending=useRef(false);
 const [slideTransition,setSlideTransition]=useState<{id:number;from:number[]}|undefined>();
 const loadedModel=useRef<Model|null>(null);
 const agentActionsRef=useRef<AgentActions|null>(null);
 const agentInspectRef=useRef<(concept:Concept)=>void>(()=>{});
 const followingRef=useRef(false),replayingUiRef=useRef(false);
 const {state,setState:updateSceneState,loadState,undo,redo}=useViewHistory(initial,()=>cameraRef.current);
 const selectedGuest=useMemo(()=>{const view=guestDefinition?.views?.find(item=>item.id===state.guestView);return view?{...guestDefinition!,roots:view.roots}:guestDefinition;},[guestDefinition,state.guestView]);
 useEffect(()=>{if(guestDefinition?.id==='genes')for(const id of state.guestExtensions??[])void guest.extend('genes',id).catch(()=>{});},[guestDefinition,state.guestExtensions,guest.extend]);
 const setState=(update:SceneState|((state:SceneState)=>SceneState))=>{if(!followingRef.current)updateSceneState(update);};
 const [addSelection,setAddSelection]=useState(false),[shareLink,setShareLink]=useState(''),[exportStatus,setExportStatus]=useState('');
 const [covering,setCovering]=useState<string[]>([]);
 const [relationshipExpansion,setRelationshipExpansion]=useState<RelationshipExpansion|null>(null);
 const peelBaseline=useRef<string[]>([]);
 const advancedPresentationRef=useRef<AdvancedPresentation|null>(null),scenePoseRef=useRef<ScenePose|null>(null),connectSnapshotRef=useRef<()=>ConnectSnapshot>(()=>{throw new Error('Viewer is still starting.');});
 const [slidesRunning,setSlidesRunning]=useState(false);
 const [advancedToolType,setAdvancedToolType]=useState('connect');
 useEffect(()=>{if(panel!=='advanced')setSlidesRunning(false);},[panel]);
 const connection=useAtlasConnection(()=>connectSnapshotRef.current());followingRef.current=connection.following;
 useEffect(()=>{if(connection.role==='host'&&connection.controlsOpen){setPanel('advanced');setFrontPanel('advanced');}},[connection.role,connection.controlsOpen]);
 useEffect(()=>{
  const remote=connection.remoteSnapshot;if(!connection.following||!remote)return;
  if(remote.model!==model){pendingUrl.current='';setModel(parseModel(remote.model));return;}
  if(!atlas||loadedModel.current!==model)return;
  const parameters=new URLSearchParams();parameters.set('tree',remote.hierarchy);remote.guestSources.forEach(source=>parameters.append('guest',source));guest.restore(parameters.toString());
  setHierarchy(layerHierarchy(remote.hierarchy,model));cameraRef.current=remote.pose?.camera??remote.state.camera;loadState(remote.state);setSlideTransition(undefined);
  setChosen(remote.choice?{id:remote.choice.id,name:remote.choice.name,elements:remote.choice.elements,terminology:remote.choice.terminology,...(remote.choice.group?{children:remote.choice.children?.map(child=>({id:child.id,name:child.name,elements:child.elements,terminology:child.terminology,...(child.group?{children:[]}:{} )}))??[]}:{} )}:null);
  setDetails(remote.ui.details);setPanel(remote.ui.panel);setFrontPanel(remote.ui.frontPanel);setLayersVisible(remote.ui.layersVisible);setMobileLayersOpen(remote.ui.mobileLayersOpen);setQuery(remote.ui.query);setAddSelection(remote.ui.addSelection);setCovering(remote.covering);setRelationshipExpansion(remote.relationshipExpansion??null);
 },[connection.remoteSnapshot,connection.following,model,atlas,loadState,guest.restore]);
 useEffect(()=>{if(loadedModel.current===model&&!pendingUrl.current)return;const abort=new AbortController();setProgress(0);setError('');setAtlas(null);setChosen(null);setRelationshipExpansion(null);setDetails(false);setPanel(current=>current==='advanced'?current:null);setMobileLayersOpen(false);setQuery('');cameraRef.current=undefined;fetchCatalogue(assetUrl(cataloguePath(model)),{signal:abort.signal,...(model.startsWith('local-')||import.meta.env.DEV&&model==='embryo'?{cache:'no-store' as const}:{})},model==='local-female-ta98'?'Female TA98 review is missing. Run npm run build:ta98-female-review.':model==='local-ta98'?'TA98 review is missing. Run npm run build:ta98-male-review.':model==='local-reference'?'Reference GLBs are missing. Run npm run import:reference-glb.':model.startsWith('local-')?'Local models are missing. Run npm run import:local-models.':'The anatomy catalogue could not be loaded.').then(data=>{setAtlas(data as Atlas);const base:SceneState={...initial,camera:defaultModelCamera(model),skinOpacity:model==='cell'?.18:model==='embryo'?1:model==='male-full'||model==='local-reference'?1:model.startsWith('local-')?0:.1,visible:defaultLayers(model,sex)};const restored=pendingUrl.current?readViewUrl(pendingUrl.current,data as Atlas,base):base;const terms=viewParameters(pendingUrl.current).getAll('select');setChosen(urlChoice(data as Atlas,restored,terms));loadState(restored);setDetails(restored.selected.length>0||!!restored.dermatomeSelected?.length);if(restored.sections?.some(section=>section.enabled)){setPanel(current=>current==='advanced'?current:'sections');setFrontPanel(current=>current==='advanced'?current:'sections');}loadedModel.current=model;pendingUrl.current='';}).catch(e=>{if(e.name!=='AbortError')setError(e.message);});return()=>abort.abort();},[model]);
 useEffect(()=>{const key=(e:KeyboardEvent)=>{if(!followingRef.current&&e.key==='/'&&!(e.target instanceof HTMLInputElement)&&!(e.target instanceof HTMLTextAreaElement)){e.preventDefault();setPanel('search');setFrontPanel('search');setDetails(false);}};window.addEventListener('keydown',key);return()=>window.removeEventListener('keydown',key);},[]);
 useEffect(()=>{const pop=()=>{if(!slideNavigationPending.current)setSlideTransition(undefined);slideNavigationPending.current=false;pendingUrl.current=location.href;const params=viewParameters(location.href),next=parseModel(params.get('model'));setHierarchy(layerHierarchy(params.get('tree'),next));guest.restore(viewParameters(location.href).toString());setRelationshipExpansion(null);if(next!==model)setModel(next);else if(atlas){const terms=params.getAll('select'),restored=readViewUrl(location.href,atlas,model==='embryo'?{...initial,camera:defaultModelCamera(model),skinOpacity:1,visible:defaultLayers(model,sex)}:initial);setChosen(urlChoice(atlas,restored,terms));loadState(restored);setDetails(restored.selected.length>0||!!restored.dermatomeSelected?.length);if(restored.sections?.some(section=>section.enabled)){setPanel(current=>current==='advanced'?current:'sections');setFrontPanel(current=>current==='advanced'?current:'sections');}pendingUrl.current='';}};window.addEventListener('popstate',pop);return()=>window.removeEventListener('popstate',pop);},[model,atlas,loadState,guest.restore]);
 const currentViewUrl=()=>{const url=new URL(viewUrl(location.href,model,state,cameraRef.current??state.camera,atlas??undefined));if(hierarchy!=='systems')url.searchParams.set('tree',hierarchy);guest.sources.forEach(source=>url.searchParams.append('guest',source));if(embedParams.has('ui'))url.searchParams.set('ui',embedParams.get('ui')!);return compactViewUrl(url.href);};
 const openSavedView=(url:string,animate:boolean)=>{
  url=compactViewUrl(url);
  const from=cameraRef.current;
  const params=viewParameters(url),target=(params.get('frame')??params.get('camera'))?.split(',').map(Number);
  setSlideTransition(animate&&!window.matchMedia('(prefers-reduced-motion: reduce)').matches&&parseModel(params.get('model'))===model&&from&&validCamera(from)&&target&&validCamera(target)?{id:++slideTransitionSerial.current,from:[...from]}:undefined);
  slideNavigationPending.current=true;
  history.pushState(null,'',url);window.dispatchEvent(new PopStateEvent('popstate'));
 };
 const share=()=>{const url=currentViewUrl();history.replaceState(null,'',url);setShareLink(url);};
 const downloadPng=async()=>{try{setExportStatus('Preparing PNG…');const blob=await captureRef.current?.();if(!blob)throw new Error('The 3D view is still loading.');const url=URL.createObjectURL(blob),link=document.createElement('a');link.href=url;link.download=`human-atlas-${model}-${Date.now()}.png`;document.body.appendChild(link);link.click();link.remove();setExportStatus('PNG download started.');setTimeout(()=>{URL.revokeObjectURL(url);setExportStatus('');},60_000);return {ok:true};}catch(error){const message=error instanceof Error?error.message:'Could not export the image.';setExportStatus(message);return {ok:false,error:message};}};
 const parts=useMemo(()=>new Map(atlas?.parts.map(p=>[p.id,p])),[atlas]);
 const conceptsById=useMemo(()=>new Map(atlas?.concepts.map(concept=>[concept.id,concept])),[atlas]);
 const counts=useMemo(()=>Object.fromEntries(SYSTEMS.map(s=>[s.id,atlas?.parts.filter(p=>p.system===s.id).length??0])),[atlas]);
 const selectionKey=state.selected.join('|');
 const activeExpansion=relationshipExpansion?.selectedKey===selectionKey?relationshipExpansion:null;
 const selectedParts=state.selected.map(id=>parts.get(id)).filter(p=>!!p),selected=selectedParts[0];
 const expansionAnchors=activeExpansion?.anchorIds.map(id=>parts.get(id)).filter(p=>!!p);
 const currentChoice=expansionAnchors?.length?chosen&&chosen.elements.length===expansionAnchors.length&&chosen.elements.every(id=>activeExpansion?.anchorIds.includes(id))?chosen:{id:expansionAnchors[0].conceptId,name:expansionAnchors[0].name,elements:activeExpansion!.anchorIds}:chosen&&chosen.elements.length===state.selected.length&&chosen.elements.every(id=>state.selected.includes(id))?chosen:selectedParts.length===1?{id:selected!.conceptId,name:selected!.name,elements:state.selected}:{id:'selection-set',name:`${state.selected.length} selected pieces`,elements:state.selected};
 const navigation=useMemo(()=>atlas?hierarchyNavigation(atlas,hierarchy,selectedGuest):null,[atlas,hierarchy,selectedGuest]);
 const detailAncestors=useMemo(()=>navigation?hierarchyAncestors(navigation,currentChoice):[],[navigation,currentChoice.id,selectionKey]);
 useEffect(()=>{peelBaseline.current=state.hidden??[];setCovering(connection.following?connection.remoteSnapshot?.covering??[]:[]);},[selectionKey]);
 const peelDepth=covering.filter(id=>(state.hidden??[]).includes(id)).length;
 const restoreInspectionHidden=(hidden:string[]|undefined)=>{const peeled=new Set(covering.filter(id=>!peelBaseline.current.includes(id)));return (hidden??[]).filter(id=>!peeled.has(id));};
 const setPeelDepth=(depth:number)=>setState(s=>({...s,hidden:[...new Set([...restoreInspectionHidden(s.hidden),...covering.slice(0,depth)])]}));
 const clearSelection=()=>{setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden),selected:[],dermatomeSelected:[],isolate:false,camera:undefined}));setChosen(null);setRelationshipExpansion(null);setDetails(false);setCovering([]);};
 const hideSelection=()=>{setState(s=>({...s,hidden:[...new Set([...restoreInspectionHidden(s.hidden),...s.selected])],selected:[],dermatomeSelected:[],isolate:false,camera:undefined}));setChosen(null);setRelationshipExpansion(null);setDetails(false);setCovering([]);};
 const hidePart=(id:string,peel=false)=>{
  if(!parts.has(id))return;
  const remaining=state.selected.filter(selectedId=>selectedId!==id);
  setState(s=>({...s,hidden:[...new Set([...(peel?s.hidden??[]:restoreInspectionHidden(s.hidden)),id])],selected:s.selected.filter(selectedId=>selectedId!==id),isolate:remaining.length?s.isolate:false,camera:peel||remaining.length?s.camera:undefined}));
  if(!state.selected.includes(id))return;
  setRelationshipExpansion(null);
  const next=parts.get(remaining[0]);
  setChosen(next?{id:next.conceptId,name:structureName(next.name),elements:[next.id]}:null);
  if(!next){setDetails(false);setCovering([]);}
 };
 useEffect(()=>{const key=(e:KeyboardEvent)=>{if((e.key!=='Delete'&&e.key!=='Backspace')||e.defaultPrevented||e.metaKey||e.ctrlKey||e.altKey||!state.selected.length)return;const target=e.target;if(target instanceof HTMLElement&&(target.isContentEditable||target.closest('input, textarea, select, [role="textbox"]')))return;e.preventDefault();hideSelection();};window.addEventListener('keydown',key);return()=>window.removeEventListener('keydown',key);},[state.selected.length,covering]);
 const results=useMemo(()=>{if(!atlas)return[];const term=query.toLowerCase().trim(),source=term.replace(/^atlas:/,'');if(!term)return ['heart','brain','liver','stomach','spleen','pancreas','urinary bladder','trachea'].map(name=>atlas.concepts.find(c=>c.name.toLowerCase()===name)).filter((x):x is Concept=>!!x);return atlas.concepts.filter(c=>{const ids=terminologyForConcept(c.id,c.ta98Term,c.ta98Kind);return c.name.toLowerCase().includes(term)||c.id.toLowerCase().includes(source)||[ids.ta98&&`ta98:${ids.ta98}`,ids.tha,ids.fma,ids.ontology,ids.latin&&`la:${ids.latin}`,ids.latin].some(value=>value?.toLowerCase().includes(term));}).sort((a,b)=>a.name.length-b.name.length).slice(0,80);},[atlas,query]);
 const sensorySelectionKey=(state.dermatomeSelected??[]).join(',');
 useEffect(()=>{if(selectedGuest?.id!=='dermatomes-myotomes')return;if(state.dermatomeSelected?.length){setChosen(dermatomeChoice(state.dermatomeSelected));setDetails(true);}else if(chosen?.id.startsWith('guest:dermatomes-myotomes:')&&/(?:DERMATOME|TRIGEMINAL)/.test(chosen.id)){setChosen(null);setDetails(false);}},[sensorySelectionKey]);
 const choose=(c:HierarchyChoice,toggle=false)=>{
  if(followingRef.current)return;
  const sensoryKeys=/(?:DERMATOME|TRIGEMINAL)(?::|$)/.test(c.id)?dermatomeTerritoriesFor(c.id):[];
  if(selectedGuest?.id==='dermatomes-myotomes'&&sensoryKeys.length){
   const keys=updateDermatomeSelection(state.dermatomeSelected??[],sensoryKeys,toggle||addSelection),next=dermatomeChoice(keys);
   setChosen(next);setRelationshipExpansion(null);setState(s=>({...s,dermatomeSelected:keys,guestQuery:keys.length===1?next?.id.replace('guest:dermatomes-myotomes:',''):undefined,selected:[],isolate:false,rotate:false,labels:true}));setDetails(!!next);setFrontPanel('details');setCovering([]);setPanel(current=>current==='sections'||current==='advanced'?current:null);return;
  }
  if(state.dermatomeSelected?.length)setState(s=>({...s,dermatomeSelected:[]}));
  if(!c.elements.length){
   setChosen(c);setRelationshipExpansion(null);
   setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden),selected:[],isolate:false,camera:undefined}));
   setDetails(true);setFrontPanel('details');setCovering([]);
   setPanel(current=>current==='sections'||current==='advanced'?current:null);
   return;
  }
  if(toggle&&c.elements.length&&c.elements.every(id=>state.selected.includes(id))){
   const removed=new Set(c.elements),remaining=state.selected.filter(id=>!removed.has(id)),next=parts.get(remaining[0]);
   setChosen(remaining.length===1&&next?{id:next.conceptId,name:structureName(next.name),elements:[next.id]}:null);
   setRelationshipExpansion(null);
   setState(s=>{const selected=s.selected.filter(id=>!removed.has(id));return {...s,selected,isolate:selected.length?s.isolate:false,camera:selected.length?s.camera:undefined};});
   setDetails(remaining.length>0);
   if(!remaining.length)setCovering([]);
   return;
  }
  setChosen(c);setRelationshipExpansion(null);
  setState(s=>{
   const group=!!c.children,memberIds=new Set(c.elements),members=c.elements.map(id=>parts.get(id)).filter(part=>!!part);
   return {...s,hidden:restoreInspectionHidden(s.hidden).filter(id=>!memberIds.has(id)),depthHidden:s.depth!==undefined?s.depthHidden?.filter(id=>!members.some(part=>depthLayerFor(part)===id)):s.depthHidden,selected:toggle||addSelection?[...new Set([...s.selected,...c.elements])]:c.elements,isolate:false,contextOpacity:s.contextOpacity===1?.18:s.contextOpacity,rotate:false,camera:undefined,
    ...(group?{visible:[...new Set([...s.visible,...members.map(part=>part.system)])],depthHidden:s.depthHidden?.filter(id=>!members.some(part=>depthLayerFor(part)===id)),skinOpacity:members.some(part=>depthLayerFor(part)==='skin')&&!(s.skinOpacity??1)?.1:s.skinOpacity}:{}),
   };
  });
  setDetails(true);setFrontPanel('details');setPanel(current=>current==='sections'||current==='advanced'?current:null);
 };
 agentInspectRef.current=choose;
 useEffect(()=>{if(!atlas)return;return registerAtlasTools(atlas,c=>flushSync(()=>agentInspectRef.current(c)),{
  snapshot:()=>agentActionsRef.current?.snapshot(),update:change=>agentActionsRef.current?.update(change),command:action=>agentActionsRef.current?.command(action),url:()=>agentActionsRef.current?.url()??location.href,openView:url=>agentActionsRef.current?.openView?.(url),
 });},[atlas]);
 const choosePart=(id:string,toggle=false)=>{
  const p=parts.get(id);if(!p)return;
  if(toggle&&state.selected.includes(id)){
   const remaining=state.selected.filter(selectedId=>selectedId!==id),next=parts.get(remaining[0]);
   setChosen(remaining.length===1&&next?{id:next.conceptId,name:structureName(next.name),elements:[next.id]}:null);
   setRelationshipExpansion(null);
   setState(s=>{const selected=s.selected.filter(selectedId=>selectedId!==id);return {...s,selected,isolate:selected.length?s.isolate:false,camera:selected.length?s.camera:undefined};});
   setDetails(remaining.length>0);
   if(!remaining.length)setCovering([]);
   return;
  }
  setChosen({id:p.conceptId,name:structureName(p.name),elements:[id]});setRelationshipExpansion(null);
  setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden).filter(hiddenId=>hiddenId!==id),depthHidden:s.depth!==undefined?s.depthHidden?.filter(layer=>depthLayerFor(p)!==layer):s.depthHidden,selected:toggle||addSelection?[...new Set([...s.selected,id])]:[id],isolate:false,contextOpacity:s.contextOpacity===1?.18:s.contextOpacity,rotate:false,camera:undefined}));
  setDetails(true);setFrontPanel('details');setPanel(current=>current==='sections'||current==='advanced'?current:null);
 };
 const showComponents=(ids:string[])=>{const members=ids.map(id=>parts.get(id)).filter((part):part is Part=>!!part&&!part.suppressed);if(!members.length)return;const selected=members.map(part=>part.id);setRelationshipExpansion(null);setChosen(members.length===1?{...(atlas?.concepts.find(concept=>concept.id===members[0].conceptId)),id:members[0].conceptId,name:structureName(members[0].name),elements:selected}:{...currentChoice,elements:selected});setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden).filter(id=>!selected.includes(id)),selected,isolate:true,labels:true,rotate:false,camera:undefined,visible:[...new Set([...s.visible,...members.map(part=>part.system)])],depthHidden:s.depthHidden?.filter(layer=>!members.some(part=>depthLayerFor(part)===layer))}));};
 const expandRelationships=()=>{if(!state.selected.length||activeExpansion?.exhausted)return;const functional=!activeExpansion?anatomicalConceptFunctionalRelations(currentChoice,selectedParts,atlas?.parts??[]).filter(relation=>!relation.target.suppressed).map(relation=>relation.target.id):[];const result=expandRelationshipSelection(state.selected,activeExpansion,id=>{const part=parts.get(id),choice=part?conceptsById.get(part.conceptId)??{id:part.conceptId,name:structureName(part.name),elements:[id]}:undefined;return [...(part&&choice?anatomicalComponentRelations(choice,anatomicalRelations(part,atlas?.parts??[]),parts,atlas?.concepts??[]).filter(relation=>relation.kind!=='partOf'&&!relation.target.suppressed).map(relation=>relation.target.id):[]),...(id===state.selected[0]?functional:[])];});setRelationshipExpansion(result.expansion);setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden),selected:result.selected,isolate:true,labels:true,rotate:false,camera:undefined}));};
 const reset=()=>{setState(s=>({...initial,camera:defaultModelCamera(model),hidden:atlas?.defaultHidden??[],skinOpacity:model==='cell'?.18:model==='embryo'?1:model==='male-full'||model==='local-reference'?1:model.startsWith('local-')?0:.1,visible:defaultLayers(model,sex),...(dermatomeReference?{depth:0,depthHidden:[],skinOpacity:1}:{}),reset:s.reset+1}));setChosen(null);setRelationshipExpansion(null);setDetails(false);setPanel(null);};
 const openPanel=(next:'search'|'advanced'|'sections')=>{
  if(panel!==next)setFrontPanel(next);
  if(next!=='advanced')setDetails(false);
  if(next==='sections'&&panel!=='sections'){
   setLayersVisible(true);
   if(window.innerWidth<768||window.innerHeight<=600)setMobileLayersOpen(true);
   setState(s=>{const sections=sectionStack(s);sections[0]={...sections[0],enabled:true};return {...s,sections,activeSection:0,section:sections[0],camera:undefined};});
  }
  setPanel(next);
 };
 const toggleLayers=()=>{setFrontPanel('layers');if(embedded||window.innerWidth<768||window.innerHeight<=600){setMobileLayersOpen(value=>!value);return;}setLayersVisible(value=>!value);};
 const setDepth=(value:number,magnetic=false)=>{const depth=sliderToExplosion(value/100,DEPTH_LAYERS.length,magnetic);setState(s=>setDepthPosition(s,depth,atlas?.parts??[]));};
 const setExplosion=(value:number,magnetic=false)=>{const explode=sliderToExplosion(value/100,explodeSteps,magnetic);setState(s=>({...s,explode,rotate:false,camera:undefined}));};
 const resolveAgentPieces=(terms:string[])=>[...new Set(terms.flatMap(term=>{const part=atlas?.parts.find(p=>!p.suppressed&&(p.id===term||p.sourceId===term));if(part)return [part.id];const matches=atlas?.concepts.filter(c=>c.id===term||c.name.toLowerCase()===term.toLowerCase())??[];if(!matches.length)throw new Error(`No structure matches "${term}" in this model.`);return matches.flatMap(c=>c.elements);} ))];
 agentActionsRef.current={
  snapshot:()=>({model,hierarchy,guestHierarchies:guest.hierarchies.map(item=>({id:item.id,name:item.name})),selected:state.selected,systems:state.visible,availableSystems:Object.entries(counts).filter(([,count])=>count>0).map(([id,count])=>({id,pieces:count})),depthHidden:state.depthHidden??[],availableDepthLayers:model==='cell'?[]:DEPTH_LAYERS,hidden:state.hidden??[],view:state.view,region:state.region??'all',context:state.contextOpacity??1,skinOpacity:state.skinOpacity??.1,depth:depthPosition(state),explode:state.explode,labels:state.labels!==false,isolate:state.isolate,rotate:state.rotate,section:state.section??null,sections:sectionStack(state),activeSection:state.activeSection??0,camera:cameraRef.current??state.camera??null,url:currentViewUrl()}),
  url:currentViewUrl,
  openView:url=>openSavedView(url,true),
  update:(change:AgentViewChange)=>{
   if(followingRef.current)throw new Error('This guest view is controlled by its host.');
   if(change.model&&change.model!==model){if(!modelChoices.some(choice=>choice.value===change.model))throw new Error('Unknown model.');setModel(change.model as Model);setHierarchy('systems');return {model:change.model,loading:true};}
   if(change.hierarchy?.startsWith('guest:')&&guest.menuHierarchies.some(item=>`guest:${item.id}`===change.hierarchy)&&!guest.hierarchies.some(item=>`guest:${item.id}`===change.hierarchy))return guest.ensure(change.hierarchy.slice(6)).then(()=>{setHierarchy(change.hierarchy as LayerHierarchy);const rest={...change};delete rest.hierarchy;return Object.keys(rest).length?agentActionsRef.current?.update(rest):{hierarchy:change.hierarchy,loaded:true};});
   if(change.guestUrl)return guest.load(change.guestUrl).then(loaded=>{if(change.hierarchy&&change.hierarchy!==`guest:${loaded.id}`)throw new Error('Guest URL and hierarchy ID disagree.');setHierarchy(`guest:${loaded.id}`);const rest={...change};delete rest.guestUrl;delete rest.hierarchy;if(Object.keys(rest).length)agentActionsRef.current?.update(rest);return {guestUrl:change.guestUrl,hierarchy:`guest:${loaded.id}`,loaded:true};});
   const selected=change.structures===undefined?undefined:resolveAgentPieces(change.structures),hidden=change.hidden===undefined?undefined:resolveAgentPieces(change.hidden);
   if(change.isolate===true&&!(selected??state.selected).length)throw new Error('Isolation requires a selected structure.');
   if(change.rotate===true&&(change.explode??state.explode)>.04)throw new Error('Auto rotation is available only when the model is assembled.');
   if(change.isolate===true&&change.explode!==undefined&&change.explode>0)throw new Error('Isolation and Explode cannot be active together.');
   if(change.activeSection!==undefined&&change.activeSection>=(change.sections?.length??sectionStack(state).length))throw new Error('Add a second section before selecting its tab.');
   if(change.hierarchy){if(change.hierarchy.startsWith('guest:')&&!guest.hierarchies.some(item=>`guest:${item.id}`===change.hierarchy))throw new Error('That guest hierarchy is not loaded.');setHierarchy(change.hierarchy as LayerHierarchy);}
   if(selected){setChosen(atlas?requestedSelectionChoice(atlas,change.structures!,selected):null);setRelationshipExpansion(null);setDetails(selected.length>0);if(selected.length)setFrontPanel('details');setCovering([]);}
   setState(s=>{const sections=change.sections?change.sections.map(item=>item?{enabled:true,...item}:{enabled:false,axis:'axial' as const,position:.38,flip:true}):sectionStack(s),active=Math.min(change.activeSection??s.activeSection??0,sections.length-1);if(change.section!==undefined)sections[active]=change.section?{enabled:true,...change.section}:{enabled:false,axis:'axial',position:.38,flip:true};const depthState=change.depth!==undefined?setDepthPosition(s,change.depth,atlas?.parts??[]):s;return {...depthState,...(selected?{selected}:{}),...(hidden?{hidden}:{}),...(change.systems?{visible:change.systems as SceneState['visible']}:{}),...(change.depthHidden?{depthHidden:change.depthHidden}:{}),...(change.view?{view:change.view as View}:{}),...(change.region?{region:change.region}:{}),...(change.context!==undefined?{contextOpacity:change.context}:selected?.length&&s.contextOpacity===1?{contextOpacity:.18}:{}),...(change.skinOpacity!==undefined?{skinOpacity:change.skinOpacity,depth:undefined}:{}),...(change.explode!==undefined?{explode:change.explode,rotate:change.explode>.04?false:s.rotate}:{}),...(change.labels!==undefined?{labels:change.labels}:{}),...(change.isolate!==undefined?{isolate:change.isolate}:selected?.length===0?{isolate:false}:{}),...(change.rotate!==undefined?{rotate:change.rotate}:{}),...(change.section!==undefined||change.sections!==undefined||change.activeSection!==undefined?{sections,activeSection:active,section:sections[active],camera:change.camera??undefined}:{}),...(change.camera?{camera:change.camera}:{}),...(change.focus?{focus:(s.focus??0)+1}:{}),...(change.view?{reset:s.reset+1}:{} )};});
   return {updated:true,...change};
  },
  command:async action=>{if(followingRef.current&&action!=='download-png')throw new Error('This guest view is controlled by its host.');switch(action){case 'undo':undo();break;case 'redo':redo();break;case 'reset':reset();break;case 'focus':setState(s=>({...s,focus:(s.focus??0)+1,camera:undefined}));break;case 'expand-relationships':expandRelationships();break;case 'hide-selection':hideSelection();break;case 'clear-selection':clearSelection();break;case 'download-png':return downloadPng();default:throw new Error('Unknown viewer action.');}return {action,done:true};},
 };
 connectSnapshotRef.current=()=>({version:1,model,hierarchy,guestSources:[...guest.sources],state:{...state,camera:cameraRef.current??state.camera},choice:chosen?{id:chosen.id,name:chosen.name,elements:chosen.elements,terminology:chosen.terminology,group:!!chosen.children,...(chosen.children?{children:chosen.children.map(child=>({id:child.id,name:child.name,elements:child.elements,terminology:child.terminology,group:!!child.children}))}:{})}:null,ui:{panel,frontPanel,details,layersVisible,mobileLayersOpen,addSelection,query,expanded:expandedControls()},covering,relationshipExpansion,pose:scenePoseRef.current?.model===model?scenePoseRef.current:null,advanced:advancedPresentationRef.current?{...advancedPresentationRef.current,document:{...advancedPresentationRef.current.document,views:advancedPresentationRef.current.document.views.map(({scene,...view})=>view)}}:null});
 useEffect(()=>{connection.publishSnapshot();},[model,hierarchy,guest.sources.join('|'),state,chosen,panel,frontPanel,details,layersVisible,mobileLayersOpen,addSelection,query,covering,relationshipExpansion,connection.publishSnapshot]);
 const sceneState=useMemo(()=>({...state,selectionGroup:!!currentChoice?.children&&!activeExpansion,inspectorOpen:show('details')&&details&&selectedParts.length>0}),[state,!!currentChoice?.children,!!activeExpansion,show('details'),details,selectedParts.length]);
 const dermatomeReference=import.meta.env.DEV&&selectedGuest?.id==='dermatomes-myotomes'&&(!chosen||/(?:^|:)guest:dermatomes-myotomes:/.test(chosen.id))&&!/(?:^|guest:dermatomes-myotomes:)(?:NERVE|MYOTOME|MUSCLE|AXIAL|SEGMENT|UNSEGMENTED|ANATOMY|INNERVATION|CRANIAL|AUTONOMIC):?/.test(chosen?.id??state.guestQuery??'');
 useEffect(()=>{if(dermatomeReference)setState(s=>s.depth===undefined?{...s,depth:0,depthHidden:[],skinOpacity:1}:s);},[dermatomeReference]);
 const blockGuestAction=(event:React.SyntheticEvent)=>{if(connection.following&&!replayingUiRef.current&&event.target instanceof Element&&!event.target.closest('.connect-tools,.connection-session,.download-panel,[aria-label="Advanced tools"],[aria-label="Download"],[aria-label="Stop recording"],[aria-label="Cancel PNG"]')){event.preventDefault();event.stopPropagation();}else if(!connection.following)connection.publishSnapshot();};
 useEffect(()=>{
  const remote=connection.remoteSnapshot;if(!connection.following||!remote||remote.model!==model||!atlas)return;
  let frame=0,attempt=0;const sync=()=>{frame=0;replayingUiRef.current=true;let pending=false;try{pending=applyExpandedControls(remote.ui.expanded??[]);}finally{replayingUiRef.current=false;}if(pending&&++attempt<12)frame=requestAnimationFrame(sync);};
  // Relationships and guest hierarchies can finish loading after the snapshot.
  const watchPanels=()=>studioRef.current?.querySelectorAll('#anatomy-browser,.detail-sheet').forEach(panel=>observer.observe(panel,{childList:true,subtree:true}));
  const observer=new MutationObserver(()=>{watchPanels();if(!frame){attempt=0;frame=requestAnimationFrame(sync);}});if(studioRef.current)observer.observe(studioRef.current,{childList:true});watchPanels();
  frame=requestAnimationFrame(sync);return()=>{cancelAnimationFrame(frame);observer.disconnect();};
 },[connection.remoteSnapshot,connection.following,model,atlas,hierarchy,selectionKey,guest.hierarchies.map(item=>item.id).join('|')]);
 return <main ref={studioRef} className={`studio ${embedded?'embedded':''}`} data-front-panel={frontPanel} data-connect-following={connection.following} data-presenting={slidesRunning||connection.active} onPointerDownCapture={blockGuestAction} onClickCapture={blockGuestAction} onDoubleClickCapture={blockGuestAction} onKeyDownCapture={blockGuestAction}>
  {connection.active&&(panel!=='advanced'||advancedToolType!=='connect')&&<div className="connection-session"><ConnectSessionButton role={connection.role} waiting={connection.role==='guest'&&!connection.following} onOpen={()=>{connection.openControls();setPanel('advanced');setFrontPanel('advanced');}} onStop={connection.stop}/></div>}
  {atlas&&dermatomeReference?<DermatomeSurface onBounds={setDermatomeBounds} sectionTool={show('sections')&&panel==='sections'} onSectionPosition={position=>setState(s=>{const sections=sectionStack(s),active=Math.min(s.activeSection??0,sections.length-1);sections[active]={...sections[active],enabled:true,position};return {...s,sections,activeSection:active,section:sections[active]};})} sex={sex} state={sceneState} selection={chosen?.id??state.guestQuery} onCapture={capture=>{captureRef.current=capture;frameCaptureRef.current=null;}} onTerritory={(id,toggle)=>{choose({id:`guest:dermatomes-myotomes:${id}`,name:`${id.split(":")[1]} ${id.startsWith("TRIGEMINAL")?"facial sensory territory":"dermatome"} (${id.split(":")[2]})`,elements:[]},toggle);}} onProgress={setProgress} onError={setError}/>:atlas&&<AnatomyScene atlas={atlas} vrModelUrl={cataloguePath(model).replace('.json','.vr.json')} hierarchy={selectedGuest?'guest':hierarchy==='regions'||hierarchy==='depth'?hierarchy:'systems'} guestHierarchy={selectedGuest} state={sceneState} remotePoseRef={connection.following?connection.poseRef:undefined} remoteModel={model} onPose={pose=>{const value={...pose,model:loadedModel.current??model};scenePoseRef.current=value;connection.publishPose(value);}} cameraTransition={slideTransition} sectionTool={show('sections')&&panel==='sections'} onSectionPosition={position=>setState(s=>{const sections=sectionStack(s),active=Math.min(s.activeSection??0,sections.length-1),previous=sections[active];sections[active]={...previous,enabled:true,position};return {...s,sections,activeSection:active,section:sections[active],...(previous.enabled?{}:{camera:undefined})};})} onSelect={(id,toggle)=>{if(!connection.following)choosePart(id,toggle);}} onHidePart={(id,peel)=>{if(!connection.following)hidePart(id,peel);}} onCamera={camera=>{cameraRef.current=camera;}} onCovering={ids=>{if(!connection.following)setCovering(current=>current.some(id=>(state.hidden??[]).includes(id))?current:ids.filter(id=>!peelBaseline.current.includes(id)));}} onExplosionSteps={steps=>setExplodeSteps(current=>current===steps?current:steps)} onCapture={capture=>{captureRef.current=capture;}} onFrameCapture={capture=>{frameCaptureRef.current=capture;}} onProgress={n=>{setProgress(n);if(n===100)setError('');}} onError={setError}/>}
  <div className="vignette"/>
  {(show('model')||show('explode'))&&<header className="identity"><div className="anatomy-choice">{show('model')&&<Select value={model} onValueChange={(value:string|null)=>{if(value==='about'){setDetails(false);setPanel(null);setAbout(true);}else if(modelChoices.some(choice=>choice.value===value)){setModel(value as Model);setHierarchy('systems');}}} items={[...modelChoices,{value:'about',label:'About model coverage'}]}><SelectTrigger aria-label="Choose anatomy model"><SelectValue>{modelChoices.find(choice=>choice.value===model)?.label} · {atlas?atlas.parts.length.toLocaleString():'…'}</SelectValue></SelectTrigger><SelectContent className="anatomy-choice-menu">{modelChoices.map(choice=><SelectItem key={choice.value} value={choice.value}>{choice.label}</SelectItem>)}<SelectItem value="about">About model coverage</SelectItem></SelectContent></Select>}{show('systems')&&<button type="button" className="model-action" aria-label="Layers" title="Layers" aria-controls="anatomy-browser" aria-expanded={embedded||window.innerWidth<768||window.innerHeight<=600?mobileLayersOpen:layersVisible} onClick={toggleLayers}><Layers3 size={17}/></button>}<button type="button" className="model-action" aria-label="Assemble and reset" title="Reset" onClick={reset}><RotateCcw size={17}/></button></div></header>}
  {(!embedded||['search','study','sections','systems','open','download'].some(id=>embedUi?.has(id as EmbedUi)))&&<nav className="top-actions" aria-label="Explorer panels">{show('search')&&<Button variant="ghost" className={`icon-button ${panel==='search'?'active':''}`} onClick={()=>openPanel('search')} aria-label="Find a structure" title="Find a structure"><Search size={18}/></Button>}{show('sections')&&<Button variant="ghost" className={`icon-button ${panel==='sections'?'active':''}`} onClick={()=>openPanel('sections')} aria-label="Sections" title="Sections"><ScanLine size={18}/></Button>}{show('study')&&<Button variant="ghost" className={`icon-button ${panel==='advanced'?'active':''}`} aria-label="Advanced tools" title="Advanced tools" onClick={()=>openPanel('advanced')}><SlidersHorizontal size={18}/></Button>}{embedded&&show('systems')&&!show('model')&&<Button variant="ghost" className="icon-button" aria-label="Layers" title="Layers" onClick={toggleLayers}><Layers3 size={18}/></Button>}{embedded&&show('open')&&<a className="open-full" href={currentViewUrl()} onClick={e=>{e.currentTarget.href=currentViewUrl();}} target="_blank" rel="noopener noreferrer" aria-label="Open full Human Atlas viewer" title="Open full viewer"><ArrowUpRight size={18}/></a>}{!embedded&&<Button variant="ghost" className="icon-button" aria-label="Save or open views" title="Save or open views" disabled={!atlas||progress<100} onClick={share}><BookmarkPlus size={18}/></Button>}{show('download')&&<ExportTools studioRef={studioRef} captureRef={captureRef} frameCaptureRef={frameCaptureRef} currentViewUrl={currentViewUrl} model={model} ready={!!atlas&&progress===100}/>}</nav>}
  {exportStatus&&<p className="export-status glass" role="status">{exportStatus}</p>}
  {panel==='advanced'&&show('study')&&<AdvancedTools connection={connection} onPresentationChange={value=>{advancedPresentationRef.current=value;setAdvancedToolType(value.document.type);setSlidesRunning(value.document.type==='slides'&&value.slideIndex!==null);connection.publishSnapshot();}} remotePresentation={connection.following?connection.remoteSnapshot?.advanced:null} close={()=>setPanel(null)} currentViewUrl={currentViewUrl} captureScene={()=>({model,hierarchy,guestSources:[...guest.sources],state:{...state,camera:cameraRef.current??state.camera}})} openView={openSavedView} model={model} viewName={state.selected.length?structureName(currentChoice.name):modelChoices.find(choice=>choice.value===model)?.label??'Anatomy view'} ready={!!atlas}/>}
  {panel==='sections'&&atlas&&show('sections')&&<SectionTools referenceBounds={dermatomeReference?dermatomeBounds:undefined} showPreview={!dermatomeReference} atlas={atlas} state={state} setState={setState} close={()=>setPanel(null)}/>}
  {shareLink&&<ShareView viewUrl={shareLink} defaultName={state.selected.length?structureName(currentChoice.name):model==='cell'?'Human cell view':`${sex==='female'?'Female':'Male'} anatomy view`} openView={url=>openSavedView(url,true)} close={()=>setShareLink('')}/>}
  {(show('systems')||show('explode'))&&<section id="anatomy-browser" className={`layers-panel glass ${!show('systems')?'explode-only':''} ${mobileLayersOpen?'mobile-open':''} ${layersVisible?'':'is-hidden'}`} aria-label={model==='cell'?'Cell browser':'Anatomy browser'}>
   {show('systems')&&<><div className="hierarchy-header"><div className="hierarchy-tabs" role="tablist" aria-label="Browse anatomy by"><button type="button" role="tab" aria-selected={hierarchy==='systems'} className={hierarchy==='systems'?'active':''} onClick={()=>setHierarchy('systems')}>{model==='cell'?'Functions':'Systems'}</button><button type="button" role="tab" aria-selected={hierarchy==='regions'} className={hierarchy==='regions'?'active':''} onClick={()=>setHierarchy('regions')}>{model==='cell'?'Compartments':'Regions'}</button>{model!=='cell'&&<button type="button" role="tab" aria-selected={hierarchy==='depth'} className={hierarchy==='depth'?'active':''} onClick={()=>setHierarchy('depth')}>Depth</button>}<GuestHierarchyMenu hierarchies={guest.menuHierarchies} active={hierarchy} onSelect={id=>{setHierarchy(`guest:${id}`);void guest.ensure(id).catch(()=>{});setState(s=>({...s,guestQuery:undefined,guestView:undefined,...(id==='dermatomes-myotomes'?{depth:0,depthHidden:[],skinOpacity:1}:{})}));}} onLoad={async url=>{const loaded=await guest.load(url);setHierarchy(`guest:${loaded.id}`);}}/></div><Button variant="ghost" className="mobile-only icon-button" onClick={()=>setMobileLayersOpen(false)} aria-label="Close anatomy browser"><X size={18}/></Button></div>
   {atlas&&(selectedGuest?<GuestTree key={`${model}:${selectedGuest.id}`} atlas={atlas} hierarchy={selectedGuest} state={state} setState={setState} dermatomeReference={dermatomeReference} onChoose={choose} onLoadExtension={async id=>{await guest.extend(selectedGuest.id,id);setState(s=>({...s,guestExtensions:[...new Set([...(s.guestExtensions??[]),id])]}));}} onNavigate={link=>{setHierarchy(`guest:${link.hierarchy}`);void guest.ensure(link.hierarchy).catch(()=>{});setState(s=>({...s,guestQuery:link.id,guestView:undefined}));}}/>:hierarchy.startsWith('guest:')?<div className="guest-tree-loading" role="status">{guest.error||'Loading hierarchy…'}</div>:hierarchy==='depth'&&model!=='cell'?<DepthTree key={model} atlas={atlas} state={state} setState={setState} onChoose={choose}/>:<SystemTree key={model} atlas={atlas} mode={hierarchy==='regions'?'regions':'systems'} state={state} setState={setState} onChoose={choose}/>)}
   {model!=='cell'?<div className="skin-control depth-control explode-control"><div className="explode-label"><label id="depth-control-label">Depth</label><output>{Math.round(depthPosition(state)*100)}%</output></div><div className="explode-slider-row"><Slider aria-labelledby="depth-control-label" min={0} max={100} step={.1} largeStep={2} value={[explosionToSlider(depthPosition(state),DEPTH_LAYERS.length)*100]} getAriaValueText={(_formatted,value)=>{const progress=sliderToExplosion(value/100,DEPTH_LAYERS.length)*DEPTH_LAYERS.length,index=Math.min(DEPTH_LAYERS.length-1,Math.floor(progress));return progress>=DEPTH_LAYERS.length?'All depth layers hidden':DEPTH_LAYERS[index].name+': '+Math.round((1-(progress-index))*100)+'% visible';}} onValueChange={value=>setDepth(Array.isArray(value)?value[0]:value)} onValueCommitted={(value,details)=>{if(details.reason==='drag'||details.reason==='track-press')setDepth(Array.isArray(value)?value[0]:value,true);}}/><div className="explode-ticks" aria-hidden="true">{explosionMarks(DEPTH_LAYERS.length).map((position,index)=><span key={index} style={{left:position*100+'%'}}/>)}</div></div></div>:counts['cell-boundary']>0&&<div className="skin-control"><label id="skin-label">Cell membrane <output>{state.visible.includes('cell-boundary')?Math.round((state.skinOpacity??.1)*100):0}%</output></label><Slider aria-labelledby="skin-label" min={0} max={100} step={5} value={[state.visible.includes('cell-boundary')?(state.skinOpacity??.1)*100:0]} onValueChange={value=>{const opacity=(Array.isArray(value)?value[0]:value)/100;setState(s=>({...s,skinOpacity:opacity,visible:opacity>0?[...s.visible.filter(id=>id!=='cell-boundary'),'cell-boundary']:s.visible.filter(id=>id!=='cell-boundary')}));}}/></div>}
   </>}
  {show('explode')&&<div className="layers-explode"><div className="explode-control"><div className="explode-label"><label id="explode-label">Explode</label><output>{Math.round(state.explode*100)}%</output></div><div className="explode-slider-row"><Slider disabled={dermatomeReference} aria-labelledby="explode-label" min={0} max={100} step={.1} largeStep={2} value={[explosionToSlider(state.explode,explodeSteps)*100]} getAriaValueText={(_formatted,value)=>`${Math.round(sliderToExplosion(value/100,explodeSteps)*100)}% exploded`} onValueChange={value=>setExplosion(Array.isArray(value)?value[0]:value)} onValueCommitted={(value,details)=>{if(details.reason==='drag'||details.reason==='track-press')setExplosion(Array.isArray(value)?value[0]:value,true);}}/><div className="explode-ticks" aria-hidden="true">{explosionMarks(explodeSteps).map((position,index)=><span key={index} style={{left:`${position*100}%`}}/>)}</div></div></div></div>}
  </section>}
  {panel==='search'&&show('search')&&<section className="search-panel glass" aria-label={model==='cell'?'Find cell component':'Find anatomy'}><div className="panel-heading"><span>Find a {model==='cell'?'component':'structure'}</span><Button variant="ghost" className="icon-button" onClick={()=>setPanel(null)} aria-label="Close search"><X size={18}/></Button></div><label className="check-label selection-mode"><input type="checkbox" checked={addSelection} onChange={e=>setAddSelection(e.target.checked)}/>Add to selection set</label><Combobox<Concept> items={results} value={null} onValueChange={value=>{if(value)choose(value);}} inputValue={query} onInputValueChange={setQuery} itemToStringLabel={c=>c.name} filter={null} open onOpenChange={open=>{if(!open)setPanel(null);}}><ComboboxInput autoFocus placeholder={model==='cell'?'Nucleolus, mitochondria, Golgi apparatus…':'Heart, femur, cranial nerve…'} aria-label={model==='cell'?'Search named cell components':'Search named anatomical structures'} showTrigger={false}/><ComboboxContent className="anatomy-search-results"><ComboboxEmpty>No {model==='cell'?'components':'structures'} match your search.</ComboboxEmpty><ComboboxList>{(c:Concept)=>{const display=displayLaterality(c.name);return <ComboboxItem key={c.id} value={c} aria-label={c.name} title={c.name}><span className={`search-result-name ${lateralityClass(display.side)}`}>{display.label}</span><span className="small-number">{c.elements.length} {c.elements.length===1?'piece':'pieces'}</span></ComboboxItem>;}}</ComboboxList></ComboboxContent></Combobox><p className="search-note">{query?'Showing up to 80 matches. Refine your search to find smaller structures.':model==='cell'?'Search a component or organelle.':'Start with a major organ, or search every named structure.'}</p></section>}
  {show('camera')&&<nav className="view-controls glass" aria-label="Camera controls">{(['three-quarter','front','side','back'] as View[]).map((v,i)=><Button variant="ghost" key={v} className={state.view===v?'active':''} aria-pressed={state.view===v} disabled={state.explode>.04&&v!=='front'} onClick={()=>setState(s=>({...s,view:v,reset:s.reset+1,rotate:false,camera:undefined}))} title={`${v} view`} aria-label={`${v} view`}><span>{['¾','F','S','B'][i]}</span></Button>)}<i/><Button variant="ghost" disabled={state.explode>.04} aria-label={state.rotate?'Pause rotation':model==='cell'?'Rotate cell':'Rotate body'} title="Auto rotate" className={state.rotate?'active':''} onClick={()=>setState(s=>({...s,rotate:!s.rotate}))}>{state.rotate?<Pause size={17}/>:<RotateCw size={18}/>}</Button><Button variant="ghost" aria-label="Reset view and layers" title="Reset" onClick={reset}><RotateCcw size={17}/></Button></nav>}
  <div className="scene-caption"><span className="caption-line"/><span>{state.isolate?(currentChoice?.name??'SELECTED STRUCTURE'):state.explode>.95?(model==='cell'?'CELL COMPONENTS':'ANATOMICAL INVENTORY'):state.explode>.05?'SEPARATED STRUCTURES':model==='cell'?'HUMAN EUKARYOTIC CELL':model==='embryo'?(atlas?.version.startsWith('Carnegie stage 23')?'CARNEGIE STAGE 23 · PARTIAL ASSEMBLY':'EMBRYO · ASSOCIATED TISSUES'):sex==='female'?'FEMALE · REFERENCE ANATOMY':'ADULT HUMAN · MALE'}</span><span className="caption-line"/></div>
  {error&&<div className="loading glass error" role="alert"><p>{error}</p><Button variant="ghost" onClick={()=>location.reload()}>Reload viewer</Button></div>}
  {show('details')&&<Sheet open={details&&(selectedParts.length>0||!!chosen)} modal={false} disablePointerDismissal onOpenChange={open=>{if(open)setDetails(true);else clearSelection();}}><SheetContent portalContainer={studioRef} initialFocus={detailTitle} className={`detail-sheet glass ${state.isolate?'is-isolated':''}`} showCloseButton={true}><SelectionInspector titleRef={detailTitle} choice={currentChoice} ancestors={detailAncestors} selectedParts={selectedParts} anchorParts={expansionAnchors} relationshipDepth={activeExpansion?.depth??0} relationshipExhausted={activeExpansion?.exhausted??false} partById={parts} concepts={atlas?.concepts??[]} scope={atlas?.scope} state={state} viewUrl={currentViewUrl()} covering={covering} depth={peelDepth} onDepth={setPeelDepth} onCenter={()=>setState(s=>({...s,focus:(s.focus??0)+1,camera:undefined}))} onIsolate={()=>setState(s=>({...s,isolate:!s.isolate}))} onExpandRelationships={expandRelationships} onShowComponents={showComponents} onHide={hideSelection} onHidePart={hidePart} onChoosePart={choosePart} onChooseChild={choose}/></SheetContent></Sheet>}
  <Sheet open={about} onOpenChange={setAbout}><SheetContent className="about-sheet glass"><div className="eyebrow">MODEL COVERAGE</div><SheetTitle className="structure-title">About this anatomy</SheetTitle><SheetDescription>Explore male and female reference models with different levels of detail.</SheetDescription><div className="about-copy"><p><strong>Male · detailed</strong><br/>Z-Anatomy mesh and curve anatomy, original Open 3D Model upper-limb structures, and selected digestive and pelvic structures from original BodyParts3D 4.0 geometry, including a closed replacement stomach surface. Muscle attachments and fascia are hidden by default. The detailed and full-resolution menu options use the same geometry; the latter initially shows the outer surface.</p>{import.meta.env.DEV&&<p><strong>Female · TA98 review</strong><br/>An independent local female source frame. Missing-term HRA donors are hidden until female source registration and packing review. TA98 material coverage and image textures remain incomplete. This review model has unverified redistribution rights and stays local.</p>}{import.meta.env.DEV&&<p><strong>Male · reference GLBs</strong><br/>Local import of Anatomy Atlas adapted geometry and its embedded textures. Source credits and separate NonCommercial component licenses are in the downloaded NOTICE.md. This model stays outside deployed assets.</p>}<p><strong>Embryo 3month</strong><br/>Carnegie stage 23, specimen 9226: 166 source structures, retaining original source faces and material colors, alongside eight separately licensed Human Reference Atlas placental structures. The source dates to 56–60 days after fertilization; the display name does not establish a three-month specimen. The placental reference is resized and flattened toward the earlier volume target. Adjoining surfaces remain separate objects; shared tissue walls and continuous vessel lumens are unverified. CS23 material retains CC BY-NC-ND 4.0; placental context retains CC BY 4.0. See ATTRIBUTION.md for credits and license terms.</p><p><strong>Male · standard</strong><br/>2,233 individual meshes and 3,429 named concepts built from BodyParts3D 4.0.</p><p><strong>Female anatomy</strong><br/>1,030 meshes and 1,234 named concepts built from original Human Reference Atlas v1.10 and v1.5 files plus 74 original University of Denver lower-limb muscle meshes aligned to HRA bones. The placenta and associated tissues are in the Embryo model. Upper-body muscle and bone coverage remains incomplete.</p><p>These models have different coverage. None contains every human structure or variation. Named concepts can contain multiple pieces; each mesh is rendered once.</p><p>Colors and system groupings are designed for exploration. Short explanations provide general educational context. This is an anatomical reference, not a diagnostic or surgical tool.</p></div></SheetContent></Sheet>
 </main>;
}
