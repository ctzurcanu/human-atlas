import {assetUrl} from './asset-url';
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
import {readEmbedUi,type EmbedUi} from './embed';
import ShareView from './share-view';
import SelectionInspector from './selection-inspector';
import {anatomicalRelations} from './anatomical-relations';
import {expandRelationshipSelection,type RelationshipExpansion} from './relationship-selection';
import {displayLaterality,lateralityClass} from './laterality';
import SystemTree,{type Hierarchy} from './system-tree';
import DepthTree from './depth-tree';
import GuestTree from './guest-tree';
import GuestHierarchyMenu from './guest-hierarchy-menu';
import {useGuestHierarchies} from './use-guest-hierarchies';
import {DEPTH_LAYERS} from './depth-layers';
import {depthLayerFor} from './depth-layers';
import type {HierarchyChoice} from './hierarchy-choice';
import {hierarchyAncestors,hierarchyNavigation} from './hierarchy-navigation';
import {terminologyForConcept} from './anatomical-terminology';
import {flushSync} from 'react-dom';
import {registerAtlasTools,type AgentActions,type AgentViewChange} from './agent-tools';
import {useEffect,useMemo,useRef,useState} from 'react';
import {ArrowUpRight,BookmarkPlus,Download,Layers3,Pause,RotateCcw,RotateCw,ScanLine,Search,SlidersHorizontal,X} from 'lucide-react';
import {Select,SelectTrigger,SelectValue,SelectContent,SelectItem} from '@/components/ui/select';
import {Button} from '@/components/ui/button';
import {Slider} from '@/components/ui/slider';
import {Sheet,SheetContent,SheetTitle,SheetDescription} from '@/components/ui/sheet';
import {Combobox,ComboboxInput,ComboboxContent,ComboboxList,ComboboxItem,ComboboxEmpty} from '@/components/ui/combobox';
import AnatomyScene from './scene';
import {CELL_VISIBLE,DEFAULT_VISIBLE,SYSTEMS,structureName,type Atlas,type Concept,type SceneState,type View} from './anatomy';
const initial:SceneState={focus:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:true,explode:0,visible:DEFAULT_VISIBLE,selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0};
type Model='male-full'|'male-detail'|'male'|'female'|'embryo'|'cell'|'local-male'|'local-female'|'local-reference';
type LayerHierarchy=Hierarchy|'depth'|`guest:${string}`;
const layerHierarchy=(tree:string|null,model:Model):LayerHierarchy=>tree==='regions'||tree==='depth'&&model!=='cell'?tree:tree?.startsWith('guest:')?tree as `guest:${string}`:'systems';
const modelChoices:{value:Model;label:string}[]=[
 {value:'male-full',label:'Male · full resolution + skin'},
 {value:'male-detail',label:'Male · detailed'},
 {value:'male',label:'Male · standard'},
 {value:'female',label:'Female · partial'},
 {value:'embryo',label:'Embryo · partial'},
 {value:'cell',label:'Human eukaryotic cell'},
 ...(import.meta.env.DEV?[{value:'local-reference' as const,label:'Male · reference GLBs'},{value:'local-male' as const,label:'Male · local Zygote'},{value:'local-female' as const,label:'Female · local Zygote'}]:[]),
];
const parseModel=(value:string|null):Model=>modelChoices.find(choice=>choice.value===value)?.value??'male-detail';
const cataloguePath=(model:Model)=>model==='local-reference'?'/local-models/reference.json':model==='local-male'?'/local-models/male.json':model==='local-female'?'/local-models/female.json':model==='male-full'||model==='male-detail'?'/models/atlas-male-complete.json':model==='female'?'/models/atlas-hra-female.json':model==='embryo'?'/models/atlas-embryo.json':model==='cell'?'/models/atlas-cell.json':'/models/atlas.json';
const defaultLayers=(model:Model,sex:'male'|'female')=>model==='embryo'?['pregnancy' as const]:model==='cell'?CELL_VISIBLE:sex==='female'||model==='male-full'||model.startsWith('local-')?[...DEFAULT_VISIBLE,'integumentary' as const]:DEFAULT_VISIBLE;
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
 const sex=model==='female'||model==='embryo'||model==='local-female'?'female':'male';
 const [atlas,setAtlas]=useState<Atlas|null>(null),[progress,setProgress]=useState(0),[error,setError]=useState(''),[panel,setPanel]=useState<'search'|'advanced'|'sections'|null>(()=>embedParams.has('connect')?'advanced':null),[details,setDetails]=useState(false),[about,setAbout]=useState(false),[query,setQuery]=useState(''),[chosen,setChosen]=useState<HierarchyChoice|null>(null);
 const [frontPanel,setFrontPanel]=useState<'layers'|'search'|'advanced'|'sections'|'details'>('layers');
 const cameraRef=useRef<number[]|undefined>(undefined),captureRef=useRef<(()=>Promise<Blob>)|null>(null),pendingUrl=useRef(location.href);
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
 useEffect(()=>{if(loadedModel.current===model&&!pendingUrl.current)return;const abort=new AbortController();setProgress(0);setError('');setAtlas(null);setChosen(null);setRelationshipExpansion(null);setDetails(false);setPanel(current=>current==='advanced'?current:null);setMobileLayersOpen(false);setQuery('');cameraRef.current=undefined;fetch(assetUrl(cataloguePath(model)),{signal:abort.signal}).then(r=>{if(!r.ok)throw new Error(model==='local-reference'?'Reference GLBs are missing. Run npm run import:reference-glb.':model.startsWith('local-')?'Local models are missing. Run npm run import:local-models.':'The anatomy catalogue could not be loaded.');return r.json();}).then(data=>{setAtlas(data as Atlas);const base:SceneState={...initial,skinOpacity:model==='cell'?.18:model==='male-full'||model==='local-reference'?1:model.startsWith('local-')?0:.1,visible:defaultLayers(model,sex)};const restored=pendingUrl.current?readViewUrl(pendingUrl.current,data as Atlas,base):base;const terms=viewParameters(pendingUrl.current).getAll('select');setChosen(terms.length===1?(data as Atlas).concepts.find(concept=>concept.id===terms[0])??null:null);loadState(restored);setDetails(restored.selected.length>0);if(restored.sections?.some(section=>section.enabled)){setPanel(current=>current==='advanced'?current:'sections');setFrontPanel(current=>current==='advanced'?current:'sections');}loadedModel.current=model;pendingUrl.current='';}).catch(e=>{if(e.name!=='AbortError')setError(e.message);});return()=>abort.abort();},[model]);
 useEffect(()=>{const key=(e:KeyboardEvent)=>{if(!followingRef.current&&e.key==='/'&&!(e.target instanceof HTMLInputElement)&&!(e.target instanceof HTMLTextAreaElement)){e.preventDefault();setPanel('search');setFrontPanel('search');setDetails(false);}};window.addEventListener('keydown',key);return()=>window.removeEventListener('keydown',key);},[]);
 useEffect(()=>{const pop=()=>{if(!slideNavigationPending.current)setSlideTransition(undefined);slideNavigationPending.current=false;pendingUrl.current=location.href;const params=viewParameters(location.href),next=parseModel(params.get('model'));setHierarchy(layerHierarchy(params.get('tree'),next));guest.restore(viewParameters(location.href).toString());setRelationshipExpansion(null);if(next!==model)setModel(next);else if(atlas){const terms=params.getAll('select'),restored=readViewUrl(location.href,atlas,initial);setChosen(terms.length===1?atlas.concepts.find(concept=>concept.id===terms[0])??null:null);loadState(restored);setDetails(restored.selected.length>0);if(restored.sections?.some(section=>section.enabled)){setPanel(current=>current==='advanced'?current:'sections');setFrontPanel(current=>current==='advanced'?current:'sections');}pendingUrl.current='';}};window.addEventListener('popstate',pop);return()=>window.removeEventListener('popstate',pop);},[model,atlas,loadState,guest.restore]);
 const currentViewUrl=()=>{const url=new URL(viewUrl(location.href,model,state,cameraRef.current??state.camera));if(hierarchy!=='systems')url.searchParams.set('tree',hierarchy);guest.sources.forEach(source=>url.searchParams.append('guest',source));if(embedParams.has('ui'))url.searchParams.set('ui',embedParams.get('ui')!);return compactViewUrl(url.href);};
 const openSavedView=(url:string,animate:boolean)=>{
  url=compactViewUrl(url);
  const from=cameraRef.current;
  const params=viewParameters(url),target=(params.get('frame')??params.get('camera'))?.split(',').map(Number);
  setSlideTransition(animate&&from&&validCamera(from)&&target&&validCamera(target)?{id:++slideTransitionSerial.current,from:[...from]}:undefined);
  slideNavigationPending.current=true;
  history.pushState(null,'',url);window.dispatchEvent(new PopStateEvent('popstate'));
 };
 const share=()=>{const url=currentViewUrl();history.replaceState(null,'',url);setShareLink(url);};
 const downloadPng=async()=>{try{setExportStatus('Preparing PNG…');const blob=await captureRef.current?.();if(!blob)throw new Error('The 3D view is still loading.');const url=URL.createObjectURL(blob),link=document.createElement('a');link.href=url;link.download=`human-atlas-${model}-${Date.now()}.png`;document.body.appendChild(link);link.click();link.remove();setExportStatus('PNG download started.');setTimeout(()=>{URL.revokeObjectURL(url);setExportStatus('');},60_000);return {ok:true};}catch(error){const message=error instanceof Error?error.message:'Could not export the image.';setExportStatus(message);return {ok:false,error:message};}};
 const parts=useMemo(()=>new Map(atlas?.parts.map(p=>[p.id,p])),[atlas]);
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
 const clearSelection=()=>{setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden),selected:[],isolate:false,camera:cameraRef.current??s.camera}));setChosen(null);setRelationshipExpansion(null);setDetails(false);setCovering([]);};
 const hideSelection=()=>{setState(s=>({...s,hidden:[...new Set([...restoreInspectionHidden(s.hidden),...s.selected])],selected:[],isolate:false,camera:cameraRef.current??s.camera}));setChosen(null);setRelationshipExpansion(null);setDetails(false);setCovering([]);};
 const hidePart=(id:string)=>{
  if(!parts.has(id))return;
  const remaining=state.selected.filter(selectedId=>selectedId!==id);
  setState(s=>({...s,hidden:[...new Set([...restoreInspectionHidden(s.hidden),id])],selected:s.selected.filter(selectedId=>selectedId!==id),isolate:remaining.length?s.isolate:false,camera:remaining.length?s.camera:cameraRef.current??s.camera}));
  if(!state.selected.includes(id))return;
  setRelationshipExpansion(null);
  const next=parts.get(remaining[0]);
  setChosen(next?{id:next.conceptId,name:structureName(next.name),elements:[next.id]}:null);
  if(!next){setDetails(false);setCovering([]);}
 };
 useEffect(()=>{const key=(e:KeyboardEvent)=>{if((e.key!=='Delete'&&e.key!=='Backspace')||e.defaultPrevented||e.metaKey||e.ctrlKey||e.altKey||!state.selected.length)return;const target=e.target;if(target instanceof HTMLElement&&(target.isContentEditable||target.closest('input, textarea, select, [role="textbox"]')))return;e.preventDefault();hideSelection();};window.addEventListener('keydown',key);return()=>window.removeEventListener('keydown',key);},[state.selected.length,covering]);
 const results=useMemo(()=>{if(!atlas)return[];const term=query.toLowerCase().trim(),source=term.replace(/^atlas:/,'');if(!term)return ['heart','brain','liver','stomach','spleen','pancreas','urinary bladder','trachea'].map(name=>atlas.concepts.find(c=>c.name.toLowerCase()===name)).filter((x):x is Concept=>!!x);return atlas.concepts.filter(c=>{const ids=terminologyForConcept(c.id);return c.name.toLowerCase().includes(term)||c.id.toLowerCase().includes(source)||[ids.ta98&&`ta98:${ids.ta98}`,ids.tha,ids.fma,ids.ontology,ids.latin&&`la:${ids.latin}`,ids.latin].some(value=>value?.toLowerCase().includes(term));}).sort((a,b)=>a.name.length-b.name.length).slice(0,80);},[atlas,query]);
 const choose=(c:HierarchyChoice,toggle=false)=>{
  if(followingRef.current)return;
  if(toggle&&c.elements.length&&c.elements.every(id=>state.selected.includes(id))){
   const removed=new Set(c.elements),remaining=state.selected.filter(id=>!removed.has(id)),next=parts.get(remaining[0]);
   setChosen(remaining.length===1&&next?{id:next.conceptId,name:structureName(next.name),elements:[next.id]}:null);
   setRelationshipExpansion(null);
   setState(s=>{const selected=s.selected.filter(id=>!removed.has(id));return {...s,selected,isolate:selected.length?s.isolate:false,camera:selected.length?s.camera:cameraRef.current??s.camera};});
   setDetails(remaining.length>0);
   if(!remaining.length)setCovering([]);
   return;
  }
  setChosen(c);setRelationshipExpansion(null);
  setState(s=>{
   const group=!!c.children,memberIds=new Set(c.elements),members=c.elements.map(id=>parts.get(id)).filter(part=>!!part);
   return {...s,hidden:restoreInspectionHidden(s.hidden).filter(id=>!memberIds.has(id)),selected:toggle||addSelection?[...new Set([...s.selected,...c.elements])]:c.elements,isolate:false,contextOpacity:s.contextOpacity===1?.18:s.contextOpacity,rotate:false,camera:undefined,
    ...(group?{visible:[...new Set([...s.visible,...members.map(part=>part.system)])],depthHidden:s.depthHidden?.filter(id=>!members.some(part=>depthLayerFor(part)===id)),skinOpacity:members.some(part=>depthLayerFor(part)==='skin')&&!(s.skinOpacity??1)?.1:s.skinOpacity}:{}),
   };
  });
  setDetails(true);setFrontPanel('details');setPanel(current=>current==='sections'||current==='advanced'?current:null);
 };
 agentInspectRef.current=choose;
 useEffect(()=>{if(!atlas)return;return registerAtlasTools(atlas,c=>flushSync(()=>agentInspectRef.current(c)),{
  snapshot:()=>agentActionsRef.current?.snapshot(),update:change=>agentActionsRef.current?.update(change),command:action=>agentActionsRef.current?.command(action),url:()=>agentActionsRef.current?.url()??location.href,
 });},[atlas]);
 const choosePart=(id:string,toggle=false)=>{
  const p=parts.get(id);if(!p)return;
  if(toggle&&state.selected.includes(id)){
   const remaining=state.selected.filter(selectedId=>selectedId!==id),next=parts.get(remaining[0]);
   setChosen(remaining.length===1&&next?{id:next.conceptId,name:structureName(next.name),elements:[next.id]}:null);
   setRelationshipExpansion(null);
   setState(s=>{const selected=s.selected.filter(selectedId=>selectedId!==id);return {...s,selected,isolate:selected.length?s.isolate:false,camera:selected.length?s.camera:cameraRef.current??s.camera};});
   setDetails(remaining.length>0);
   if(!remaining.length)setCovering([]);
   return;
  }
  setChosen({id:p.conceptId,name:structureName(p.name),elements:[id]});setRelationshipExpansion(null);
  setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden).filter(hiddenId=>hiddenId!==id),selected:toggle||addSelection?[...new Set([...s.selected,id])]:[id],isolate:false,contextOpacity:s.contextOpacity===1?.18:s.contextOpacity,rotate:false,camera:undefined}));
  setDetails(true);setFrontPanel('details');setPanel(current=>current==='sections'||current==='advanced'?current:null);
 };
 const expandRelationships=()=>{if(!state.selected.length||activeExpansion?.exhausted)return;const result=expandRelationshipSelection(state.selected,activeExpansion,id=>{const part=parts.get(id);return part?anatomicalRelations(part,atlas?.parts??[]).map(relation=>relation.target).filter(target=>!target.suppressed).map(target=>target.id):[];});setRelationshipExpansion(result.expansion);setState(s=>({...s,hidden:restoreInspectionHidden(s.hidden),selected:result.selected,isolate:true,labels:true,explode:0,rotate:false,camera:undefined}));};
 const reset=()=>{setState(s=>({...initial,skinOpacity:model==='cell'?.18:model==='male-full'||model==='local-reference'?1:model.startsWith('local-')?0:.1,visible:defaultLayers(model,sex),reset:s.reset+1}));setChosen(null);setRelationshipExpansion(null);setDetails(false);setPanel(null);};
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
 const setExplosion=(value:number,magnetic=false)=>{const interval=100/Math.max(1,explodeSteps),nearest=Math.round(value/interval)*interval,adjusted=magnetic&&Math.abs(value-nearest)<=Math.min(4,interval*.28)?nearest:value;setState(s=>({...s,explode:adjusted/100,rotate:false,camera:undefined}));};
 const resolveAgentPieces=(terms:string[])=>[...new Set(terms.flatMap(term=>{const part=atlas?.parts.find(p=>p.id===term||p.sourceId===term);if(part)return [part.id];const matches=atlas?.concepts.filter(c=>c.id===term||c.name.toLowerCase()===term.toLowerCase())??[];if(!matches.length)throw new Error(`No structure matches "${term}" in this model.`);return matches.flatMap(c=>c.elements);} ))];
 agentActionsRef.current={
  snapshot:()=>({model,hierarchy,guestHierarchies:guest.hierarchies.map(item=>({id:item.id,name:item.name})),selected:state.selected,systems:state.visible,availableSystems:Object.entries(counts).filter(([,count])=>count>0).map(([id,count])=>({id,pieces:count})),depthHidden:state.depthHidden??[],availableDepthLayers:model==='cell'?[]:DEPTH_LAYERS,hidden:state.hidden??[],view:state.view,region:state.region??'all',context:state.contextOpacity??1,skinOpacity:state.skinOpacity??.1,explode:state.explode,labels:state.labels!==false,isolate:state.isolate,rotate:state.rotate,section:state.section??null,sections:sectionStack(state),activeSection:state.activeSection??0,camera:cameraRef.current??state.camera??null,url:currentViewUrl()}),
  url:currentViewUrl,
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
   if(selected){setChosen(null);setRelationshipExpansion(null);setDetails(selected.length>0);if(selected.length)setFrontPanel('details');setCovering([]);}
   setState(s=>{const sections=change.sections?change.sections.map(item=>item?{enabled:true,...item}:{enabled:false,axis:'axial' as const,position:.38,flip:true}):sectionStack(s),active=Math.min(change.activeSection??s.activeSection??0,sections.length-1);if(change.section!==undefined)sections[active]=change.section?{enabled:true,...change.section}:{enabled:false,axis:'axial',position:.38,flip:true};return {...s,...(selected?{selected}:{}),...(hidden?{hidden}:{}),...(change.systems?{visible:change.systems as SceneState['visible']}:{}),...(change.depthHidden?{depthHidden:change.depthHidden}:{}),...(change.view?{view:change.view as View}:{}),...(change.region?{region:change.region}:{}),...(change.context!==undefined?{contextOpacity:change.context}:selected?.length&&s.contextOpacity===1?{contextOpacity:.18}:{}),...(change.skinOpacity!==undefined?{skinOpacity:change.skinOpacity}:{}),...(change.explode!==undefined?{explode:change.explode,rotate:change.explode>.04?false:s.rotate}:{}),...(change.labels!==undefined?{labels:change.labels}:{}),...(change.isolate!==undefined?{isolate:change.isolate,explode:change.isolate?0:s.explode}:selected?.length===0?{isolate:false}:{}),...(change.rotate!==undefined?{rotate:change.rotate}:{}),...(change.section!==undefined||change.sections!==undefined||change.activeSection!==undefined?{sections,activeSection:active,section:sections[active],camera:change.camera??undefined}:{}),...(change.camera?{camera:change.camera}:{}),...(change.focus?{focus:(s.focus??0)+1}:{}),...((change.view||change.region)?{reset:s.reset+1}:{} )};});
   return {updated:true,...change};
  },
  command:async action=>{if(followingRef.current&&action!=='download-png')throw new Error('This guest view is controlled by its host.');switch(action){case 'undo':undo();break;case 'redo':redo();break;case 'reset':reset();break;case 'focus':setState(s=>({...s,focus:(s.focus??0)+1,camera:undefined}));break;case 'expand-relationships':expandRelationships();break;case 'hide-selection':hideSelection();break;case 'clear-selection':clearSelection();break;case 'download-png':return downloadPng();default:throw new Error('Unknown viewer action.');}return {action,done:true};},
 };
 connectSnapshotRef.current=()=>({version:1,model,hierarchy,guestSources:[...guest.sources],state:{...state,camera:cameraRef.current??state.camera},choice:chosen?{id:chosen.id,name:chosen.name,elements:chosen.elements,terminology:chosen.terminology,group:!!chosen.children,...(chosen.children?{children:chosen.children.map(child=>({id:child.id,name:child.name,elements:child.elements,terminology:child.terminology,group:!!child.children}))}:{})}:null,ui:{panel,frontPanel,details,layersVisible,mobileLayersOpen,addSelection,query,expanded:expandedControls()},covering,relationshipExpansion,pose:scenePoseRef.current?.model===model?scenePoseRef.current:null,advanced:advancedPresentationRef.current?{...advancedPresentationRef.current,document:{...advancedPresentationRef.current.document,views:advancedPresentationRef.current.document.views.map(({scene,...view})=>view)}}:null});
 useEffect(()=>{connection.publishSnapshot();},[model,hierarchy,guest.sources.join('|'),state,chosen,panel,frontPanel,details,layersVisible,mobileLayersOpen,addSelection,query,covering,relationshipExpansion,connection.publishSnapshot]);
 const sceneState=useMemo(()=>({...state,selectionGroup:!!currentChoice?.children&&!activeExpansion,inspectorOpen:show('details')&&details&&selectedParts.length>0}),[state,!!currentChoice?.children,!!activeExpansion,show('details'),details,selectedParts.length]);
 const blockGuestAction=(event:React.SyntheticEvent)=>{if(connection.following&&!replayingUiRef.current&&event.target instanceof Element&&!event.target.closest('.connect-tools,.connection-session,[aria-label="Advanced tools"]')){event.preventDefault();event.stopPropagation();}else if(!connection.following)connection.publishSnapshot();};
 useEffect(()=>{
  const remote=connection.remoteSnapshot;if(!connection.following||!remote||remote.model!==model||!atlas)return;
  let frame=0,attempt=0;const sync=()=>{frame=0;replayingUiRef.current=true;let pending=false;try{pending=applyExpandedControls(remote.ui.expanded??[]);}finally{replayingUiRef.current=false;}if(pending&&++attempt<12)frame=requestAnimationFrame(sync);};
  // Relationships and guest hierarchies can finish loading after the snapshot.
  const watchPanels=()=>studioRef.current?.querySelectorAll('#anatomy-browser,.detail-sheet').forEach(panel=>observer.observe(panel,{childList:true,subtree:true}));
  const observer=new MutationObserver(()=>{watchPanels();if(!frame){attempt=0;frame=requestAnimationFrame(sync);}});if(studioRef.current)observer.observe(studioRef.current,{childList:true});watchPanels();
  frame=requestAnimationFrame(sync);return()=>{cancelAnimationFrame(frame);observer.disconnect();};
 },[connection.remoteSnapshot,connection.following,model,atlas,hierarchy,selectionKey,guest.hierarchies.map(item=>item.id).join('|')]);
 return <main ref={studioRef} className={`studio ${embedded?'embedded':''}`} data-front-panel={frontPanel} data-connect-following={connection.following} data-presenting={slidesRunning||connection.active} onPointerDownCapture={blockGuestAction} onClickCapture={blockGuestAction} onDoubleClickCapture={blockGuestAction} onKeyDownCapture={blockGuestAction}>
  {connection.active&&(panel!=='advanced'||advancedToolType!=='connect')&&<div className="connection-session"><ConnectSessionButton role={connection.role} onOpen={()=>{connection.openControls();setPanel('advanced');setFrontPanel('advanced');}} onStop={connection.stop}/></div>}
  {atlas&&<AnatomyScene atlas={atlas} vrModelUrl={cataloguePath(model).replace('.json','.vr.json')} hierarchy={selectedGuest?'guest':hierarchy==='regions'||hierarchy==='depth'?hierarchy:'systems'} guestHierarchy={selectedGuest} state={sceneState} remotePoseRef={connection.following?connection.poseRef:undefined} remoteModel={model} onPose={pose=>{const value={...pose,model:loadedModel.current??model};scenePoseRef.current=value;connection.publishPose(value);}} cameraTransition={slideTransition} sectionTool={show('sections')&&panel==='sections'} onSectionPosition={position=>setState(s=>{const sections=sectionStack(s),active=Math.min(s.activeSection??0,sections.length-1),previous=sections[active];sections[active]={...previous,enabled:true,position};return {...s,sections,activeSection:active,section:sections[active],...(previous.enabled?{}:{camera:undefined})};})} onSelect={(id,toggle)=>{if(!connection.following)choosePart(id,toggle);}} onHidePart={id=>{if(!connection.following)hidePart(id);}} onCamera={camera=>{cameraRef.current=camera;}} onCovering={ids=>{if(!connection.following)setCovering(current=>current.some(id=>(state.hidden??[]).includes(id))?current:ids.filter(id=>!peelBaseline.current.includes(id)));}} onExplosionSteps={steps=>setExplodeSteps(current=>current===steps?current:steps)} onCapture={capture=>{captureRef.current=capture;}} onProgress={n=>{setProgress(n);if(n===100)setError('');}} onError={setError}/>}
  <div className="vignette"/>
  {(show('model')||show('explode'))&&<header className="identity"><div className="anatomy-choice">{show('model')&&<Select value={model} onValueChange={(value:string|null)=>{if(value==='about'){setDetails(false);setPanel(null);setAbout(true);}else if(modelChoices.some(choice=>choice.value===value)){setModel(value as Model);setHierarchy('systems');}}} items={[...modelChoices,{value:'about',label:'About model coverage'}]}><SelectTrigger aria-label="Choose anatomy model"><SelectValue>{modelChoices.find(choice=>choice.value===model)?.label} · {atlas?atlas.parts.length.toLocaleString():'…'}</SelectValue></SelectTrigger><SelectContent className="anatomy-choice-menu">{modelChoices.map(choice=><SelectItem key={choice.value} value={choice.value}>{choice.label}</SelectItem>)}<SelectItem value="about">About model coverage</SelectItem></SelectContent></Select>}{show('systems')&&<button type="button" className="model-action" aria-label="Layers" title="Layers" aria-controls="anatomy-browser" aria-expanded={embedded||window.innerWidth<768||window.innerHeight<=600?mobileLayersOpen:layersVisible} onClick={toggleLayers}><Layers3 size={17}/></button>}<button type="button" className="model-action" aria-label="Assemble and reset" title="Reset" onClick={reset}><RotateCcw size={17}/></button></div></header>}
  {(!embedded||['search','study','sections','systems','open','download'].some(id=>embedUi?.has(id as EmbedUi)))&&<nav className="top-actions" aria-label="Explorer panels">{show('search')&&<Button variant="ghost" className={`icon-button ${panel==='search'?'active':''}`} onClick={()=>openPanel('search')} aria-label="Find a structure" title="Find a structure"><Search size={18}/></Button>}{show('sections')&&<Button variant="ghost" className={`icon-button ${panel==='sections'?'active':''}`} onClick={()=>openPanel('sections')} aria-label="Sections" title="Sections"><ScanLine size={18}/></Button>}{show('study')&&<Button variant="ghost" className={`icon-button ${panel==='advanced'?'active':''}`} aria-label="Advanced tools" title="Advanced tools" onClick={()=>openPanel('advanced')}><SlidersHorizontal size={18}/></Button>}{embedded&&show('systems')&&!show('model')&&<Button variant="ghost" className="icon-button" aria-label="Layers" title="Layers" onClick={toggleLayers}><Layers3 size={18}/></Button>}{embedded&&show('open')&&<a className="open-full" href={currentViewUrl()} onClick={e=>{e.currentTarget.href=currentViewUrl();}} target="_blank" rel="noopener noreferrer" aria-label="Open full Human Atlas viewer" title="Open full viewer"><ArrowUpRight size={18}/></a>}{!embedded&&<Button variant="ghost" className="icon-button" aria-label="Save or open views" title="Save or open views" disabled={!atlas||progress<100} onClick={share}><BookmarkPlus size={18}/></Button>}{show('download')&&<Button variant="ghost" className="icon-button" aria-label="Download PNG of 3D view" title="Download PNG" onClick={()=>void downloadPng()}><Download size={18}/></Button>}</nav>}
  {exportStatus&&<p className="export-status glass" role="status">{exportStatus}</p>}
  {panel==='advanced'&&show('study')&&<AdvancedTools connection={connection} onPresentationChange={value=>{advancedPresentationRef.current=value;setAdvancedToolType(value.document.type);setSlidesRunning(value.document.type==='slides'&&value.slideIndex!==null);connection.publishSnapshot();}} remotePresentation={connection.following?connection.remoteSnapshot?.advanced:null} close={()=>setPanel(null)} currentViewUrl={currentViewUrl} captureScene={()=>({model,hierarchy,guestSources:[...guest.sources],state:{...state,camera:cameraRef.current??state.camera}})} openView={openSavedView} model={model} viewName={state.selected.length?structureName(currentChoice.name):modelChoices.find(choice=>choice.value===model)?.label??'Anatomy view'} ready={!!atlas}/>}
  {panel==='sections'&&atlas&&show('sections')&&<SectionTools atlas={atlas} state={state} setState={setState} close={()=>setPanel(null)}/>}
  {shareLink&&<ShareView viewUrl={shareLink} defaultName={state.selected.length?structureName(currentChoice.name):model==='cell'?'Human cell view':`${sex==='female'?'Female':'Male'} anatomy view`} close={()=>setShareLink('')}/>}
  {show('systems')&&<section id="anatomy-browser" className={`layers-panel glass ${mobileLayersOpen?'mobile-open':''} ${layersVisible?'':'is-hidden'}`} aria-label={model==='cell'?'Cell browser':'Anatomy browser'}>
   <div className="hierarchy-header"><div className="hierarchy-tabs" role="tablist" aria-label="Browse anatomy by"><button type="button" role="tab" aria-selected={hierarchy==='systems'} className={hierarchy==='systems'?'active':''} onClick={()=>setHierarchy('systems')}>{model==='cell'?'Functions':'Systems'}</button><button type="button" role="tab" aria-selected={hierarchy==='regions'} className={hierarchy==='regions'?'active':''} onClick={()=>setHierarchy('regions')}>{model==='cell'?'Compartments':'Regions'}</button>{model!=='cell'&&<button type="button" role="tab" aria-selected={hierarchy==='depth'} className={hierarchy==='depth'?'active':''} onClick={()=>setHierarchy('depth')}>Depth</button>}<GuestHierarchyMenu hierarchies={guest.menuHierarchies} active={hierarchy} onSelect={id=>{setHierarchy(`guest:${id}`);void guest.ensure(id).catch(()=>{});setState(s=>({...s,guestQuery:undefined,guestView:undefined}));}} onLoad={async url=>{const loaded=await guest.load(url);setHierarchy(`guest:${loaded.id}`);}}/></div><Button variant="ghost" className="mobile-only icon-button" onClick={()=>setMobileLayersOpen(false)} aria-label="Close anatomy browser"><X size={18}/></Button></div>
   {atlas&&(selectedGuest?<GuestTree key={`${model}:${selectedGuest.id}`} atlas={atlas} hierarchy={selectedGuest} state={state} setState={setState} onChoose={choose} onLoadExtension={async id=>{await guest.extend(selectedGuest.id,id);setState(s=>({...s,guestExtensions:[...new Set([...(s.guestExtensions??[]),id])]}));}} onNavigate={link=>{setHierarchy(`guest:${link.hierarchy}`);void guest.ensure(link.hierarchy).catch(()=>{});setState(s=>({...s,guestQuery:link.id,guestView:undefined}));}}/>:hierarchy.startsWith('guest:')?<div className="guest-tree-loading" role="status">{guest.error||'Loading hierarchy…'}</div>:hierarchy==='depth'&&model!=='cell'?<DepthTree key={model} atlas={atlas} state={state} setState={setState} onChoose={choose}/>:<SystemTree key={model} atlas={atlas} mode={hierarchy==='regions'?'regions':'systems'} state={state} setState={setState} onChoose={choose}/>)}
   {hierarchy!=='depth'&&(counts.integumentary>0||counts['cell-boundary']>0)&&<div className="skin-control"><label id="skin-label">{model==='cell'?'Cell membrane':'Skin / body surface'} <output>{state.visible.includes(model==='cell'?'cell-boundary':'integumentary')?Math.round((state.skinOpacity??.1)*100):0}%</output></label><Slider aria-labelledby="skin-label" min={0} max={100} step={5} value={[state.visible.includes(model==='cell'?'cell-boundary':'integumentary')?(state.skinOpacity??.1)*100:0]} onValueChange={v=>{const opacity=(Array.isArray(v)?v[0]:v)/100,surface=model==='cell'?'cell-boundary':'integumentary';setState(s=>({...s,skinOpacity:opacity,visible:opacity>0?[...s.visible.filter(x=>x!==surface),surface]:s.visible.filter(x=>x!==surface),depthHidden:opacity>0?s.depthHidden?.filter(x=>x!=='skin'):s.depthHidden}));}}/>{model!=='cell'&&<p>Outer surface reference; not separate skin layers.</p>}</div>}
  </section>}
  {panel==='search'&&show('search')&&<section className="search-panel glass" aria-label={model==='cell'?'Find cell component':'Find anatomy'}><div className="panel-heading"><span>Find a {model==='cell'?'component':'structure'}</span><Button variant="ghost" className="icon-button" onClick={()=>setPanel(null)} aria-label="Close search"><X size={18}/></Button></div><label className="check-label selection-mode"><input type="checkbox" checked={addSelection} onChange={e=>setAddSelection(e.target.checked)}/>Add to selection set</label><Combobox<Concept> items={results} value={null} onValueChange={value=>{if(value)choose(value);}} inputValue={query} onInputValueChange={setQuery} itemToStringLabel={c=>c.name} filter={null} open onOpenChange={open=>{if(!open)setPanel(null);}}><ComboboxInput autoFocus placeholder={model==='cell'?'Nucleolus, mitochondria, Golgi apparatus…':'Heart, femur, cranial nerve…'} aria-label={model==='cell'?'Search named cell components':'Search named anatomical structures'} showTrigger={false}/><ComboboxContent className="anatomy-search-results"><ComboboxEmpty>No {model==='cell'?'components':'structures'} match your search.</ComboboxEmpty><ComboboxList>{(c:Concept)=>{const display=displayLaterality(c.name);return <ComboboxItem key={c.id} value={c} aria-label={c.name} title={c.name}><span className={`search-result-name ${lateralityClass(display.side)}`}>{display.label}</span><span className="small-number">{c.elements.length} {c.elements.length===1?'piece':'pieces'}</span></ComboboxItem>;}}</ComboboxList></ComboboxContent></Combobox><p className="search-note">{query?'Showing up to 80 matches. Refine your search to find smaller structures.':model==='cell'?'Search a component or organelle.':'Start with a major organ, or search every named structure.'}</p></section>}
  {show('camera')&&<nav className="view-controls glass" aria-label="Camera controls">{(['three-quarter','front','side','back'] as View[]).map((v,i)=><Button variant="ghost" key={v} className={state.view===v?'active':''} aria-pressed={state.view===v} disabled={state.explode>.04&&v!=='front'} onClick={()=>setState(s=>({...s,view:v,reset:s.reset+1,rotate:false,camera:undefined}))} title={`${v} view`} aria-label={`${v} view`}><span>{['¾','F','S','B'][i]}</span></Button>)}<i/><Button variant="ghost" disabled={state.explode>.04} aria-label={state.rotate?'Pause rotation':model==='cell'?'Rotate cell':'Rotate body'} title="Auto rotate" className={state.rotate?'active':''} onClick={()=>setState(s=>({...s,rotate:!s.rotate}))}>{state.rotate?<Pause size={17}/>:<RotateCw size={18}/>}</Button><Button variant="ghost" aria-label="Reset view and layers" title="Reset" onClick={reset}><RotateCcw size={17}/></Button></nav>}
  <div className="scene-caption"><span className="caption-line"/><span>{state.isolate?(currentChoice?.name??'SELECTED STRUCTURE'):state.explode>.95?(model==='cell'?'CELL COMPONENTS':'ANATOMICAL INVENTORY'):state.explode>.05?'SEPARATED STRUCTURES':model==='cell'?'HUMAN EUKARYOTIC CELL':model==='embryo'?'EMBRYO · ASSOCIATED TISSUES':sex==='female'?'FEMALE · REFERENCE ANATOMY':'ADULT HUMAN · MALE'}</span><span className="caption-line"/></div>
  {show('explode')&&<div className="bottom-dock glass"><div className="explode-control"><div className="explode-label"><label id="explode-label">Explode</label><output>{String(Math.round(state.explode*100)).padStart(2,'0')}%</output></div><div className="explode-slider-row"><Slider aria-labelledby="explode-label" min={0} max={100} step={1} value={[state.explode*100]} onValueChange={(value,details)=>setExplosion(Array.isArray(value)?value[0]:value,details.reason==='drag'||details.reason==='track-press')} onValueCommitted={(value,details)=>{if(details.reason==='drag'||details.reason==='track-press')setExplosion(Array.isArray(value)?value[0]:value,true);}}/><div className="explode-ticks" aria-hidden="true">{Array.from({length:explodeSteps+1},(_,index)=><span key={index} style={{left:`${index/Math.max(1,explodeSteps)*100}%`}}/>)}</div></div></div></div>}
  {error&&<div className="loading glass error" role="alert"><p>{error}</p><Button variant="ghost" onClick={()=>location.reload()}>Reload viewer</Button></div>}
  {show('details')&&<Sheet open={details&&selectedParts.length>0} modal={false} disablePointerDismissal onOpenChange={open=>{if(open)setDetails(true);else clearSelection();}}><SheetContent portalContainer={studioRef} initialFocus={detailTitle} className={`detail-sheet glass ${state.isolate?'is-isolated':''}`} showCloseButton={true}><SelectionInspector titleRef={detailTitle} choice={currentChoice} ancestors={detailAncestors} selectedParts={selectedParts} anchorParts={expansionAnchors} relationshipDepth={activeExpansion?.depth??0} relationshipExhausted={activeExpansion?.exhausted??false} partById={parts} scope={atlas?.scope} state={state} covering={covering} depth={peelDepth} onOpacity={value=>setState(s=>({...s,contextOpacity:value/100}))} onDepth={setPeelDepth} onCenter={()=>setState(s=>({...s,focus:(s.focus??0)+1,camera:undefined}))} onIsolate={()=>setState(s=>({...s,isolate:!s.isolate,explode:0}))} onExpandRelationships={expandRelationships} onHide={hideSelection} onHidePart={hidePart} onChoosePart={choosePart} onChooseChild={choose}/></SheetContent></Sheet>}
  <Sheet open={about} onOpenChange={setAbout}><SheetContent className="about-sheet glass"><div className="eyebrow">MODEL COVERAGE</div><SheetTitle className="structure-title">About this anatomy</SheetTitle><SheetDescription>Explore male and female reference models with different levels of detail.</SheetDescription><div className="about-copy"><p><strong>Male · detailed</strong><br/>Z-Anatomy mesh and curve anatomy, original Open 3D Model upper-limb structures, and selected digestive and pelvic structures from original BodyParts3D 4.0 geometry, including a closed replacement stomach surface. Muscle attachments and fascia are hidden by default. The detailed and full-resolution menu options use the same geometry; the latter initially shows the outer surface.</p>{import.meta.env.DEV&&<p><strong>Male · reference GLBs</strong><br/>Local import of Anatomy Atlas adapted geometry and its embedded textures. Source credits and separate NonCommercial component licenses are in the downloaded NOTICE.md. This model stays outside deployed assets.</p>}<p><strong>Male · standard</strong><br/>2,233 individual meshes and 3,429 named concepts built from BodyParts3D 4.0.</p><p><strong>Female anatomy</strong><br/>1,030 meshes and 1,234 named concepts built from original Human Reference Atlas v1.10 and v1.5 files plus 74 original University of Denver lower-limb muscle meshes aligned to HRA bones. The placenta and associated tissues are in the Embryo model. Upper-body muscle and bone coverage remains incomplete.</p><p>These models have different coverage. None contains every human structure or variation. Named concepts can contain multiple pieces; each mesh is rendered once.</p><p>Colors and system groupings are designed for exploration. Short explanations provide general educational context. This is an anatomical reference, not a diagnostic or surgical tool.</p></div></SheetContent></Sheet>
 </main>;
}
