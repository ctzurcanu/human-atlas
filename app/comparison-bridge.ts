import type {Atlas,Part,SceneState} from './anatomy';
import type {ScenePose} from './connect-state';
import {terminologyForConcept} from './anatomical-terminology';
import {validCamera} from '../shared/camera-frame.mjs';
import {sectionStack} from './section-stack';
import type {SectionState} from './section-plane';
import {viewParameters} from './view-url';
import {validSceneWorldPose} from './scene-world-pose';

export const COMPARISON_CHANNEL='atlas-compare/v1';
export type ComparisonSide='left'|'right';
export type ComparisonSession={session:string;side:ComparisonSide};
export type ComparisonSelection={id:string;conceptId:string;name:string;ta98?:string;fma?:string;ta98Kind?:'exact'|'parent';side?:ComparisonSide};
export type ComparisonSnapshot={selection:ComparisonSelection[];hidden:ComparisonSelection[];systems:string[];hierarchy:string;context:number;skinOpacity:number;labels:boolean;isolate:boolean;depth?:number;depthHidden:string[];explode?:number;sections?:SectionState[];activeSection?:number};

/** A standalone viewer never accepts comparison commands. */
export function comparisonSession(query:URLSearchParams,framed:boolean):ComparisonSession|null{
 const session=query.get('compareSession'),side=query.get('compareSide');
 return framed&&session&&/^[A-Za-z0-9_-]{8,128}$/.test(session)&&(side==='left'||side==='right')?{session,side}:null;
}

export function comparisonControlVisible(side:ComparisonSide|null,control:string):boolean{
 if(!side)return true;
 if(control==='model')return true;
 if(control==='camera')return false;
 return control==='systems'||control==='explode'?side==='left':side==='right';
}

/** A model choice loads that pane's source, without the old catalogue's indices. */
export function comparisonModelUrl(current:string,model:string,comparison:ComparisonSession):string{
 if(!/^[a-z0-9-]{1,80}$/.test(model))throw Error('Choose a valid anatomy model.');
 const parameters=viewParameters(current),url=new URL(current);url.search='';url.hash='';
 const presentation=['ui','embed','tree','guest','layers','l','skin','context','dp','depth','labels','explode','view','region','cut','slice','flip','azimuth','elevation','sections','cut2','slice2','flip2','azimuth2','elevation2','sectionTab'];
 for(const key of presentation)for(const value of parameters.getAll(key))url.searchParams.append(key,value);
 url.searchParams.set('model',model);url.searchParams.set('compareSession',comparison.session);url.searchParams.set('compareSide',comparison.side);
 return url.href;
}

function laterality(part:Pick<Part,'id'|'name'|'conceptId'>):ComparisonSide|undefined{
 const value=`${part.name} ${part.id} ${part.conceptId}`;
 const left=/(?:\bleft\b|\.l(?:\b|$))/i.test(value),right=/(?:\bright\b|\.r(?:\b|$))/i.test(value);
 return left!==right?left?'left':'right':undefined;
}

/** Exact terminology is portable; a parent-term placement is not an identity. */
const descriptorCache=new WeakMap<Atlas,{byId:Map<string,ComparisonSelection>;byConcept:Map<string,ComparisonSelection[]>;byTa:Map<string,ComparisonSelection[]>;byFma:Map<string,ComparisonSelection[]>}>();
function descriptorsFor(atlas:Atlas){
 const cached=descriptorCache.get(atlas);if(cached)return cached;
 const concepts=new Map(atlas.concepts.map(concept=>[concept.id,concept]));
 const hidden=new Set(atlas.defaultHidden??[]),active=new Set(atlas.parts.filter(part=>!part.suppressed&&!part.referenceOnly&&!hidden.has(part.id)).map(part=>part.id));
 const descriptors=atlas.parts.map(part=>{
  const concept=concepts.get(part.conceptId),term=terminologyForConcept(part.conceptId,part.ta98Term??concept?.ta98Term,concept?.ta98Kind),side=laterality(part);
  return {id:part.id,conceptId:part.conceptId,name:part.name,...(term.ta98?{ta98:term.ta98}:{}),...(term.fma?{fma:term.fma}:{}),...(concept?.ta98Kind?{ta98Kind:concept.ta98Kind}:{}),...(side?{side}:{})};
 });
 const index=(key:'conceptId'|'ta98'|'fma')=>{const map=new Map<string,ComparisonSelection[]>();for(const item of descriptors){const value=item[key];if(value&&active.has(item.id)&&(key!=='conceptId'||concepts.get(item.conceptId)?.ta98Kind!=='parent'))map.set(value,[...map.get(value)??[],item]);}return map;};
 const result={byId:new Map(descriptors.map(item=>[item.id,item])),byConcept:index('conceptId'),byTa:index('ta98'),byFma:index('fma')};descriptorCache.set(atlas,result);return result;
}
export function describeComparisonSelection(atlas:Atlas,ids:string[]):ComparisonSelection[]{
 const index=descriptorsFor(atlas);return ids.flatMap(id=>{const item=index.byId.get(id);return item?[item]:[];});
}

/** Never guess by a similar name or substitute a different side. */
export function resolveComparisonSelection(atlas:Atlas,selection:ComparisonSelection[]):string[]{
 const {byId,byConcept,byTa,byFma}=descriptorsFor(atlas);
 const output=new Set<string>();
 for(const item of selection){
  const exact=byId.get(item.id);if(exact){output.add(exact.id);continue;}
  if(item.ta98Kind==='parent')continue;
  const sameSide=(other:ComparisonSelection)=>!item.side||other.side===item.side;
  let matches=(byConcept.get(item.conceptId)??[]).filter(sameSide);
  if(!matches.length&&item.ta98)matches=(byTa.get(item.ta98)??[]).filter(sameSide);
  if(!matches.length&&item.fma)matches=(byFma.get(item.fma)??[]).filter(sameSide);
  for(const other of matches)output.add(other.id);
 }
 return [...output];
}

export function comparisonSnapshot(atlas:Atlas,state:SceneState,hierarchy:string):ComparisonSnapshot{
 const defaults=new Set(atlas.defaultHidden??[]);
 return {selection:describeComparisonSelection(atlas,state.selected),hidden:describeComparisonSelection(atlas,(state.hidden??[]).filter(id=>!defaults.has(id))),systems:[...state.visible],hierarchy,context:state.contextOpacity??1,skinOpacity:state.skinOpacity??.1,labels:state.labels!==false,isolate:state.isolate,...(state.depth===undefined?{}:{depth:state.depth}),depthHidden:[...state.depthHidden??[]],explode:state.explode,sections:sectionStack(state).map(section=>({...section})),activeSection:Math.min(state.activeSection??0,sectionStack(state).length-1)};
}

/** Inspection synchronizes without changing either catalogue or camera intent. */
export function applyComparisonSnapshot(atlas:Atlas,state:SceneState,snapshot:ComparisonSnapshot):SceneState{
 const selected=resolveComparisonSelection(atlas,snapshot.selection),selectedSet=new Set(selected),hidden=[...new Set([...(atlas.defaultHidden??[]),...resolveComparisonSelection(atlas,snapshot.hidden)])].filter(id=>!selectedSet.has(id)),available=new Set(atlas.parts.map(part=>part.system));
 const sections=snapshot.sections?.map(section=>({...section})),activeSection=sections?Math.min(snapshot.activeSection??0,sections.length-1):undefined;
 return {...state,selected,hidden,visible:snapshot.systems.filter(system=>available.has(system as Part['system'])) as SceneState['visible'],contextOpacity:snapshot.context,skinOpacity:snapshot.skinOpacity,labels:snapshot.labels,isolate:snapshot.isolate&&selected.length>0,depth:snapshot.depth,depthHidden:[...snapshot.depthHidden],...(snapshot.explode===undefined?{}:{explode:snapshot.explode,...(snapshot.explode>0?{rotate:false}:{})}),...(sections?{sections,activeSection,section:sections[activeSection!]}:{})};
}

const record=(value:unknown):value is Record<string,unknown>=>!!value&&typeof value==='object'&&!Array.isArray(value);
const strings=(value:unknown)=>Array.isArray(value)&&value.length<=30000&&value.every(item=>typeof item==='string'&&item.length<=2048);
const fraction=(value:unknown)=>typeof value==='number'&&Number.isFinite(value)&&value>=0&&value<=1;
function validSelection(value:unknown):value is ComparisonSelection[]{
 return Array.isArray(value)&&value.length<=30000&&value.every(item=>record(item)&&typeof item.id==='string'&&item.id.length<=2048&&typeof item.conceptId==='string'&&item.conceptId.length<=2048&&typeof item.name==='string'&&item.name.length<=2048&&[item.ta98,item.fma].every(code=>code===undefined||typeof code==='string'&&code.length<=100)&&(item.ta98Kind===undefined||item.ta98Kind==='exact'||item.ta98Kind==='parent')&&(item.side===undefined||item.side==='left'||item.side==='right'));
}
export function validComparisonSnapshot(value:unknown):value is ComparisonSnapshot{
 return record(value)&&validSelection(value.selection)&&validSelection(value.hidden)&&strings(value.systems)&&typeof value.hierarchy==='string'&&value.hierarchy.length<=100&&fraction(value.context)&&fraction(value.skinOpacity)&&typeof value.labels==='boolean'&&typeof value.isolate==='boolean'&&(value.depth===undefined||fraction(value.depth))&&strings(value.depthHidden)&&(value.explode===undefined||fraction(value.explode))&&(value.sections===undefined||Array.isArray(value.sections)&&value.sections.length>=1&&value.sections.length<=2&&value.sections.every(validSection))&&(value.activeSection===undefined||Number.isInteger(value.activeSection)&&Number(value.activeSection)>=0&&Number(value.activeSection)<(Array.isArray(value.sections)?value.sections.length:2));
}
function validSection(value:unknown):value is SectionState{
 return record(value)&&typeof value.enabled==='boolean'&&['axial','sagittal','coronal','oblique'].includes(String(value.axis))&&fraction(value.position)&&typeof value.flip==='boolean'&&(value.azimuth===undefined||typeof value.azimuth==='number'&&Number.isFinite(value.azimuth)&&value.azimuth>=-180&&value.azimuth<=180)&&(value.elevation===undefined||typeof value.elevation==='number'&&Number.isFinite(value.elevation)&&value.elevation>=-90&&value.elevation<=90);
}
export function validComparisonPose(value:unknown):value is ScenePose{
 if(!record(value)||typeof value.model!=='string'||!validCamera(value.camera)||!Array.isArray(value.up)||value.up.length!==3||!value.up.every(Number.isFinite)||!Array.isArray(value.rotation)||value.rotation.length!==4||!value.rotation.every(Number.isFinite))return false;
 if(Math.hypot(...value.up)<1e-12||Math.hypot(...value.rotation)<1e-12)return false;
 return value.world===undefined||validSceneWorldPose(value.world);
}

const metreBodyModels=new Set(['male-ta98','male-full','male-detail','male','female','local-male','local-female','local-reference','local-ta98','local-female-ta98']);
export const comparisonSharesWorld=(first:string,second:string)=>first===second||metreBodyModels.has(first)&&metreBodyModels.has(second);

export function comparisonMessage(event:Pick<MessageEvent,'origin'|'source'|'data'>,expected:{origin:string;parent:MessageEventSource|null;session:string;side:ComparisonSide;instance?:string}):Record<string,unknown>|null{
 const data=event.data;
 return !!expected.parent&&event.origin===expected.origin&&event.source===expected.parent&&record(data)&&data.channel===COMPARISON_CHANNEL&&data.session===expected.session&&data.side===expected.side&&(!expected.instance||data.instance===expected.instance||data.type==='request-snapshot'&&data.instance===undefined)?data:null;
}

const documentInstances=new WeakMap<Window,{instance:string;seq:number}>();
function documentInstance(win:Window):{instance:string;seq:number}{
 const cached=documentInstances.get(win);if(cached)return cached;
 const crypto=win.crypto??globalThis.crypto;
 const token=crypto?.randomUUID?.()??`${Date.now().toString(36)}-${Math.random().toString(36).slice(2)}-${Math.random().toString(36).slice(2)}`;
 const state={instance:token,seq:0};documentInstances.set(win,state);return state;
}

type BridgeOptions=ComparisonSession&{window:Window;identity:{model:string;source:string};getSnapshot:()=>ComparisonSnapshot;getViewUrl:()=>string;getPose?:()=>ScenePose|null|undefined;applySnapshot:(snapshot:ComparisonSnapshot)=>ComparisonSnapshot;applyPose:(pose:ScenePose)=>void};
export type ComparisonBridge={publishSnapshot:(snapshot:ComparisonSnapshot)=>void;publishPose:(pose:ScenePose)=>void;dispose:()=>void};

/** Parent messages carry view intent only: they cannot navigate or replace data. */
export function createComparisonBridge(options:BridgeOptions):ComparisonBridge{
 const win=options.window,origin=win.location.origin,scope=documentInstance(win),instance=scope.instance;let disposed=false,stateRevision=-1,poseRevision=-1,lastSnapshot='',expectedSnapshot:string|null=null,poseOwned=true,interacting=false,lastPose='',poseFrame=0,pendingPose:ScenePose|null=null,currentPose:ScenePose|null=null,interactionTimer=0;
 const post=(type:string,payload:Record<string,unknown>={})=>{if(!disposed)win.parent.postMessage({channel:COMPARISON_CHANNEL,session:options.session,side:options.side,instance,seq:++scope.seq,type,...payload},origin);};
 const snapshotKey=(snapshot:ComparisonSnapshot)=>JSON.stringify(snapshot);
 const publishSnapshot=(snapshot:ComparisonSnapshot)=>{
  const key=snapshotKey(snapshot);
  if(expectedSnapshot!==null){if(key===expectedSnapshot){expectedSnapshot=null;lastSnapshot=key;}return;}
  if(key===lastSnapshot)return;lastSnapshot=key;post('snapshot',{snapshot});
 };
 const flushPose=()=>{poseFrame=0;const pose=pendingPose;pendingPose=null;if(!pose||!poseOwned)return;const key=JSON.stringify(pose.world?{world:pose.world,rotation:pose.rotation}:pose);if(key===lastPose)return;lastPose=key;post('pose',{pose});};
 const publishPose=(pose:ScenePose)=>{if(disposed)return;currentPose=pose;if(!poseOwned)return;pendingPose=pose;if(!poseFrame)poseFrame=win.requestAnimationFrame(flushPose);};
 const interactionStart=()=>{interacting=true;poseOwned=true;expectedSnapshot=null;post('interaction',{phase:'start'});};
 const interactionEnd=()=>{interacting=false;post('interaction',{phase:'end'});};
 const wheel=()=>{if(interactionTimer)win.clearTimeout(interactionTimer);interactionStart();interactionTimer=win.setTimeout(interactionEnd,180);};
 const receive=(event:MessageEvent)=>{
  const data=comparisonMessage(event,{origin,parent:win.parent,session:options.session,side:options.side,instance});if(!data)return;
  if(data.type==='request-view'&&typeof data.requestId==='string'&&/^[A-Za-z0-9_-]{1,128}$/.test(data.requestId)){post('view',{requestId:data.requestId,url:options.getViewUrl()});return;}
  if(data.type==='request-snapshot'){post('ready',{identity:options.identity});lastSnapshot='';expectedSnapshot=null;publishSnapshot(options.getSnapshot());const pose=options.getPose?.()??currentPose;if(pose&&validComparisonPose(pose))post('pose',{pose});return;}
  if(data.type!=='apply'||interacting||!Number.isSafeInteger(data.revision)||Number(data.revision)<0)return;
  if(data.kind==='state'&&Number(data.revision)>stateRevision&&validComparisonSnapshot(data.snapshot)){
   stateRevision=Number(data.revision);expectedSnapshot=snapshotKey(options.applySnapshot(data.snapshot));
  }else if(data.kind==='pose'&&Number(data.revision)>poseRevision&&validComparisonPose(data.pose)){
   poseRevision=Number(data.revision);poseOwned=false;pendingPose=null;
   const {world,...normalized}=data.pose;options.applyPose({...(comparisonSharesWorld(data.pose.model,options.identity.model)?data.pose:normalized),model:options.identity.model});
  }
 };
 win.addEventListener('message',receive);win.addEventListener('pointerdown',interactionStart,true);win.addEventListener('pointerup',interactionEnd,true);win.addEventListener('pointercancel',interactionEnd,true);win.addEventListener('keydown',interactionStart,true);win.addEventListener('keyup',interactionEnd,true);win.addEventListener('blur',interactionEnd);win.addEventListener('wheel',wheel,{capture:true,passive:true});
 post('ready',{identity:options.identity});publishSnapshot(options.getSnapshot());
 return {publishSnapshot,publishPose,dispose(){disposed=true;win.removeEventListener('message',receive);win.removeEventListener('pointerdown',interactionStart,true);win.removeEventListener('pointerup',interactionEnd,true);win.removeEventListener('pointercancel',interactionEnd,true);win.removeEventListener('keydown',interactionStart,true);win.removeEventListener('keyup',interactionEnd,true);win.removeEventListener('blur',interactionEnd);win.removeEventListener('wheel',wheel,true);if(poseFrame)win.cancelAnimationFrame(poseFrame);if(interactionTimer)win.clearTimeout(interactionTimer);}};
}
