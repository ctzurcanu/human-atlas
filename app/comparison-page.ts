import './comparison-page.css';
export const COMPARE_CHANNEL='atlas-compare/v1';
export type CompareSide='left'|'right';
type ChildType='ready'|'snapshot'|'pose'|'interaction'|'view';
export type CompareEnvelope={channel:typeof COMPARE_CHANNEL;session:string;side:CompareSide;instance:string;type:ChildType;seq:number;[key:string]:unknown};
type ParentMessage={channel:typeof COMPARE_CHANNEL;session:string;side:CompareSide;instance?:string;type:'apply'|'request-snapshot'|'request-view';kind?:'state'|'pose';snapshot?:unknown;pose?:unknown;revision?:number;requestId?:string};
const other=(side:CompareSide):CompareSide=>side==='left'?'right':'left';
const object=(value:unknown):value is Record<string,unknown>=>!!value&&typeof value==='object'&&!Array.isArray(value);

/** A comparison embeds the atlas viewer, never an arbitrary page or remote asset. */
export function normalizeViewerUrl(value:string,baseHref:string):string{
 const base=new URL(baseHref),viewerBase=new URL('.',base);
 if(!value.trim())throw Error('Enter an atlas viewer link.');
 let url:URL;try{url=new URL(value.trim(),viewerBase);}catch{throw Error('Enter a valid atlas viewer link.');}
 if(!/^https?:$/.test(url.protocol)||url.origin!==base.origin||url.username||url.password)throw Error('Use an atlas viewer link from this site.');
 if(url.pathname!==viewerBase.pathname&&url.pathname!==viewerBase.pathname+'index.html')throw Error('Use the atlas viewer page, not a document or comparison page.');
 const source=url.searchParams.get('atlasSource');
 if(source){let asset:URL;try{asset=new URL(source,url);}catch{throw Error('Use a local atlas catalogue.');}if(asset.origin!==base.origin||!/^https?:$/.test(asset.protocol)||asset.username||asset.password||!asset.pathname.endsWith('.json'))throw Error('Use a local JSON atlas catalogue.');}
 url.searchParams.delete('compareSession');url.searchParams.delete('compareSide');
 return url.href;
}

export function comparisonFrameUrl(viewer:string,session:string,side:CompareSide):string{
 const url=new URL(viewer);url.searchParams.set('compareSession',session);url.searchParams.set('compareSide',side);url.searchParams.set('ui','model,systems,explode,search,study,sections,details,open,download');return url.href;
}

export function makeComparisonLink(first:string,second:string,baseHref:string):string{
 const url=new URL('compare.html',new URL('.',baseHref));url.searchParams.set('first',normalizeViewerUrl(first,baseHref));url.searchParams.set('second',normalizeViewerUrl(second,baseHref));return url.href;
}

export function isCompareEnvelope(data:unknown,session:string,side:CompareSide):data is CompareEnvelope{
 return object(data)&&data.channel===COMPARE_CHANNEL&&data.session===session&&data.side===side&&typeof data.instance==='string'&&/^[\w-]{1,128}$/.test(data.instance)&&typeof data.seq==='number'&&Number.isSafeInteger(data.seq)&&data.seq>0&&['ready','snapshot','pose','interaction','view'].includes(String(data.type));
}

export function acceptCompareEvent(event:{origin:string;source:unknown;data:unknown},expectedWindow:unknown,origin:string,session:string,side:CompareSide):boolean{
 return expectedWindow!==null&&expectedWindow!==undefined&&event.source===expectedWindow&&event.origin===origin&&isCompareEnvelope(event.data,session,side);
}

export function comparisonSourceLabel(model:string,source:string):string{
 if(!source)return model;
 let pieces:string[];try{pieces=new URL(source,'https://atlas.invalid/').pathname.split('/').filter(Boolean).map(p=>decodeURIComponent(p));}catch{return model;}
 const file=pieces.at(-1)??'',folder=pieces.at(-2)??'';
 const name=folder&&!['assets','models','local-models','.local-models'].includes(folder)?folder:file.replace(/\.(?:json|gltf|glb)$/i,'');
 if(!name)return model;
 const version=/(?:^|[-_])(v\d+)(?:$|[-_])/i.exec(name)?.[1];
 return version?`${version} · ${model} · ${name}`:`${name} · ${model}`;
}

/** The latest real gesture leads; programmatic follower echoes never lead. */
export class ComparisonRelay{
 readonly ready:Record<CompareSide,boolean>={left:false,right:false};
 readonly lastSeq:Record<CompareSide,number>={left:0,right:0};
 readonly instances:Record<CompareSide,string|null>={left:null,right:null};
 readonly identities:Record<CompareSide,{model:string;source:string}|null>={left:null,right:null};
 leader:CompareSide|null=null;
 revision=0;
 private latest:Record<CompareSide,{snapshot?:unknown;pose?:unknown}>={left:{},right:{}};
 private retired:Record<CompareSide,Set<string>>={left:new Set(),right:new Set()};
 constructor(readonly session:string,private send:(side:CompareSide,message:ParentMessage)=>void){}
 private post(side:CompareSide,message:ParentMessage){this.send(side,{...message,instance:this.instances[side]??undefined});}
 requestSnapshot(side:CompareSide){this.post(side,{channel:COMPARE_CHANNEL,session:this.session,side,type:'request-snapshot'});}
 requestView(side:CompareSide,requestId:string){if(!/^[\w-]{1,128}$/.test(requestId))throw Error('Invalid view request.');this.post(side,{channel:COMPARE_CHANNEL,session:this.session,side,type:'request-view',requestId});}
 resetSide(side:CompareSide){
  const hadDocument=this.instances[side]!==null||this.ready[side];
  if(this.instances[side])this.retired[side].add(this.instances[side]!);
  this.instances[side]=null;this.identities[side]=null;this.ready[side]=false;this.lastSeq[side]=0;this.latest[side]={};
  if(hadDocument&&this.leader===side)this.leader=this.ready[other(side)]?other(side):null;
 }
 private apply(side:CompareSide,kind:'state'|'pose',value:unknown){
  this.latest[side][kind==='state'?'snapshot':'pose']=value;
  this.post(side,{channel:COMPARE_CHANNEL,session:this.session,side,type:'apply',kind,[kind==='state'?'snapshot':'pose']:value,revision:++this.revision});
 }
 receive(side:CompareSide,data:unknown):boolean{
  if(!isCompareEnvelope(data,this.session,side))return false;
  if(data.type==='ready'&&data.instance!==this.instances[side]){
   if(this.retired[side].has(data.instance))return false;
   if(this.instances[side])this.resetSide(side);
   this.instances[side]=data.instance;
  }
  if(data.instance!==this.instances[side]||data.seq<=this.lastSeq[side])return false;
  if(data.type==='snapshot'&&!object(data.snapshot)||data.type==='pose'&&!object(data.pose))return false;
  if(data.type==='interaction'&&data.phase!=='start'&&data.phase!=='end')return false;
  if(data.type==='view'&&(typeof data.requestId!=='string'||!/^[\w-]{1,128}$/.test(data.requestId)||typeof data.url!=='string'))return false;
  this.lastSeq[side]=data.seq;
  if(data.type==='view')return true;
  if(data.type==='interaction'){if(data.phase==='start')this.leader=side;return true;}
  if(data.type==='ready'){
   const identity=object(data.identity)&&typeof data.identity.model==='string'?{model:data.identity.model,source:typeof data.identity.source==='string'?data.identity.source:''}:null;
   const previous=this.identities[side],changed=!!(previous&&identity&&(previous.model!==identity.model||previous.source!==identity.source));
   const alreadyReady=this.ready[side]&&!changed;
   if(changed){this.latest[side]={};if(this.leader===side)this.leader=this.ready[other(side)]?other(side):null;}
   this.identities[side]=identity;
   this.ready[side]=true;
   if(alreadyReady)return true;
   if(this.ready.left&&this.ready.right){this.leader??='left';const cached=this.latest[this.leader],target=other(this.leader);if(cached.snapshot)this.apply(target,'state',cached.snapshot);if(cached.pose)this.apply(target,'pose',cached.pose);this.requestSnapshot(this.leader);}
   return true;
  }
  const kind=data.type==='snapshot'?'state':'pose',value=data.type==='snapshot'?data.snapshot:data.pose;
  this.latest[side][data.type==='snapshot'?'snapshot':'pose']=value;
  if(this.leader===null)this.leader='left';
  if(this.leader!==side)return false;
  const target=other(side);if(this.ready[target])this.apply(target,kind,value);
  return true;
 }
}

function freshToken():string{return typeof crypto.randomUUID==='function'?crypto.randomUUID():Array.from(crypto.getRandomValues(new Uint8Array(16)),v=>v.toString(16).padStart(2,'0')).join('');}

/** Obtain saved views from the live camera/state, rather than stale iframe URLs. */
export class ComparisonViewRequests{
 private pending=new Map<string,{side:CompareSide;timer:ReturnType<typeof setTimeout>;resolve:(url:string)=>void;reject:(reason:Error)=>void}>();
 constructor(readonly session:string,private baseHref:string,private send:(side:CompareSide,message:ParentMessage)=>void,private timeoutMs=4000){}
 get pendingCount(){return this.pending.size;}
 request(side:CompareSide):Promise<string>{
  return new Promise((resolve,reject)=>{const requestId=freshToken();const timer=setTimeout(()=>{this.pending.delete(requestId);reject(Error('A viewer did not return its current view. Try copying again.'));},this.timeoutMs);this.pending.set(requestId,{side,timer,resolve,reject});try{this.send(side,{channel:COMPARE_CHANNEL,session:this.session,side,type:'request-view',requestId});}catch(reason){this.pending.delete(requestId);clearTimeout(timer);reject(reason);}});
 }
 requestPair():Promise<[string,string]>{return Promise.all([this.request('left'),this.request('right')]).catch(reason=>{this.cancel('The paired view request did not complete.');throw reason;});}
 receive(side:CompareSide,data:unknown):boolean{
  if(!isCompareEnvelope(data,this.session,side)||data.type!=='view'||typeof data.requestId!=='string'||typeof data.url!=='string')return false;
  const request=this.pending.get(data.requestId);if(!request||request.side!==side)return false;
  this.pending.delete(data.requestId);clearTimeout(request.timer);
  try{request.resolve(normalizeViewerUrl(data.url,this.baseHref));}catch{request.reject(Error('The viewer returned an invalid atlas link.'));}
  return true;
 }
 cancel(message='The views changed. Copy the new pair when it finishes loading.'){
  for(const request of this.pending.values()){clearTimeout(request.timer);request.reject(Error(message));}this.pending.clear();
 }
}

/** Continuously save live paired views, coalescing a drag without losing its end. */
export class ComparisonUrlSync{
 private timer:ReturnType<typeof setTimeout>|undefined;
 private running=false;
 private dirty=false;
 private generation=0;
 constructor(private readPair:()=>Promise<[string,string]>,private writePair:(pair:[string,string])=>void,private delayMs=80){}
 changed(){this.dirty=true;this.schedule();}
 flush(){this.dirty=true;if(this.timer!==undefined){clearTimeout(this.timer);this.timer=undefined;}void this.capture();}
 cancel(){this.generation++;this.dirty=false;if(this.timer!==undefined){clearTimeout(this.timer);this.timer=undefined;}}
 private schedule(){if(!this.running&&this.timer===undefined)this.timer=setTimeout(()=>{this.timer=undefined;void this.capture();},this.delayMs);}
 private async capture(){
  if(this.running||!this.dirty)return;
  this.running=true;this.dirty=false;const generation=this.generation;
  try{const pair=await this.readPair();if(generation===this.generation)this.writePair(pair);}catch{/* A navigating frame will re-announce readiness. */}
  finally{this.running=false;if(this.dirty)this.schedule();}
 }
}

// Local saved views retain the exact 38-vein roster and 14 translucent context meshes.
export const LOCAL_DEFAULT_VIEWER_URLS={"first":"/?model=local-ta98&p=6rv.1i5aj2l.9t7hdf&s=679.66k.676.6qu.6qv.6qw.68c.67x.65r.64l.68e.68a.681.64h.67y.65t.65v.682.687.678.66l.677.6r0.6r1.6r2.68b.683.65s.64m.68d.689.680.64i.684.65u.65w.67z.688&h=0-19i.19l-1a2.1a5-1k1.1k4-1x1.1x4-2up.2us-36t.36w-37f.37i-64g.64j.64k.64n-65q.65x-66j.66m-675.67a-67w.685.686.68f-6qt.6qx-6qz.6r3-6ru&layers=venous%2Cmuscular%2Cnervous&skin=0&context=0.18&isolate=0&labels=0&view=back&frame=0.5%2C0.509998%2C0.0001%2C0.5%2C0.5%2C0.5%2C0.5%2C0.43987%2C0.574468%2C0.5%2C1%2C0.5&tree=systems&atlasSource=%2Flocal-models%2Fgeometry-workbench%2Fregistered-collection-male-circulation-v1%2Focular-IOV-pterygoid-diagnostic-installed-v2938%2Fcirculation-atlas.json","second":"/?model=local-ta98&p=6rv.1i5aj2l.9t7hdf&s=679.66k.676.6qu.6qv.6qw.68c.67x.65r.64l.68e.68a.681.64h.67y.65t.65v.682.687.678.66l.677.6r0.6r1.6r2.68b.683.65s.64m.68d.689.680.64i.684.65u.65w.67z.688&h=0-19i.19l-1a2.1a5-1k1.1k4-1x1.1x4-2up.2us-36t.36w-37f.37i-64g.64j.64k.64n-65q.65x-66j.66m-675.67a-67w.685.686.68f-6qt.6qx-6qz.6r3-6ru&layers=venous%2Cmuscular%2Cnervous&skin=0&context=0.18&isolate=0&labels=0&view=back&frame=0.5%2C0.509998%2C0.0001%2C0.5%2C0.5%2C0.5%2C0.5%2C0.43987%2C0.574468%2C0.5%2C1%2C0.5&tree=systems&atlasSource=%2Flocal-models%2Fgeometry-workbench%2Fregistered-collection-male-circulation-v1%2Focular-vortex-uninstalled-rotation-preview-v2941%2Fcirculation-atlas.json"};
export const DEFAULT_VIEWER_URLS=import.meta.env?.DEV?LOCAL_DEFAULT_VIEWER_URLS:{first:'./?model=male-detail',second:'./?model=female-ta98'};

function mountComparison(){
 const frames:Record<CompareSide,HTMLIFrameElement>={left:document.querySelector<HTMLIFrameElement>('#first-view')!,right:document.querySelector<HTMLIFrameElement>('#second-view')!};
 let session='',relay:ComparisonRelay,viewRequests:ComparisonViewRequests|undefined;
 let documents:Record<CompareSide,Document|null>={left:null,right:null};
 const urlSync=new ComparisonUrlSync(()=>relay.ready.left&&relay.ready.right?viewRequests!.requestPair():Promise.reject(Error('View loading')),([first,second])=>{
  if(relay.ready.left&&relay.ready.right){const url=makeComparisonLink(first,second,location.href);if(url!==location.href)history.replaceState(null,'',url);}
 });
 function rememberPair(immediate=false){if(relay.ready.left&&relay.ready.right){if(immediate)urlSync.flush();else urlSync.changed();}}
 function bindDocument(side:CompareSide){
  let current:Document|null;try{current=frames[side].contentDocument;}catch{return;}
  if(current&&current!==documents[side]){documents[side]=current;relay.resetSide(side);urlSync.cancel();viewRequests?.cancel();}
 }
 function loadPair(updateHistory=true){
  const params=new URL(location.href).searchParams;
  const sources=[params.get('first')??DEFAULT_VIEWER_URLS.first,params.get('second')??DEFAULT_VIEWER_URLS.second];
  const defaults=[DEFAULT_VIEWER_URLS.first,DEFAULT_VIEWER_URLS.second];
  const urls=sources.map((value,index)=>{try{return normalizeViewerUrl(value,location.href);}catch{return normalizeViewerUrl(defaults[index],location.href);}});
  urlSync.cancel();viewRequests?.cancel();session=freshToken();documents={left:null,right:null};
  const send=(side:CompareSide,message:ParentMessage)=>frames[side].contentWindow?.postMessage({...message,instance:relay.instances[side]??undefined},location.origin);
  relay=new ComparisonRelay(session,send);viewRequests=new ComparisonViewRequests(session,location.href,send);
  frames.left.src=comparisonFrameUrl(urls[0],session,'left');frames.right.src=comparisonFrameUrl(urls[1],session,'right');
  if(updateHistory)history.replaceState(null,'',makeComparisonLink(urls[0],urls[1],location.href));
 }
 for(const side of ['left','right']as const)frames[side].addEventListener('load',()=>{bindDocument(side);relay.requestSnapshot(side);});
 window.addEventListener('message',(event:MessageEvent)=>{
  for(const side of ['left','right']as const){
   if(!acceptCompareEvent(event,frames[side].contentWindow,location.origin,session,side))continue;
   bindDocument(side);
   const data=event.data as CompareEnvelope,wasReady=relay.ready[side],previousInstance=relay.instances[side],previousIdentity=relay.identities[side];
   if(!relay.receive(side,data))return;
   if(data.type==='view'){viewRequests?.receive(side,data);return;}
   if(data.type==='ready'){
    const identity=relay.identities[side];
    frames[side].title=identity?comparisonSourceLabel(identity.model,identity.source):'Anatomy viewport';
    if(!wasReady||previousInstance!==data.instance||previousIdentity?.model!==identity?.model||previousIdentity?.source!==identity?.source)rememberPair();
   }
   if(data.type==='pose'||data.type==='snapshot')rememberPair();
   if(data.type==='interaction'&&data.phase==='end')rememberPair(true);
   return;
  }
 });
 window.addEventListener('popstate',()=>loadPair(false));
 loadPair();
}

if(typeof window!=='undefined'&&typeof document!=='undefined'&&document.querySelector('#first-view'))mountComparison();
