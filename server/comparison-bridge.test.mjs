import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
import {gzipSync,strToU8} from 'fflate';

const built=await build({entryPoints:['app/comparison-bridge.ts'],bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const {COMPARISON_CHANNEL,comparisonSession,comparisonControlVisible,comparisonModelUrl,comparisonMessage,describeComparisonSelection,resolveComparisonSelection,comparisonSnapshot,applyComparisonSnapshot,validComparisonSnapshot,validComparisonPose,createComparisonBridge}=api;
const SESSION='comparison-session-1234',ORIGIN='https://atlas.example';
const camera=[.5,.5,.9,.5,.5,.5,.5,.5,.5,.5,1,.5];
const pose={model:'male-model',camera,up:[0,1,0],rotation:[0,0,0,1],world:{position:[0,1,2],target:[0,1,0],up:[0,1,0],quaternion:[0,0,0,1],fov:40,zoom:1,near:.001,far:100,view:{enabled:true,fullWidth:1200,fullHeight:800,offsetX:80,offsetY:20,width:1120,height:780}}};

function atlas(prefix,{alternatives=false}={}){
 const definitions=[
  [`${prefix}-left`,'Vorticose veins (left)','A12.3.06.106'],
  [`${prefix}-right`,'Vorticose veins (right)','A12.3.06.106'],
  ['shared-local','Shared explicit geometry',null],
  [`${prefix}-uncoded`,'Unqualified source patch',null],
  [`${prefix}-eye`,'Eye (left)',null,'FMA54450'],
 ];
 if(alternatives)definitions.push([`${prefix}-old`,'Vorticose veins (left)','A12.3.06.106']);
 const parts=definitions.map(([id,name,term,concept])=>({id,name,conceptId:concept??`${id}-concept`,system:id.endsWith('-eye')?'sensory':'venous',...(term?{ta98Term:term}:{}),...(id.endsWith('-old')?{suppressed:true}:{})}));
 const concepts=parts.map(p=>({id:p.conceptId,name:p.name,elements:[p.id],...(p.ta98Term?{ta98Term:p.ta98Term,ta98Kind:'exact'}:{})}));
 return {parts,concepts,defaultHidden:alternatives?[`${prefix}-old`]:[]};
}
function state(selected=[]){return {model:'female-model',source:'/models/female.json',selected,hidden:[],visible:['venous'],depthHidden:[],contextOpacity:.2,skinOpacity:.1,labels:true,isolate:false,camera:[...camera],focus:17,reset:4,view:'right',region:'head',explode:.4,section:{enabled:true,axis:'coronal',position:.37,flip:true},sections:[{enabled:true,axis:'coronal',position:.37,flip:true}],activeSection:0,panel:'details',hierarchy:'systems'};}
const source=atlas('male'),target=atlas('female',{alternatives:true});
function sourceSnapshot(){return comparisonSnapshot(source,{...state(['male-left']),hidden:['male-right'],visible:['venous','not-in-target'],labels:false,isolate:true,contextOpacity:.15,skinOpacity:.3,depth:.6,depthHidden:['muscles']},'guest:ta98');}

class FakeWindow{
 constructor(){this.location={origin:ORIGIN,href:ORIGIN+'/?model=female'};this.sent=[];this.listeners=new Map;this.frames=new Map;this.timers=new Map;this.next=0;this.parent={postMessage:(data,origin)=>this.sent.push({data:structuredClone(data),origin})};}
 addEventListener(type,fn){if(!this.listeners.has(type))this.listeners.set(type,new Set);this.listeners.get(type).add(fn);}
 removeEventListener(type,fn){this.listeners.get(type)?.delete(fn);}
 fire(type,event={}){for(const fn of [...this.listeners.get(type)??[]])fn(event);}
 message(data,overrides={}){const scoped={...(['apply','request-view'].includes(data.type)?{instance:this.sent[0]?.data.instance}:{}),...data};this.fire('message',{origin:ORIGIN,source:this.parent,data:scoped,...overrides});}
 requestAnimationFrame(fn){const id=++this.next;this.frames.set(id,fn);return id;}
 cancelAnimationFrame(id){this.frames.delete(id);}
 flushFrames(){const callbacks=[...this.frames.values()];this.frames.clear();for(const fn of callbacks)fn(1);}
 setTimeout(fn){const id=++this.next;this.timers.set(id,fn);return id;}
 clearTimeout(id){this.timers.delete(id);}
}
function bridgeFixture(options={}){
 const win=new FakeWindow,appliedStates=[],appliedPoses=[];let localState={...state(),hidden:['female-old']},liveUrl=ORIGIN+'/?model=female&view=front',viewRequests=0;
 const snapshot=()=>comparisonSnapshot(target,localState,'systems');
 const bridge=createComparisonBridge({window:win,session:SESSION,side:'right',identity:{model:'female-model',source:'/models/female.json'},getSnapshot:snapshot,getViewUrl:()=>{viewRequests++;return liveUrl;},...options,applySnapshot:value=>{appliedStates.push(value);localState=applyComparisonSnapshot(target,localState,value);return snapshot();},applyPose:value=>appliedPoses.push(value)});
 return {win,bridge,appliedStates,appliedPoses,snapshot,getState:()=>localState,setLiveUrl:value=>{liveUrl=value;},viewRequests:()=>viewRequests};
}
function incoming(kind,value,revision=1){return {channel:COMPARISON_CHANNEL,session:SESSION,side:'right',type:'apply',kind,revision,[kind==='state'?'snapshot':'pose']:value};}

test('Comparison control is enabled only for a valid framed session and side',()=>{
 assert.deepEqual(comparisonSession(new URLSearchParams({compareSession:SESSION,compareSide:'left'}),true),{session:SESSION,side:'left'});
 for(const [query,framed]of [[new URLSearchParams({compareSession:SESSION,compareSide:'right'}),false],[new URLSearchParams({compareSession:'short',compareSide:'left'}),true],[new URLSearchParams({compareSession:SESSION,compareSide:'both'}),true],[new URLSearchParams(),true],[new URLSearchParams({compareSession:'bad token!',compareSide:'left'}),true]])assert.equal(comparisonSession(query,framed),null);
});

test('A comparison model choice retains live presentation and transport while clearing every old-source route and index',()=>{
 const current=new URL(ORIGIN+'/?model=local-ta98&atlasSource=%2Flocal-models%2Fmale.json&blood=biodigital-review&cardiac=review&bowel=prototype&hepatic=repair&biliary=folds&venous=confluence&caps=shared-v28&cricoid=partition-v1&cricotracheal=review&delivery=1&review=alternatives&cardiacRefStage=old&p=old-fingerprint&s=packed&h=old-hidden&select=male-mesh&hide=male-mesh&bio=male-term&frame=old-frame&camera=old-camera&focus=1&isolate=1&ui=model,systems,sections&tree=systems&layers=venous,arterial&skin=0&context=0.2&labels=0&explode=0.4&dp=0.3&cut=oblique&slice=0.6&azimuth=30&elevation=-15&flip=1#old-saved-model');
 const before=current.href,next=new URL(comparisonModelUrl(current.href,'local-female-ta98',{session:SESSION,side:'right'}));assert.equal(current.href,before);assert.equal(next.origin,current.origin);assert.equal(next.searchParams.get('model'),'local-female-ta98');assert.equal(next.searchParams.get('compareSession'),SESSION);assert.equal(next.searchParams.get('compareSide'),'right');assert.equal(next.hash,'');
 for(const key of ['atlasSource','blood','cardiac','bowel','hepatic','biliary','venous','caps','cricoid','cricotracheal','delivery','review','cardiacRefStage','p','s','h','select','hide','bio','frame','camera','focus','isolate'])assert.equal(next.searchParams.get(key),null,key+' is scoped to the old model');for(const key of ['ui','tree','layers','skin','context','labels','explode','dp','cut','slice','azimuth','elevation','flip'])assert.equal(next.searchParams.get(key),current.searchParams.get(key));assert.throws(()=>comparisonModelUrl(current.href,'not a model',{session:SESSION,side:'right'}));
});

test('Changing a model strips old model routing hidden inside a compressed saved view as well as the readable query',()=>{
 const saved=new URLSearchParams({model:'local-ta98',atlasSource:'/local-models/old-male.json',blood:'biodigital-review',p:'old-index',s:'old-selection',frame:'old-frame',layers:'venous',skin:'0',context:'0.18',cut:'coronal',slice:'0.61',ui:'model,systems,details'}),encoded=Buffer.from(gzipSync(strToU8(saved.toString()))).toString('base64url'),url=ORIGIN+'/?compareSession='+SESSION+'&compareSide=left#atlas-view=1.'+encoded,next=new URL(comparisonModelUrl(url,'female',{session:SESSION,side:'left'}));
 for(const key of ['atlasSource','blood','p','s','frame'])assert.equal(next.searchParams.get(key),null);assert.equal(next.searchParams.get('model'),'female');assert.equal(next.searchParams.get('layers'),'venous');assert.equal(next.searchParams.get('context'),'0.18');assert.equal(next.searchParams.get('cut'),'coronal');assert.equal(next.searchParams.get('ui'),'model,systems,details');assert.equal(next.hash,'');
});

test('Child message guard independently rejects wrong origin, parent window, session, side and channel',()=>{
 const parent={},expected={origin:ORIGIN,parent,session:SESSION,side:'right'},data=incoming('state',sourceSnapshot());
 assert.equal(comparisonMessage({origin:ORIGIN,source:parent,data},expected),data);
 for(const event of [{origin:'https://other.example',source:parent,data},{origin:ORIGIN,source:{},data},...['session','side','channel'].map(key=>({origin:ORIGIN,source:parent,data:{...data,[key]:'wrong'}})),{origin:ORIGIN,source:parent,data:[]}])assert.equal(comparisonMessage(event,expected),null);
});

test('Exact local identities and portable TA98 selection map across model-specific IDs without changing laterality',()=>{
 const described=describeComparisonSelection(source,['male-left']);assert.equal(described.length,1);assert.equal(described[0].ta98,'A12.3.06.106');assert.equal(described[0].side,'left');
 assert.deepEqual(resolveComparisonSelection(target,described),['female-left']);
 assert.deepEqual(resolveComparisonSelection(target,describeComparisonSelection(source,['male-right'])),['female-right']);
 assert.deepEqual(resolveComparisonSelection(target,describeComparisonSelection(source,['shared-local'])),['shared-local']);
 assert.deepEqual(resolveComparisonSelection(target,[...described,...described]),['female-left']);
});

test('FMA identity is portable even when the originating local concept ID differs',()=>{
 const local=describeComparisonSelection(target,['female-eye'])[0];assert.equal(local.fma,'FMA:54450');
 const remote={id:'other-eye-source',conceptId:'other-eye-concept',name:'Eye (left)',fma:local.fma,side:'left'};
 assert.deepEqual(resolveComparisonSelection(target,[remote]),['female-eye']);
});

test('Names alone and parent-term placements do not invent a semantic selection',()=>{
 assert.deepEqual(resolveComparisonSelection(target,describeComparisonSelection(source,['male-uncoded'])),[]);
 const parent={parts:[{id:'parent-part',name:'Vorticose veins (left)',conceptId:'parent-concept',system:'venous',ta98Term:'A12.3.06.106'}],concepts:[{id:'parent-concept',name:'Incomplete source',elements:['parent-part'],ta98Term:'A12.3.06.106',ta98Kind:'parent'}]};
 const described=describeComparisonSelection(parent,['parent-part']);assert.equal(described[0].ta98,undefined);assert.deepEqual(resolveComparisonSelection(target,described),[]);
});

test('Semantic selection uses the target active source instead of exposing its suppressed alternative',()=>{
 assert.deepEqual(resolveComparisonSelection(target,describeComparisonSelection(source,['male-left'])),['female-left']);
 const explicit=describeComparisonSelection(target,['female-old']);assert.deepEqual(resolveComparisonSelection(target,explicit),['female-old'],'An exact explicitly chosen local ID remains intentional');
});

test('State synchronization preserves model/source/camera/sections/focus and both input objects',()=>{
 const before=state(['female-right']),snapshot={...sourceSnapshot(),model:'remote-model',source:'/untrusted-catalogue.json',camera:[0,0,0],focus:99,reset:99,view:'back',section:{axis:'sagittal',position:.8,flip:false}},saved=structuredClone(before),remote=structuredClone(snapshot),next=applyComparisonSnapshot(target,before,snapshot);
 assert.deepEqual(next.selected,['female-left']);assert.deepEqual(new Set(next.hidden),new Set(['female-old','female-right']));assert.deepEqual(next.visible,['venous']);assert.equal(next.contextOpacity,.15);assert.equal(next.skinOpacity,.3);assert.equal(next.labels,false);assert.equal(next.depth,.6);
 for(const key of ['model','source','camera','focus','reset','view','region','explode','section','sections','activeSection','panel','hierarchy'])assert.deepEqual(next[key],saved[key],`${key} is local identity or camera/navigation intent`);
 assert.deepEqual(before,saved);assert.deepEqual(snapshot,remote);
});

test('Source default hidden guides are not broadcast as manual hides and target defaults remain local',()=>{
 const from=atlas('from',{alternatives:true}),to=atlas('to',{alternatives:true}),snap=comparisonSnapshot(from,{...state(['from-left']),hidden:['from-old','from-right']},'systems');
 assert.deepEqual(snap.hidden.map(x=>x.id),['from-right']);
 const next=applyComparisonSnapshot(to,{...state(),hidden:['to-old']},snap);assert.deepEqual(new Set(next.hidden),new Set(['to-old','to-right']));assert.ok(!next.hidden.includes('to-left'));
 const unmatched={...snap,selection:[{id:'missing',conceptId:'missing',name:'No counterpart'}],isolate:true};assert.equal(applyComparisonSnapshot(to,state(),unmatched).isolate,false,'Unmatched anatomy does not isolate an empty viewport');
});

test('Shared Explode and Sections controls synchronize presentation without replacing either camera or catalogue',()=>{
 const sections=[{enabled:true,axis:'oblique',position:.61,flip:false,azimuth:-47,elevation:28},{enabled:false,axis:'sagittal',position:.23,flip:true}],remoteState={...state(['male-left']),explode:.83,sections,activeSection:1,section:sections[1]},snapshot=comparisonSnapshot(source,remoteState,'systems'),before={...state(['female-right']),rotate:true},saved=structuredClone(before),next=applyComparisonSnapshot(target,before,snapshot);
 assert.equal(snapshot.explode,.83);assert.deepEqual(snapshot.sections,sections);assert.equal(snapshot.activeSection,1);assert.equal(next.explode,.83);assert.equal(next.rotate,false,'Explode pauses rotation in both panes');assert.deepEqual(next.sections,sections);assert.deepEqual(next.section,sections[1]);assert.equal(next.activeSection,1);
 for(const key of ['model','source','camera','focus','reset','view','region','panel'])assert.deepEqual(next[key],saved[key]);assert.notEqual(snapshot.sections,sections);assert.notEqual(snapshot.sections[0],sections[0]);assert.notEqual(next.sections,snapshot.sections);assert.notEqual(next.sections[0],snapshot.sections[0]);assert.deepEqual(before,saved);
 const imploded=applyComparisonSnapshot(target,next,{...snapshot,explode:0});assert.equal(imploded.explode,0);assert.deepEqual(imploded.camera,saved.camera);assert.equal(imploded.reset,saved.reset);
});

test('Legacy comparison state leaves local explode/sections intact and serialization normalizes a single cut',()=>{
 const snapshot=sourceSnapshot();delete snapshot.explode;delete snapshot.sections;delete snapshot.activeSection;const before=state(),next=applyComparisonSnapshot(target,before,snapshot);for(const key of ['explode','sections','activeSection','section'])assert.deepEqual(next[key],before[key]);
 const legacy={...state(),sections:undefined,activeSection:99,section:{enabled:true,axis:'axial',position:.42,flip:false}};const normalized=comparisonSnapshot(source,legacy,'systems');assert.deepEqual(normalized.sections,[legacy.section]);assert.equal(normalized.activeSection,0);assert.ok(validComparisonSnapshot(normalized));
});

test('Shared state and pose validation reject malformed or nonfinite payloads',()=>{
 const snap=sourceSnapshot();assert.ok(validComparisonSnapshot(snap));assert.ok(validComparisonPose(pose));
 for(const bad of [{...snap,context:NaN},{...snap,skinOpacity:2},{...snap,selection:[{id:'x',conceptId:'x',name:'x',side:'both'}]},{...snap,systems:'venous'},null])assert.equal(validComparisonSnapshot(bad),false);
 for(const bad of [{...snap,explode:-.1},{...snap,explode:1.1},{...snap,explode:NaN},{...snap,sections:[]},{...snap,sections:Array(3).fill(snap.sections[0])},{...snap,sections:[{...snap.sections[0],enabled:'true'}]},{...snap,sections:[{...snap.sections[0],axis:'sideways'}]},{...snap,sections:[{...snap.sections[0],position:Infinity}]},{...snap,sections:[{...snap.sections[0],azimuth:181}]},{...snap,sections:[{...snap.sections[0],elevation:-91}]},{...snap,activeSection:1},{...snap,activeSection:-1},{...snap,activeSection:.5}])assert.equal(validComparisonSnapshot(bad),false);
 for(const bad of [{...pose,camera:[1,1,1,1,1,1]},{...pose,up:[0,NaN,0]},{...pose,up:[0,0,0]},{...pose,rotation:[0,0,1]},{...pose,world:{...pose.world,position:[0,Infinity,0]}},{...pose,world:{...pose.world,zoom:0}},{...pose,world:{...pose.world,fov:180}},{...pose,world:{...pose.world,far:.0001}},{...pose,world:{...pose.world,quaternion:[0,0,0,0]}},{...pose,world:{...pose.world,target:[...pose.world.position]}},null])assert.equal(validComparisonPose(bad),false);
});

test('The bridge starts with scoped identity and rejects unauthorized state changes',()=>{
 const f=bridgeFixture();assert.deepEqual(f.win.sent.map(x=>x.data.type),['ready','snapshot']);assert.deepEqual(f.win.sent[0].data.identity,{model:'female-model',source:'/models/female.json'});assert.ok(f.win.sent.every(x=>x.origin===ORIGIN&&x.data.session===SESSION&&x.data.side==='right'&&x.data.instance===f.win.sent[0].data.instance));assert.match(f.win.sent[0].data.instance,/^[A-Za-z0-9_-]{8,128}$/);
 const data=incoming('state',sourceSnapshot());for(const override of [{origin:'https://other.example'},{source:{}},{data:{...data,session:'other-session'}},{data:{...data,side:'left'}}])f.win.message(data,override);
 assert.equal(f.appliedStates.length,0);f.bridge.dispose();
});

test('Each document instance rejects old apply/copy packets while an initial snapshot probe renews its ready identity',()=>{
 const first=bridgeFixture(),second=bridgeFixture();assert.notEqual(first.win.sent[0].data.instance,second.win.sent[0].data.instance);const old=first.win.sent[0].data.instance;
 second.win.message({...incoming('pose',pose,999),instance:old});second.win.message({...incoming('state',sourceSnapshot(),999),instance:old});second.win.message({channel:COMPARISON_CHANNEL,session:SESSION,side:'right',instance:old,type:'request-view',requestId:'old-copy'});assert.equal(second.appliedPoses.length,0);assert.equal(second.appliedStates.length,0);assert.equal(second.viewRequests(),0);
 second.win.fire('message',{origin:ORIGIN,source:second.win.parent,data:incoming('state',sourceSnapshot())});assert.equal(second.appliedStates.length,0,'An apply without a document token is rejected');
 const before=second.win.sent.length;second.win.message({channel:COMPARISON_CHANNEL,session:SESSION,side:'right',type:'request-snapshot'});assert.deepEqual(second.win.sent.slice(before).map(x=>x.data.type),['ready','snapshot']);assert.equal(second.win.sent.at(-2).data.instance,second.win.sent[0].data.instance);second.win.message(incoming('state',sourceSnapshot(),1));assert.equal(second.appliedStates.length,1,'Old document high revisions cannot block the new document');first.bridge.dispose();second.bridge.dispose();
});

test('A reload handshake explicitly captures the latest live follower pose without enabling unsolicited camera echo',()=>{
 let live={...structuredClone(pose),model:'female-model'};const f=bridgeFixture({getPose:()=>live});f.win.message(incoming('pose',pose));const before=f.win.sent.length;live={...live,camera:live.camera.map((v,i)=>i===2?.83:v),world:{...live.world,position:[0,1,3]}};f.bridge.publishPose(live);f.win.flushFrames();assert.equal(f.win.sent.length,before,'A follower still cannot broadcast pose updates by itself');
 f.win.message({channel:COMPARISON_CHANNEL,session:SESSION,side:'right',type:'request-snapshot'});const replies=f.win.sent.slice(before);assert.deepEqual(replies.map(x=>x.data.type),['ready','snapshot','pose']);assert.deepEqual(replies[2].data.pose,live);const count=f.win.sent.length;f.bridge.publishPose({...live,world:{...live.world,position:[0,1,4]}});f.win.flushFrames();assert.equal(f.win.sent.length,count,'The explicit probe does not transfer gesture ownership');f.bridge.dispose();
});

test('Restarting the React bridge effect in one document retains its epoch and monotonic outbound sequence',()=>{
 const f=bridgeFixture(),previous=f.win.sent.at(-1).data;f.bridge.dispose();const bridge=createComparisonBridge({window:f.win,session:SESSION,side:'right',identity:{model:'female-model',source:'/models/female.json'},getSnapshot:f.snapshot,getViewUrl:()=>ORIGIN+'/?model=female',applySnapshot:value=>value,applyPose:()=>{}}),ready=f.win.sent.at(-2).data;assert.equal(ready.instance,previous.instance);assert.ok(ready.seq>previous.seq);assert.ok(f.win.sent.at(-1).data.seq>ready.seq);bridge.dispose();
});

test('Copy requests read the live viewer URL at request time and require the authenticated parent',()=>{
 const f=bridgeFixture(),request={channel:COMPARISON_CHANNEL,session:SESSION,side:'right',type:'request-view',requestId:'copy-first'};
 for(const override of [{origin:'https://other.example'},{source:{}},{data:{...request,session:'wrong-session'}},{data:{...request,side:'left'}},{data:{...request,requestId:'bad request!'}}])f.win.message(request,override);
 assert.equal(f.viewRequests(),0);f.win.message(request);assert.equal(f.viewRequests(),1);assert.equal(f.win.sent.at(-1).data.type,'view');assert.equal(f.win.sent.at(-1).data.requestId,'copy-first');assert.equal(f.win.sent.at(-1).data.url,ORIGIN+'/?model=female&view=front');
 const latest=ORIGIN+'/?model=female&select=just-selected&frame=latest-camera#atlas-view=current';f.setLiveUrl(latest);f.win.message({...request,requestId:'copy-second'});assert.equal(f.viewRequests(),2);assert.equal(f.win.sent.at(-1).data.url,latest);f.bridge.dispose();
});

test('Applied state uses monotonic revisions, suppresses its mapped echo and resumes after user input',()=>{
 const f=bridgeFixture(),data=incoming('state',sourceSnapshot(),4);f.win.message(data);assert.equal(f.appliedStates.length,1);assert.deepEqual(f.getState().selected,['female-left']);
 for(const revision of [4,3,-1])f.win.message({...data,revision});assert.equal(f.appliedStates.length,1);
 const count=f.win.sent.length;f.bridge.publishSnapshot(f.snapshot());assert.equal(f.win.sent.length,count,'The target-local mapped selection must not echo back as a new source update');
 f.win.fire('pointerdown');f.bridge.publishSnapshot({...f.snapshot(),labels:true});assert.equal(f.win.sent.at(-1).data.type,'snapshot');assert.equal(f.win.sent.at(-1).data.snapshot.labels,true);f.bridge.dispose();
});

test('Remote presentation applies both cuts and explosion once without snapshot or camera feedback',()=>{
 const f=bridgeFixture(),snapshot={...sourceSnapshot(),explode:.9,sections:[{enabled:true,axis:'oblique',position:.66,flip:false,azimuth:32,elevation:-15},{enabled:true,axis:'axial',position:.2,flip:true}],activeSection:1},cameraBefore=structuredClone(f.getState().camera),count=f.win.sent.length;
 f.win.message(incoming('state',snapshot,7));assert.equal(f.appliedStates.length,1);assert.equal(f.getState().explode,.9);assert.deepEqual(f.getState().sections,snapshot.sections);assert.deepEqual(f.getState().camera,cameraBefore);assert.equal(f.appliedPoses.length,0);f.bridge.publishSnapshot(f.snapshot());assert.equal(f.win.sent.length,count);f.bridge.dispose();
});

test('Explicit remote pose keeps target identity, blocks feedback, and a real gesture can take camera ownership',()=>{
 const f=bridgeFixture(),original=structuredClone(pose);f.win.message(incoming('state',sourceSnapshot(),20));f.win.message(incoming('pose',pose,1));assert.equal(f.appliedPoses.length,1,'Pose revision is independent of state revision');assert.equal(f.appliedPoses[0].model,'female-model');assert.equal(f.appliedPoses[0].world,undefined,'World units of different models cannot replace the target world frame');assert.deepEqual(f.appliedPoses[0].camera,pose.camera);assert.deepEqual(pose,original);assert.deepEqual(f.getState().camera,camera,'State sync did not move the target camera');
 const count=f.win.sent.length;f.bridge.publishPose(f.appliedPoses[0]);f.win.flushFrames();assert.equal(f.win.sent.length,count,'Follower camera notifications cannot echo the received pose');
 f.win.fire('pointerdown');f.bridge.publishPose(f.appliedPoses[0]);f.win.flushFrames();assert.equal(f.win.sent.at(-1).data.type,'pose');const sent=f.win.sent.length;f.bridge.publishPose(f.appliedPoses[0]);f.win.flushFrames();assert.equal(f.win.sent.length,sent,'Identical owned pose is deduplicated');f.bridge.dispose();
});

test('Same-model comparisons preserve the exact explicit world camera including projection offsets',()=>{
 const f=bridgeFixture(),sameModel={...structuredClone(pose),model:'female-model'};f.win.message(incoming('pose',sameModel));assert.deepEqual(f.appliedPoses,[sameModel]);f.win.message(incoming('pose',sameModel));assert.equal(f.appliedPoses.length,1,'Repeated pose revisions cannot replay camera motion');f.bridge.dispose();
});

test('Male and female metre-body models share actual world POV instead of outlier-normalized bounds',()=>{
 const f=bridgeFixture({identity:{model:'local-female-ta98',source:'/female.json'}}),male={...structuredClone(pose),model:'local-ta98'};f.win.message(incoming('pose',male));assert.equal(f.appliedPoses[0].model,'local-female-ta98');assert.deepEqual(f.appliedPoses[0].world,male.world);assert.ok(api.comparisonSharesWorld('male-detail','female'));assert.equal(api.comparisonSharesWorld('male-detail','cell'),false);assert.equal(api.comparisonSharesWorld('local-ta98','embryo'),false);f.bridge.dispose();
});

test('An in-flight apply cannot steal camera leadership during an active local gesture',()=>{
 const f=bridgeFixture();f.win.fire('pointerdown');f.win.message(incoming('pose',pose,1));f.win.message(incoming('state',sourceSnapshot(),1));assert.equal(f.appliedPoses.length,0);assert.equal(f.appliedStates.length,0);
 const request={channel:COMPARISON_CHANNEL,session:SESSION,side:'right',type:'request-view',requestId:'during-drag'};f.win.message(request);assert.equal(f.viewRequests(),1,'Copy still reads the current view while dragging');f.bridge.publishPose(pose);f.win.flushFrames();assert.equal(f.win.sent.at(-1).data.type,'pose');
 f.win.fire('pointerup');f.win.message(incoming('pose',pose,2));assert.equal(f.appliedPoses.length,1,'A subsequent remote gesture can lead after the local gesture ends');f.bridge.dispose();
});

test('Disposing a viewer bridge removes every listener and pending animation/timer',()=>{
 const f=bridgeFixture();f.win.fire('wheel');f.bridge.publishPose(pose);assert.ok(f.win.frames.size&&f.win.timers.size);f.bridge.dispose();assert.equal(f.win.frames.size,0);assert.equal(f.win.timers.size,0);assert.ok([...f.win.listeners.values()].every(set=>set.size===0));
 const count=f.win.sent.length;f.win.message(incoming('state',sourceSnapshot()));f.bridge.publishSnapshot(sourceSnapshot());f.bridge.publishPose(pose);f.win.flushFrames();assert.equal(f.win.sent.length,count);
});

test('Controls are partitioned between frames while standalone controls remain intact',()=>{
 for(const control of ['systems','explode']){assert.equal(comparisonControlVisible('left',control),true);assert.equal(comparisonControlVisible('right',control),false);}
 for(const control of ['search','study','sections','details','open','download']){assert.equal(comparisonControlVisible('left',control),false);assert.equal(comparisonControlVisible('right',control),true);}
 for(const side of ['left','right']){assert.equal(comparisonControlVisible(side,'model'),true);assert.equal(comparisonControlVisible(side,'camera'),false);}
 for(const control of ['systems','details','model','sections','depth','labels'])assert.equal(comparisonControlVisible(null,control),true);
});
