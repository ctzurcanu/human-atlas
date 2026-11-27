import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile,mkdtemp,mkdir,writeFile,rm,stat,utimes} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {build} from 'esbuild';
import postcss from 'postcss';
import {prepareComparisonPage} from './comparison-page.mjs';

const built=await build({entryPoints:['app/comparison-page.ts'],bundle:true,platform:'node',format:'esm',write:false,loader:{'.css':'empty'}});
const {COMPARE_CHANNEL,normalizeViewerUrl,comparisonFrameUrl,makeComparisonLink,isCompareEnvelope,acceptCompareEvent,ComparisonRelay,ComparisonViewRequests,ComparisonUrlSync,LOCAL_DEFAULT_VIEWER_URLS,DEFAULT_VIEWER_URLS}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const ORIGIN='https://atlas.example',BASE=ORIGIN+'/human-atlas/compare.html',SESSION='comparison-session-1234';
const left=ORIGIN+'/human-atlas/?model=male-detail&atlasSource=%2Fmodels%2Fmale.json&select=left-artery#atlas-view=saved-left';
const right=ORIGIN+'/human-atlas/?model=female&atlasSource=%2Fmodels%2Ffemale.json&select=right-vein#atlas-view=saved-right';
const envelope=(side,type,seq,extra={})=>({channel:COMPARE_CHANNEL,session:SESSION,side,instance:side+'-document-1234',type,seq,...(type==='ready'?{identity:{model:side+'-model',source:'/models/'+side+'.json'}}:{}),...extra});
function relayFixture(){const sent=[];const relay=new ComparisonRelay(SESSION,(side,message)=>sent.push({side,message:structuredClone(message)}));return {sent,relay};}
function ready(relay){relay.receive('left',envelope('left','ready',1));relay.receive('right',envelope('right','ready',1));}

test('Viewer URL normalization preserves each model, catalogue, selection and saved view while dropping transport tokens',()=>{
 const parsed=new URL(normalizeViewerUrl(left,BASE));
 assert.equal(parsed.origin,ORIGIN);assert.equal(parsed.searchParams.get('model'),'male-detail');assert.equal(parsed.searchParams.get('atlasSource'),'/models/male.json');assert.equal(parsed.searchParams.get('select'),'left-artery');assert.ok(parsed.hash.startsWith('#atlas-view=saved-left'));
 const transport=new URL(left);transport.searchParams.set('compareSession','old-session');transport.searchParams.set('compareSide','left');const clean=new URL(normalizeViewerUrl(transport.href,BASE));assert.equal(clean.searchParams.get('compareSession'),null);assert.equal(clean.searchParams.get('compareSide'),null);assert.equal(clean.hash,new URL(left).hash);
 assert.equal(new URL(normalizeViewerUrl('./?model=female',BASE)).pathname,'/human-atlas/');assert.equal(new URL(normalizeViewerUrl('./index.html?model=female',BASE)).pathname,'/human-atlas/index.html');
});

test('Only same-origin atlas viewers and same-origin JSON catalogue references may be embedded',()=>{
 for(const value of ['','https://remote.example/human-atlas/','https://user:password@atlas.example/human-atlas/','javascript:alert(1)','data:text/html,anything','./compare.html','./guide.pdf','/?model=female','./?atlasSource=https%3A%2F%2Fremote.example%2Fmodel.json','./?atlasSource=javascript%3Aalert(1)','./?atlasSource=%2Fmodel.html','./?atlasSource=https%3A%2F%2Fu%3Ap%40atlas.example%2Fmodel.json'])assert.throws(()=>normalizeViewerUrl(value,BASE),undefined,value);
});

test('Comparison links preserve independent viewer identity and frame transport cannot replace either source',()=>{
 const shared=new URL(makeComparisonLink(left,right,BASE));assert.equal(shared.pathname,'/human-atlas/compare.html');assert.equal(shared.searchParams.get('first'),left);assert.equal(shared.searchParams.get('second'),right);
 for(const [url,side]of [[left,'left'],[right,'right']]){const frame=new URL(comparisonFrameUrl(url,SESSION,side)),original=new URL(url);for(const key of ['model','atlasSource','select'])assert.equal(frame.searchParams.get(key),original.searchParams.get(key));assert.equal(frame.hash,original.hash);assert.equal(frame.searchParams.get('compareSession'),SESSION);assert.equal(frame.searchParams.get('compareSide'),side);assert.equal(frame.searchParams.getAll('compareSession').length,1);assert.equal(frame.searchParams.get('ui'),'model,systems,explode,search,study,sections,details,open,download');}
 assert.notEqual(new URL(left).searchParams.get('model'),new URL(right).searchParams.get('model'));assert.notEqual(new URL(left).searchParams.get('atlasSource'),new URL(right).searchParams.get('atlasSource'));
});

test('Parent guard independently checks exact child window, origin, session, side, channel and sequence shape',()=>{
 const child={},data=envelope('left','ready',1),event={origin:ORIGIN,source:child,data};assert.ok(acceptCompareEvent(event,child,ORIGIN,SESSION,'left'));
 for(const invalid of [{...event,origin:'https://remote.example'},{...event,source:{}},...['session','side','channel'].map(key=>({...event,data:{...data,[key]:'wrong'}})),...[0,-1,1.5,NaN,Number.MAX_SAFE_INTEGER+1].map(seq=>({...event,data:{...data,seq}})),{...event,data:{...data,type:'apply'}},{...event,data:[]}])assert.equal(acceptCompareEvent(invalid,child,ORIGIN,SESSION,'left'),false);
 assert.equal(acceptCompareEvent(event,null,ORIGIN,SESSION,'left'),false);assert.equal(acceptCompareEvent(event,undefined,ORIGIN,SESSION,'left'),false);assert.ok(isCompareEnvelope(data,SESSION,'left'));
});

test('Readiness tolerates arrival order and transfers only the cached leader state and pose',()=>{
 const {relay,sent}=relayFixture(),snapshot={selection:[{id:'left-mesh'}]},pose={model:'left-model'};
 assert.ok(relay.receive('left',envelope('left','ready',1)));assert.ok(relay.receive('left',envelope('left','snapshot',2,{snapshot})));assert.ok(relay.receive('left',envelope('left','pose',3,{pose})));assert.equal(sent.length,0);
 relay.receive('right',envelope('right','ready',1));assert.equal(relay.leader,'left');assert.deepEqual(sent.map(x=>[x.side,x.message.type,x.message.kind]),[['right','apply','state'],['right','apply','pose'],['left','request-snapshot',undefined]]);assert.deepEqual(sent[0].message.snapshot,snapshot);assert.deepEqual(sent[1].message.pose,pose);assert.deepEqual(sent.filter(x=>x.message.type==='apply').map(x=>x.message.revision),[1,2]);assert.ok(sent.every(x=>x.message.instance===relay.instances[x.side]));
});

test('First document binding preserves the default left leader when the right viewer becomes ready first',()=>{
 const {relay,sent}=relayFixture();relay.receive('right',envelope('right','ready',1));assert.equal(relay.receive('right',envelope('right','pose',2,{pose:{model:'right-model',camera:['right-initial']}})),false);assert.equal(relay.leader,'left');
 relay.resetSide('left');assert.equal(relay.leader,'left');relay.receive('left',envelope('left','ready',1));assert.deepEqual(sent.map(x=>[x.side,x.message.type]),[['left','request-snapshot']]);
 const leftPose={model:'left-model',camera:['left-explicit-view']};relay.receive('left',envelope('left','pose',2,{pose:leftPose}));assert.equal(sent.at(-1).side,'right');assert.deepEqual(sent.at(-1).message.pose,leftPose);assert.ok(!sent.some(x=>x.side==='left'&&x.message.kind==='pose'),'An unsolicited right initial pose cannot replace the left explicit view');
 const gesture=relayFixture();gesture.relay.receive('right',envelope('right','ready',1));gesture.relay.receive('right',envelope('right','interaction',2,{phase:'start'}));gesture.relay.resetSide('left');assert.equal(gesture.relay.leader,'right','A real right-hand gesture still takes precedence over the default');
});

test('A model reload after many poses resets only its document epoch and rejects all late old-document traffic',()=>{
 const {relay,sent}=relayFixture();ready(relay);sent.length=0;relay.receive('left',envelope('left','snapshot',2,{snapshot:{selection:['old-left']}}));for(let seq=3;seq<=50;seq++)relay.receive('left',envelope('left','pose',seq,{pose:{model:'left-model',camera:[seq]}}));
 const rightSnapshot={selection:['right-current']},rightPose={model:'right-model',camera:[900]};relay.receive('right',envelope('right','snapshot',2,{snapshot:rightSnapshot}));relay.receive('right',envelope('right','pose',3,{pose:rightPose}));const rightIdentity=structuredClone(relay.identities.right),rightSequence=relay.lastSeq.right,revision=relay.revision,oldInstance=relay.instances.left;relay.resetSide('left');assert.equal(relay.leader,'right');assert.equal(relay.ready.right,true);assert.equal(relay.lastSeq.right,rightSequence);assert.deepEqual(relay.identities.right,rightIdentity);assert.equal(relay.lastSeq.left,0);
 sent.length=0;const nextInstance='new-left-document-5678',identity={model:'female',source:'/models/female.json'};assert.ok(relay.receive('left',envelope('left','ready',1,{instance:nextInstance,identity})));assert.equal(relay.instances.left,nextInstance);assert.deepEqual(relay.identities.left,identity);assert.deepEqual(relay.identities.right,rightIdentity);assert.deepEqual(sent.filter(x=>x.message.type==='apply').map(x=>[x.side,x.message.kind]),[['left','state'],['left','pose']]);assert.deepEqual(sent[0].message.snapshot,rightSnapshot);assert.deepEqual(sent[1].message.pose,rightPose);assert.ok(sent.filter(x=>x.message.type==='apply').every(x=>x.message.instance===nextInstance&&x.message.revision>revision));
 const count=sent.length;for(const data of [envelope('left','pose',9999,{instance:oldInstance,pose:{camera:[-1]}}),envelope('left','ready',9999,{instance:oldInstance}),envelope('left','interaction',9999,{instance:oldInstance,phase:'start'})])assert.equal(relay.receive('left',data),false);assert.equal(sent.length,count);assert.equal(relay.leader,'right');assert.equal(relay.lastSeq.left,1);assert.equal(relay.lastSeq.right,rightSequence);assert.deepEqual(relay.identities.right,rightIdentity);
 assert.ok(relay.receive('left',envelope('left','ready',2,{instance:nextInstance,identity})));assert.equal(sent.length,count,'A duplicate ready response cannot cause an infinite request-snapshot loop');
});

test('Reloading the former leader preserves the other pane last-applied shared view even before its next pose notification',()=>{
 const {relay,sent}=relayFixture();ready(relay);const snapshot={selection:['shared-current']},pose={model:'left-model',camera:[88]};relay.receive('left',envelope('left','snapshot',2,{snapshot}));relay.receive('left',envelope('left','pose',3,{pose}));relay.resetSide('left');sent.length=0;relay.receive('left',envelope('left','ready',1,{instance:'new-left-document',identity:{model:'female',source:'/models/female.json'}}));assert.deepEqual(sent.filter(x=>x.message.type==='apply').map(x=>[x.side,x.message.kind]),[['left','state'],['left','pose']]);assert.deepEqual(sent[0].message.snapshot,snapshot);assert.deepEqual(sent[1].message.pose,pose);assert.deepEqual(relay.identities.right,{model:'right-model',source:'/models/right.json'});
});

test('Follower echoes never lead; an actual gesture explicitly transfers leadership in either direction',()=>{
 const {relay,sent}=relayFixture();ready(relay);sent.length=0;
 const state={selection:['source']};assert.ok(relay.receive('left',envelope('left','snapshot',2,{snapshot:state})));assert.equal(sent.length,1);assert.equal(sent[0].side,'right');
 assert.equal(relay.receive('right',envelope('right','snapshot',2,{snapshot:{selection:['mapped-target']}})),false);assert.equal(relay.receive('right',envelope('right','pose',3,{pose:{camera:[1]}})),false);assert.equal(sent.length,1);assert.equal(relay.leader,'left');
 assert.ok(relay.receive('right',envelope('right','interaction',4,{phase:'start'})));assert.equal(relay.leader,'right');assert.ok(relay.receive('right',envelope('right','snapshot',5,{snapshot:{selection:['user-chosen-target']}})));assert.equal(sent.at(-1).side,'left');assert.equal(sent.at(-1).message.kind,'state');
 relay.receive('right',envelope('right','interaction',6,{phase:'end'}));assert.equal(relay.leader,'right','Gesture end does not give a programmatic echo leadership');relay.receive('left',envelope('left','interaction',3,{phase:'start'}));relay.receive('left',envelope('left','pose',4,{pose:{camera:[2]}}));assert.equal(sent.at(-1).side,'right');assert.equal(sent.at(-1).message.kind,'pose');assert.deepEqual(sent.map(x=>x.message.revision),[1,2,3]);
});

test('Stale, repeated or malformed relay messages cannot advance the source sequence or revision',()=>{
 const {relay,sent}=relayFixture();ready(relay);sent.length=0;
 for(const bad of [envelope('left','snapshot',2,{snapshot:null}),envelope('left','pose',2,{pose:[]}),envelope('left','interaction',2,{phase:'sideways'}),envelope('left','snapshot',0,{snapshot:{}}),envelope('right','snapshot',2,{snapshot:{}})])assert.equal(relay.receive('left',bad),false);
 assert.equal(relay.lastSeq.left,1);assert.equal(relay.revision,0);assert.ok(relay.receive('left',envelope('left','snapshot',2,{snapshot:{labels:true}})));assert.equal(relay.revision,1);
 for(const seq of [2,1])assert.equal(relay.receive('left',envelope('left','pose',seq,{pose:{}})),false);assert.equal(sent.length,1);assert.equal(relay.lastSeq.left,2);
});

test('Parent forwarding leaves viewer model/source identity intact and does not navigate either frame',()=>{
 const {relay,sent}=relayFixture();ready(relay);sent.length=0;const snapshot={selection:[{id:'source-left'}],model:'untrusted-model',source:'/not-the-target.json'},saved=structuredClone(snapshot);relay.receive('left',envelope('left','snapshot',2,{snapshot}));assert.deepEqual(snapshot,saved);assert.equal(sent[0].message.type,'apply');assert.equal(sent[0].message.kind,'state');assert.equal(sent[0].message.model,undefined);assert.equal(sent[0].message.source,undefined);assert.equal(sent[0].message.url,undefined);
});

test('Copying requests both live view URLs and preserves pane order when replies arrive in reverse order',async()=>{
 const sent=[],requests=new ComparisonViewRequests(SESSION,BASE,(side,message)=>sent.push({side,message}),1000),promise=requests.requestPair();assert.equal(requests.pendingCount,2);assert.deepEqual(sent.map(x=>x.side),['left','right']);assert.ok(sent.every(x=>x.message.type==='request-view'&&x.message.session===SESSION));assert.notEqual(sent[0].message.requestId,sent[1].message.requestId);
 const liveLeft=left.replace('left-artery','latest-left-artery'),liveRight=right.replace('right-vein','latest-right-vein'),reply=(index,url)=>envelope(sent[index].side,'view',5,{requestId:sent[index].message.requestId,url});
 assert.ok(requests.receive('right',reply(1,liveRight)));assert.equal(requests.pendingCount,1);assert.ok(requests.receive('left',reply(0,liveLeft)));assert.deepEqual(await promise,[liveLeft,liveRight]);assert.equal(requests.pendingCount,0);assert.equal(requests.receive('left',reply(0,liveLeft)),false,'A duplicate response cannot resolve another copy request');
});

test('A live view reply must match the active session, pane and request and pass URL confinement',async()=>{
 const sent=[],requests=new ComparisonViewRequests(SESSION,BASE,(side,message)=>sent.push({side,message}),1000),promise=requests.request('left'),valid=envelope('left','view',2,{requestId:sent[0].message.requestId,url:left});
 for(const bad of [{...valid,session:'stale-session'},{...valid,side:'right'},{...valid,requestId:'not-requested'},{...valid,url:null},{...valid,seq:0},{...valid,type:'snapshot'}])assert.equal(requests.receive('left',bad),false);assert.equal(requests.pendingCount,1);
 assert.ok(requests.receive('left',{...valid,url:'https://remote.example/human-atlas/?model=other'}));await assert.rejects(promise,/invalid atlas link/);assert.equal(requests.pendingCount,0);
});

test('Copy requests time out or cancel without retaining stale pending work',async()=>{
 const timed=new ComparisonViewRequests(SESSION,BASE,()=>{},10),pending=timed.request('right');await assert.rejects(pending,/did not return its current view/);assert.equal(timed.pendingCount,0);
 const sent=[],cancelled=new ComparisonViewRequests(SESSION,BASE,(side,message)=>sent.push({side,message}),1000),pair=cancelled.requestPair();const rejected=assert.rejects(pair,/views changed/i);cancelled.cancel();await rejected;assert.equal(cancelled.pendingCount,0);assert.equal(cancelled.receive('left',envelope('left','view',1,{requestId:sent[0].message.requestId,url:left})),false);
 const failed=new ComparisonViewRequests(SESSION,BASE,()=>{throw Error('Frame unavailable');},1000);await assert.rejects(failed.request('left'),/Frame unavailable/);assert.equal(failed.pendingCount,0);
});

test('Live-copy replies are accepted for sequencing without becoming pose or selection updates',()=>{
 const {relay,sent}=relayFixture();ready(relay);sent.length=0;assert.ok(relay.receive('left',envelope('left','view',2,{requestId:'copy-1',url:left})));assert.equal(relay.lastSeq.left,2);assert.equal(relay.revision,0);assert.equal(sent.length,0);assert.equal(relay.leader,'left');
});

test('Continuous URL saves read the latest paired selections and each model own live POV for reload',async()=>{
 let live=[left,right],href=BASE;const writes=[],sync=new ComparisonUrlSync(async()=>[...live],pair=>{href=makeComparisonLink(...pair,BASE);writes.push(href);},5);
 sync.changed();sync.changed();await new Promise(resolve=>setTimeout(resolve,15));assert.equal(writes.length,1);
 const male=new URL(left),female=new URL(right);male.searchParams.set('frame','male-own-live-camera');female.searchParams.set('frame','female-own-live-camera');male.searchParams.set('select','new-left-selection');female.searchParams.set('select','mapped-right-selection');male.hash='atlas-view=male-latest';female.hash='atlas-view=female-latest';live=[male.href,female.href];sync.changed();await new Promise(resolve=>setTimeout(resolve,15));
 const restored=new URL(href),first=new URL(restored.searchParams.get('first')),second=new URL(restored.searchParams.get('second'));
 assert.equal(writes.length,2);assert.equal(first.href,male.href);assert.equal(second.href,female.href);assert.equal(first.searchParams.get('model'),'male-detail');assert.equal(second.searchParams.get('model'),'female');assert.notEqual(first.searchParams.get('frame'),second.searchParams.get('frame'));assert.equal(first.searchParams.get('atlasSource'),'/models/male.json');assert.equal(second.searchParams.get('atlasSource'),'/models/female.json');assert.equal(first.hash,'#atlas-view=male-latest');assert.equal(second.hash,'#atlas-view=female-latest');sync.cancel();
});

test('A change during a live capture schedules the latest trailing capture and gesture end flushes immediately',async()=>{
 const reads=[],writes=[],sync=new ComparisonUrlSync(()=>new Promise(resolve=>reads.push(resolve)),pair=>writes.push(pair),5);
 sync.flush();assert.equal(reads.length,1);sync.changed();sync.changed();reads[0]([left,right]);await new Promise(resolve=>setTimeout(resolve,15));assert.equal(reads.length,2);
 const last=[left.replace('saved-left','final-left'),right.replace('saved-right','final-right')];reads[1](last);await new Promise(resolve=>setTimeout(resolve,0));assert.deepEqual(writes.at(-1),last);
 sync.changed();sync.flush();assert.equal(reads.length,3);reads[2](last);await new Promise(resolve=>setTimeout(resolve,0));assert.equal(writes.length,3);sync.cancel();
});

test('Reloading or replacing a pane invalidates pending URL writes from its old document',async()=>{
 const reads=[],writes=[],sync=new ComparisonUrlSync(()=>new Promise(resolve=>reads.push(resolve)),pair=>writes.push(pair),5);sync.flush();sync.cancel();reads[0]([left,right]);await new Promise(resolve=>setTimeout(resolve,0));assert.equal(writes.length,0);
 sync.flush();const current=[left.replace('male-detail','female'),right];reads[1](current);await new Promise(resolve=>setTimeout(resolve,0));assert.deepEqual(writes,[current]);sync.cancel();
});

test('Landscape uses two equal columns and portrait uses two equal rows with bounded iframe panes',async()=>{
 const root=postcss.parse(await readFile('app/comparison-page.css','utf8'));
 const declarations=(selector,media)=>{const out={};root.walkRules(selector,rule=>{const parent=rule.parent;const actualMedia=parent.type==='atrule'&&parent.name==='media'?parent.params:null;if(actualMedia===media)rule.walkDecls(d=>out[d.prop]=d.value);});return out;};
 const landscape=declarations('.comparison',null),portrait=declarations('.comparison','(orientation:portrait)');assert.equal(landscape.display,'grid');assert.equal(landscape['grid-template-columns'],'minmax(0,1fr) minmax(0,1fr)');assert.equal(landscape['grid-template-rows'],'minmax(0,1fr)');assert.equal(portrait['grid-template-columns'],'minmax(0,1fr)');assert.equal(portrait['grid-template-rows'],'minmax(0,1fr) minmax(0,1fr)');assert.equal(landscape.height,'100dvh');const iframe=declarations('.comparison iframe',null);assert.equal(iframe.height,'100%');assert.equal(iframe.width,'100%');assert.equal(iframe['min-height'],'0');assert.equal(iframe.border,'0');
});

test('The comparison wrapper contains exactly two accessible viewports and no wrapper chrome',async()=>{
 const html=await readFile('app/comparison-page.html','utf8');assert.equal((html.match(/<iframe\b/g)??[]).length,2);for(const id of ['first-view','second-view'])assert.match(html,new RegExp('<iframe[^>]*id="'+id+'"[^>]*title="[^"]+"'));const body=html.split('<body>')[1].split('</body>')[0];assert.doesNotMatch(body,/<(?:header|footer|nav|form|input|button|h[1-6]|p)\b/);assert.doesNotMatch(body,/First view|Second view|comparison-status|comparison-sources/);assert.match(html,/type="module"[^>]*src="\/compare\.ts"/);
});

test('Comparison controls keep Reset beside the right model picker and hide duplicate camera navigation',async()=>{
 const page=await readFile('app/page.tsx','utf8');assert.match(page,/\{\(!comparison\|\|comparison.side==='right'\)&&<button[^>]*aria-label="Assemble and reset"/);assert.match(page,/if\(comparison&&control==='model'\)return true/);assert.match(page,/if\(comparison\)window\.location\.assign\(comparisonModelUrl\(currentViewUrl\(\),value!,comparison\)\)/);assert.match(page,/\{show\('camera'\)&&<nav className="view-controls/);const css=postcss.parse(await readFile('app/globals.css','utf8'));let cameraRule=false;css.walkRules('.studio[data-compare-side] .view-controls',rule=>rule.walkDecls('display',d=>{if(d.value==='none'&&d.important)cameraRule=true;}));assert.ok(cameraRule,'Any stale comparison camera markup is hidden without changing standalone CSS');
 const declarations=selector=>{const out={};css.walkRules(rule=>{if(rule.selectors.includes(selector))rule.walkDecls(d=>out[d.prop]=d.value);});return out;};const chooser=declarations('.studio[data-compare-side=right] .identity');assert.equal(chooser.top,'auto');assert.equal(chooser.bottom,'5px');assert.equal(chooser.right,'5px');assert.equal(chooser.width,'var(--comparison-detail-width)');const details=declarations('body .studio[data-compare-side=right] .detail-sheet');assert.equal(details.width,chooser.width);assert.equal(details.right,chooser.right);assert.equal(details.bottom,'54px');const leftChooser=declarations('.studio[data-compare-side=left] .identity'),systems=declarations('.studio[data-compare-side=left] .layers-panel.mobile-open');assert.equal(leftChooser.width,systems.width);assert.equal(systems.top,'54px');const actions=declarations('.studio[data-compare-side=right] .top-actions');assert.equal(actions.top,'5px');assert.equal(actions.right,'5px');let reserved=false;css.walkRules(rule=>{if(rule.selector.startsWith('body .studio[data-compare-side=right] :is('))rule.walkDecls('bottom',d=>{if(d.value==='54px'&&d.important)reserved=true;});});assert.ok(reserved,'Right panels leave the original five-pixel gap above the 44px model picker and Reset capsule');
});

test('The build generator emits the tracked comparison template and module shim without rewriting unchanged entries',async()=>{
 const root=await mkdtemp(join(tmpdir(),'atlas-compare-test-'));
 try{
  await mkdir(join(root,'app'));const template=await readFile('app/comparison-page.html','utf8');await writeFile(join(root,'app/comparison-page.html'),template);await writeFile(join(root,'app/comparison-page.ts'),"import './comparison-page.css';\n");
  const generated=await prepareComparisonPage(root);assert.equal(generated.html,join(root,'web/compare.html'));assert.equal(generated.entry,join(root,'web/compare.ts'));assert.equal(generated.sourceHtml,join(root,'app/comparison-page.html'));assert.equal(generated.sourceEntry,join(root,'app/comparison-page.ts'));assert.equal(await readFile(generated.html,'utf8'),template);assert.equal(await readFile(generated.entry,'utf8'),"import '../app/comparison-page';\n");
  const previous=new Date('2001-01-01T00:00:00.000Z');for(const path of [generated.html,generated.entry])await utimes(path,previous,previous);assert.deepEqual(await prepareComparisonPage(root),generated);for(const path of [generated.html,generated.entry])assert.equal((await stat(path)).mtimeMs,previous.getTime(),'Unchanged generated entries retain their timestamps and avoid dev-server reloads');assert.equal(await readFile(generated.sourceHtml,'utf8'),template);
 }finally{await rm(root,{recursive:true,force:true});}
});

test('Local comparison defaults retain matching complete saved selections while production defaults remain public',()=>{
 const first=new URL(LOCAL_DEFAULT_VIEWER_URLS.first,ORIGIN),second=new URL(LOCAL_DEFAULT_VIEWER_URLS.second,ORIGIN);for(const key of ['model','p','s','h','layers','frame','context','tree'])assert.equal(first.searchParams.get(key),second.searchParams.get(key));assert.notEqual(first.searchParams.get('atlasSource'),second.searchParams.get('atlasSource'));assert.ok(first.searchParams.get('atlasSource').includes('v2938'));assert.ok(second.searchParams.get('atlasSource').includes('v2941'));assert.deepEqual(DEFAULT_VIEWER_URLS,{first:'./?model=male-detail',second:'./?model=female'});
});
