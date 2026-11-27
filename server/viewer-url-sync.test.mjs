import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
import postcss from 'postcss';

const built=await build({stdin:{contents:"export * from './app/viewer-location';export {viewParameters,compactViewUrl} from './app/view-url';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {createViewUrlSync,preserveViewerRouting,navigateViewerUrl,viewParameters,compactViewUrl}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
function fixture(){let next=0,url='https://atlas.example/?model=male&frame=first',reads=0;const timers=new Map(),writes=[];const sync=createViewUrlSync({getUrl:()=>{reads++;return url;},replace:value=>writes.push(value),setTimeout:(fn,delay)=>{const id=++next;timers.set(id,{fn,delay});return id;},clearTimeout:id=>timers.delete(id)});return{sync,timers,writes,setUrl:value=>url=value,reads:()=>reads,tick(){const callbacks=[...timers.values()];timers.clear();callbacks.forEach(item=>item.fn());}};}

test('Live camera updates compress and replace only the latest pose at a bounded cadence',()=>{
 const f=fixture();for(let n=0;n<100;n++){f.setUrl('https://atlas.example/?frame='+n);f.sync.schedule();}assert.equal(f.timers.size,1);assert.equal(f.reads(),0,'Camera events do not repeatedly serialize or compress the atlas');assert.equal([...f.timers.values()][0].delay,300);f.tick();assert.deepEqual(f.writes,['https://atlas.example/?frame=99']);assert.equal(f.reads(),1);
 f.setUrl('https://atlas.example/?frame=still-dragging');f.sync.schedule();f.tick();assert.equal(f.writes.at(-1),'https://atlas.example/?frame=still-dragging','An ongoing drag is periodically reflected instead of waiting forever for a debounce');f.sync.schedule();f.tick();assert.equal(f.writes.length,2,'The same state does not create redundant history writes');f.sync.dispose();
});

test('Gesture-end and page-hide flushes publish the latest view without leaving a trailing old timer',()=>{
 const f=fixture();f.sync.schedule();f.setUrl('https://atlas.example/?frame=final');f.sync.flush();assert.deepEqual(f.writes,['https://atlas.example/?frame=final']);assert.equal(f.timers.size,0);f.sync.flush();assert.equal(f.writes.length,1);f.setUrl(null);f.sync.schedule();f.tick();assert.equal(f.writes.length,1,'A loading or unmounted model has no view URL to install');f.sync.dispose();
});

test('Disposing URL observation cancels its pending work and never writes an old model after navigation',()=>{
 const f=fixture();f.sync.schedule();f.sync.dispose();assert.equal(f.timers.size,0);f.sync.flush();f.sync.schedule();f.tick();assert.equal(f.reads(),0);assert.equal(f.writes.length,0);
});

test('Serialized live URLs retain exact source, review, iframe and comparison routing without replaying focus intent',()=>{
 const current=new URLSearchParams({model:'local-ta98',atlasSource:'/local-models/current.json',blood:'biodigital-review',cardiac:'review',bowel:'layers',hepatic:'completion-final',biliary:'folds',venous:'confluence',caps:'shared-v28',cricoid:'partition-v1',cricotracheal:'first-ring-v122',delivery:'1',review:'alternatives',cardiacRefStage:'current',embed:'1',ui:'model,details',compareSession:'session-1234',compareSide:'right',connect:'room',connectJoin:'1',relay:'local',focus:'1'}),url=new URL('https://atlas.example/?model=local-ta98&frame=actual-live-frame&select=current-selection');preserveViewerRouting(url,current);
 for(const [key,value]of current)if(!['model','focus'].includes(key))assert.equal(url.searchParams.get(key),value,key);assert.equal(url.searchParams.get('focus'),null,'The stored live frame replaces the one-time centering request');assert.equal(url.searchParams.get('frame'),'actual-live-frame');assert.equal(url.searchParams.get('select'),'current-selection');const restored=viewParameters(compactViewUrl(url.href));assert.equal(restored.get('atlasSource'),current.get('atlasSource'));assert.equal(restored.get('compareSide'),'right');assert.equal(restored.get('frame'),'actual-live-frame');
});

test('URL replacement changes no application state, browser events or history entries',()=>{
 const oldHistory=globalThis.history,oldDispatch=globalThis.dispatchEvent,calls=[];globalThis.history={replaceState:(...args)=>calls.push(['replace',...args]),pushState:(...args)=>calls.push(['push',...args])};globalThis.dispatchEvent=()=>{throw Error('A URL refresh must not emit navigation events');};
 try{navigateViewerUrl('https://atlas.example/?frame=latest',true);assert.deepEqual(calls,[['replace',null,'','https://atlas.example/?frame=latest']]);}finally{if(oldHistory===undefined)delete globalThis.history;else globalThis.history=oldHistory;if(oldDispatch===undefined)delete globalThis.dispatchEvent;else globalThis.dispatchEvent=oldDispatch;}
});

test('Page wiring observes live cameras and settling without re-registering WebMCP tools or replaying navigation',async()=>{
 const page=await readFile('app/page.tsx','utf8');const registration=page.slice(page.indexOf('return registerAtlasTools('),page.indexOf('const choosePart='));assert.match(registration,/\},\[atlas\]\)/);assert.match(registration,/agentActionsRef\.current/);assert.doesNotMatch(registration,/\},\[atlas,/,'URL/state values do not extend the tool-registration lifecycle');const start=page.indexOf('const sync=createViewUrlSync'),observer=page.slice(start,page.indexOf('const openSavedView',start));assert.match(observer,/navigateViewerUrl\(url,true\)/);assert.doesNotMatch(observer,/dispatchEvent|setState\(|loadState\(|setAtlas\(|setModel\(/);for(const event of ['pointerup','wheel','pagehide','visibilitychange'])assert.ok(observer.includes("'"+event+"'"));assert.match(page,/onCamera=\{camera=>\{cameraRef\.current=camera;scheduleViewUrlRef\.current\(\);\}\}/);
});

test('The left comparison never renders the standalone Views button or preset camera navigation',async()=>{
 const page=await readFile('app/page.tsx','utf8');assert.match(page,/\{!embedded&&comparison\?\.side!=='left'&&<Button[^>]*aria-label="Save or open views"/);assert.match(page,/\{show\('camera'\)&&<nav className="view-controls/);
});

test('Only the right comparison chooser moves below panels with a matching width and an attached Reset',async()=>{
 const page=await readFile('app/page.tsx','utf8');assert.match(page,/\{\(!comparison\|\|comparison\.side==='right'\)&&<button[^>]*aria-label="Assemble and reset"/);const css=postcss.parse(await readFile('app/globals.css','utf8'));const rules=selector=>{const out={};css.walkRules(selector,r=>r.walkDecls(d=>out[d.prop]=d.value));return out;},chooser=rules('.studio[data-compare-side=right] .identity');assert.equal(chooser.bottom,'5px');assert.equal(chooser.right,'5px');assert.equal(chooser.top,'auto');assert.equal(chooser.width,'var(--comparison-detail-width)');assert.equal(rules('body .studio[data-compare-side=right] .detail-sheet').width,'var(--comparison-detail-width)');let reserve=false;css.walkRules(r=>{if(r.selector.startsWith('body .studio[data-compare-side=right] :is(')&&['.detail-sheet','.search-panel','.section-panel','.advanced-panel','.study-panel','.visit-panel','.download-panel'].every(part=>r.selector.includes(part)))r.walkDecls('bottom',d=>{if(d.value==='54px'&&d.important)reserve=true;});});assert.ok(reserve,'Every panel sharing the right-bottom position reserves the chooser strip');
});
