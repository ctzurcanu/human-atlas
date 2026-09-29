import {test,before,after} from 'node:test';
import assert from 'node:assert/strict';
import {randomBytes} from 'node:crypto';
import {build} from 'esbuild';
import {Miniflare} from 'miniflare';
import {fileURLToPath} from 'node:url';

const code=()=>randomBytes(24).toString('base64url');
const viewer='https://ctzurcanu.github.io',origin='https://relay.example';
let mf;
before(async()=>{
 const bundle=await build({entryPoints:[fileURLToPath(new URL('connect-worker.mjs',import.meta.url))],bundle:true,format:'esm',platform:'browser',external:['cloudflare:workers'],loader:{'.html':'text'},write:false});
 // Only the test bundle exports this RPC method. Reconstructing the object
 // exercises exactly the constructor used when Cloudflare wakes a sleeping room.
 const script=bundle.outputFiles[0].text+`\nexport class TestAtlasRoom extends AtlasRoom {
  reconstructForTest(){Object.assign(this,new AtlasRoom(this.ctx,this.env));}
  expireForTest(){this.meta.hostMissingUntil=Date.now()-1;this.save();return this.alarm();}
 }`;
 mf=new Miniflare({modules:true,script,compatibilityDate:'2026-04-07',bindings:{ALLOWED_VIEWER_ORIGINS:viewer},durableObjects:{ROOMS:{className:'TestAtlasRoom',useSQLite:true}}});
});
after(async()=>{await mf?.dispose();});
async function fixture(t){
 const sockets=[];
 t.after(()=>{for(const socket of sockets)try{socket.close(1000,'Test complete');}catch{}});
 async function client({roomId=code(),role='host',clientToken=code(),ip='203.0.113.7',name=''}={}){
  const response=await mf.dispatchFetch(`${origin}/atlas-connect?room=${roomId}`,{headers:{Upgrade:'websocket',Origin:viewer,'CF-Connecting-IP':ip}});
  assert.equal(response.status,101);const socket=response.webSocket,queue=[],waiters=[];sockets.push(socket);
  socket.addEventListener('message',event=>{
   const message=JSON.parse(event.data),index=waiters.findIndex(waiter=>waiter.type===message.type&&waiter.predicate(message));
   if(index<0)queue.push(message);else{const waiter=waiters.splice(index,1)[0];clearTimeout(waiter.timer);waiter.resolve(message);}
  });
  const next=(type,predicate=()=>true)=>{
   const index=queue.findIndex(message=>message.type===type&&predicate(message));if(index>=0)return Promise.resolve(queue.splice(index,1)[0]);
   return new Promise((resolve,reject)=>{const waiter={type,predicate,resolve,timer:setTimeout(()=>reject(new Error(`Timed out waiting for ${type}`)),3000)};waiters.push(waiter);});
  };
  socket.accept();await next('hello');const send=message=>socket.send(JSON.stringify(message));send({type:'join',role,roomId,clientToken,name});
  return {socket,next,send,roomId,clientToken,queue};
 }
 return {client};
}
const pose={model:'male-detail',camera:[1,1,3,0,.8,0,0,0],up:[0,1,0],rotation:[0,0,0,1]};
test('public MCP initializes, discovers tools, searches and returns the GitHub Pages UI',async()=>{
 const rpc=async(method,params={})=>{
  const response=await mf.dispatchFetch(`${origin}/mcp`,{method:'POST',headers:{'Content-Type':'application/json',Accept:'application/json, text/event-stream'},body:JSON.stringify({jsonrpc:'2.0',id:1,method,params})});
  assert.equal(response.status,200);return (await response.json()).result;
 };
 const init=await rpc('initialize',{protocolVersion:'2025-06-18',capabilities:{},clientInfo:{name:'test',version:'1'}});
 assert.equal(init.serverInfo.name,'human-atlas');
 const tools=(await rpc('tools/list')).tools;
 assert.equal(tools.length,3);
 assert.ok(!JSON.stringify(tools.find(t=>t.name==='show_anatomy').inputSchema).includes('local-male'));
 const search=await rpc('tools/call',{name:'search_anatomy',arguments:{query:'Stomach'}});
 assert.ok(search.structuredContent.matches.some(m=>m.id==='ZA:Stomach'));
 const view=await rpc('tools/call',{name:'show_anatomy',arguments:{structure:'Stomach'}});
 assert.equal(new URL(view.structuredContent.url).origin,viewer);
 const options=await rpc('tools/call',{name:'get_anatomy_options',arguments:{}});
 assert.ok(options.structuredContent.models.every(id=>!id.startsWith('local-')));
 const resources=(await rpc('resources/list')).resources;
 const ui=await rpc('resources/read',{uri:resources[0].uri});
 assert.ok(ui.contents[0].text.includes('ui/initialize'));
 assert.deepEqual(ui.contents[0]._meta.ui.csp.frameDomains,[viewer]);
 const forbidden=await mf.dispatchFetch(`${origin}/mcp`,{method:'POST',headers:{Origin:'https://untrusted.example'}});
 assert.equal(forbidden.status,403);
});
const snapshot={version:1,model:'male-detail',hierarchy:'systems',guestSources:[],state:{focus:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:true,explode:0,visible:['muscular'],selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0},choice:null,ui:{panel:'advanced',frontPanel:'advanced',details:false,layersVisible:true,mobileLayersOpen:false,addSelection:false,query:''},covering:[],pose,advanced:{document:{schemaVersion:1,type:'connect',views:[]},slideIndex:null}};
test('health and identity return the relay origin; unrelated browser origins are rejected',async()=>{
 const health=await mf.dispatchFetch(`${origin}/health`);assert.equal((await health.json()).status,'ok');
 const identity=await mf.dispatchFetch(`${origin}/atlas-connect/identity`,{headers:{Origin:viewer,'CF-Connecting-IP':'198.51.100.4'}});
 assert.equal(identity.headers.get('Access-Control-Allow-Origin'),viewer);assert.deepEqual(await identity.json(),{ip:'198.51.100.4',origin});
 for(const path of ['/atlas-connect/identity',`/atlas-connect?room=${code()}`]){
  const rejected=await mf.dispatchFetch(origin+path,{headers:{Upgrade:'websocket',Origin:'https://unrelated.example'}});assert.equal(rejected.status,403);
 }
 assert.equal((await mf.dispatchFetch(`${origin}/atlas-connect?room=bad`,{headers:{Upgrade:'websocket',Origin:viewer}})).status,400);
});
test('waiting guests, Play cohort, live state and poses, and Stop match the existing viewer',async t=>{
 const {client}=await fixture(t),host=await client();const joined=await host.next('joined');assert.equal(joined.roomId,host.roomId);
 const guest=await client({role:'guest',roomId:host.roomId});assert.equal((await guest.next('joined')).playing,false);
 await host.next('roster',message=>message.guests.length===1);
 host.send({type:'snapshot',snapshot});host.send({type:'play'});await host.next('playing');await guest.next('playing');assert.equal(guest.queue.some(message=>message.type==='snapshot'),false);
 host.send({type:'snapshot',snapshot});assert.deepEqual((await guest.next('snapshot')).snapshot,snapshot);
 const moved={...pose,camera:[2,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:moved});assert.deepEqual((await guest.next('pose')).pose,moved);
 const late=await client({role:'guest',roomId:host.roomId});assert.equal((await late.next('joined')).playing,false);
 await host.next('roster',message=>message.guests.length===2);
 host.send({type:'pose',pose});await guest.next('pose');assert.equal(late.queue.some(message=>message.type==='pose'||message.type==='snapshot'),false);
 host.send({type:'stop'});await guest.next('ended');await late.next('ended');await host.next('stopped');await host.next('roster',message=>message.guests.length===0);
});
test('guest controls cannot broadcast or kick; host rename and kick persist across reconnection',async t=>{
 const {client}=await fixture(t),host=await client();await host.next('joined');
 const guest=await client({role:'guest',roomId:host.roomId});const joined=await guest.next('joined');await host.next('roster',message=>message.guests.length===1);
 host.send({type:'play'});await host.next('playing');await guest.next('playing');
 guest.send({type:'snapshot',snapshot});guest.send({type:'kick',id:joined.host.id});guest.send({type:'rename',id:joined.self.id,name:'Student'});assert.equal((await guest.next('self')).self.name,'Student');assert.equal(host.queue.some(message=>message.type==='snapshot'),false);
 host.send({type:'rename',id:joined.self.id,name:'Learner'});assert.equal((await guest.next('self')).self.name,'Learner');
 host.send({type:'kick',id:joined.self.id});await guest.next('ended');
 const kicked=await client({role:'guest',roomId:host.roomId,clientToken:guest.clientToken});await kicked.next('ended');
 const intruder=await client({roomId:host.roomId});assert.match((await intruder.next('error')).message,/another host/);
 const missing=await client({role:'guest'});assert.match((await missing.next('error')).message,/not running/);
});
test('the host admits a waiting guest with current state and camera, preserving live guests and hibernation',async t=>{
 const {client}=await fixture(t),host=await client();await host.next('joined');
 const guest=await client({role:'guest',roomId:host.roomId});await guest.next('joined');
 host.send({type:'play'});await host.next('playing');await guest.next('playing');
 host.send({type:'snapshot',snapshot});await guest.next('snapshot');
 const moved={...pose,camera:[4,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:moved});await guest.next('pose');
 const late=await client({role:'guest',roomId:host.roomId}),joined=await late.next('joined');assert.equal(joined.playing,false);
 guest.send({type:'admit',id:joined.self.id});guest.send({type:'rename',id:(await guest.next('roster',message=>message.guests.length===2)).guests.find(peer=>peer.inSession).id,name:'Original guest'});await guest.next('self');
 assert.equal(late.queue.some(message=>message.type==='playing'),false);
 host.send({type:'admit',id:joined.self.id});await late.next('playing');
 assert.deepEqual((await late.next('snapshot')).snapshot,snapshot);assert.deepEqual((await late.next('pose')).pose,moved);
 await host.next('snapshot-request');await host.next('roster',message=>message.guests.length===2&&message.guests.every(peer=>peer.inSession));
 const next={...moved,camera:[5,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:next});assert.deepEqual((await late.next('pose')).pose,next);assert.deepEqual((await guest.next('pose')).pose,next);
 const rooms=await mf.getDurableObjectNamespace('ROOMS');await rooms.getByName(host.roomId).reconstructForTest();
 late.socket.close(1000,'Reconnect admitted guest');await host.next('roster',message=>message.guests.some(peer=>peer.id===joined.self.id&&!peer.online));
 const resumed=await client({role:'guest',roomId:host.roomId,clientToken:late.clientToken});assert.equal((await resumed.next('joined')).playing,true);assert.deepEqual((await resumed.next('pose')).pose,next);
 host.send({type:'stop'});await host.next('stopped');await guest.next('ended');await resumed.next('ended');
});
test('host and captured guests resume the same room with the latest camera',async t=>{
 const {client}=await fixture(t),host=await client();await host.next('joined');const guest=await client({role:'guest',roomId:host.roomId});await guest.next('joined');
 host.send({type:'play'});await host.next('playing');await guest.next('playing');host.send({type:'snapshot',snapshot});await guest.next('snapshot');
 const moved={...pose,camera:[3,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:moved});await guest.next('pose');
 host.socket.close(1000,'Reconnect');await guest.next('roster',message=>!message.host.online);
 const resumed=await client({roomId:host.roomId,clientToken:host.clientToken});assert.equal((await resumed.next('joined')).playing,true);
 guest.socket.close(1000,'Reconnect');await resumed.next('roster',message=>message.guests.some(guest=>!guest.online));
 const rejoined=await client({role:'guest',roomId:host.roomId,clientToken:guest.clientToken});assert.equal((await rejoined.next('joined')).playing,true);assert.deepEqual((await rejoined.next('snapshot')).snapshot,snapshot);assert.deepEqual((await rejoined.next('pose')).pose,moved);
 rejoined.send({type:'stop'});await rejoined.next('ended');await resumed.next('roster',message=>message.guests.length===0);
 resumed.send({type:'stop'});await resumed.next('stopped');
});
test('one host command includes 100 waiting guests with a large view and streams motion to every guest',async t=>{
 const {client}=await fixture(t),host=await client();const hostJoined=await host.next('joined');
 host.send({type:'play'});await host.next('playing');
 const large={...snapshot,state:{...snapshot.state,hidden:Array.from({length:3756},(_,i)=>`ZA:Hidden anatomical structure ${i}`)}};
 host.send({type:'snapshot',snapshot:large});host.send({type:'rename',id:hostJoined.self.id,name:'Bulk host'});await host.next('self');
 const guests=await Promise.all(Array.from({length:100},()=>client({role:'guest',roomId:host.roomId})));
 const joined=await Promise.all(guests.map(guest=>guest.next('joined')));assert.ok(joined.every(message=>!message.playing));
 guests[0].send({type:'admit-all'});guests[0].send({type:'rename',id:joined[0].self.id,name:'Waiting guest'});await guests[0].next('self');assert.ok(guests.every(guest=>!guest.queue.some(message=>message.type==='playing')));
 host.send({type:'admit-all'});
 await Promise.all(guests.map(async guest=>{await guest.next('playing');assert.deepEqual((await guest.next('snapshot')).snapshot,large);assert.deepEqual((await guest.next('pose')).pose,pose);}));
 await host.next('snapshot-request');await host.next('roster',message=>message.guests.length===100&&message.guests.every(peer=>peer.inSession));
 const moved={...pose,camera:[5,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:moved});await Promise.all(guests.map(async guest=>assert.deepEqual((await guest.next('pose')).pose,moved)));
 host.send({type:'admit-all'});host.send({type:'rename',id:hostJoined.self.id,name:'Still bulk host'});await host.next('self');assert.equal(host.queue.some(message=>message.type==='snapshot-request'),false);
 host.send({type:'stop'});await host.next('stopped');await Promise.all(guests.map(guest=>guest.next('ended')));
});
test('hibernation reconstruction preserves the cohort, names, state, camera and host-only control',async t=>{
 const {client}=await fixture(t),host=await client();await host.next('joined');const guest=await client({role:'guest',roomId:host.roomId});const joined=await guest.next('joined');
 host.send({type:'play'});await host.next('playing');await guest.next('playing');host.send({type:'snapshot',snapshot});await guest.next('snapshot');
 const moved={...pose,camera:[4,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:moved});await guest.next('pose');
 host.send({type:'rename',id:joined.self.id,name:'Hibernating guest'});await guest.next('self');
 const rooms=await mf.getDurableObjectNamespace('ROOMS'),room=rooms.getByName(host.roomId);await room.reconstructForTest();
 guest.socket.close(1000,'Reconnect after hibernation');await host.next('roster',message=>message.guests.some(guest=>!guest.online));
 const rejoined=await client({role:'guest',roomId:host.roomId,clientToken:guest.clientToken});const resumed=await rejoined.next('joined');assert.equal(resumed.playing,true);assert.equal(resumed.self.name,'Hibernating guest');
 assert.deepEqual((await rejoined.next('snapshot')).snapshot,snapshot);assert.deepEqual((await rejoined.next('pose')).pose,moved);
 const intruder=await client({roomId:host.roomId});assert.match((await intruder.next('error')).message,/another host/);
 host.socket.close(1000,'Host gone');await rejoined.next('roster',message=>!message.host.online);await room.expireForTest();assert.match((await rejoined.next('ended')).reason,/host disconnected/);
 const expired=await client({role:'guest',roomId:host.roomId});assert.match((await expired.next('error')).message,/not running/);
});
