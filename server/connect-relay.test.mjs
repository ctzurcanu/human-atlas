import {test,before,after} from 'node:test';
import assert from 'node:assert/strict';
import {createServer} from 'node:https';
import {mkdtempSync,readFileSync,rmSync} from 'node:fs';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import {execFileSync} from 'node:child_process';
import {randomBytes} from 'node:crypto';
import {once} from 'node:events';
import WebSocket from 'ws';
import {attachConnectRelay} from './connect-relay.mjs';
import {validSnapshot,validPose} from '../shared/connect-protocol.mjs';

let directory,credentials;
before(()=>{
 directory=mkdtempSync(join(tmpdir(),'atlas-connect-test-'));
 // This certificate exists only in the test's temporary directory. Production
 // browsers use their normal HTTPS origin and normal certificate validation.
 execFileSync('openssl',['req','-x509','-newkey','rsa:2048','-nodes','-keyout',join(directory,'key.pem'),'-out',join(directory,'cert.pem'),'-days','1','-subj','/CN=localhost'],{stdio:'ignore'});
 credentials={key:readFileSync(join(directory,'key.pem')),cert:readFileSync(join(directory,'cert.pem'))};
});
after(()=>rmSync(directory,{recursive:true,force:true}));
const token=()=>randomBytes(24).toString('base64url');
const pose={model:'male-detail',camera:[1,1,3,0,.8,0,0,0],up:[0,1,0],rotation:[0,0,0,1]};
const snapshot=()=>({version:1,model:'male-detail',hierarchy:'systems',guestSources:[],state:{focus:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:true,explode:0,visible:['muscular'],selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0},choice:null,ui:{panel:'advanced',frontPanel:'advanced',details:false,layersVisible:true,mobileLayersOpen:false,addSelection:false,query:''},covering:[],pose,advanced:{document:{schemaVersion:1,type:'connect',views:[]},slideIndex:null}});
async function fixture(t,graceMs=500){
 const server=createServer(credentials),relay=attachConnectRelay(server,{graceMs}),clients=[];
 server.listen(0,'127.0.0.1');await once(server,'listening');
 const origin=`https://127.0.0.1:${server.address().port}`,url=origin.replace('https:','wss:')+'/atlas-connect';
 t.after(async()=>{for(const socket of clients)socket.terminate();relay.dispose();await new Promise(resolve=>server.close(resolve));});
 async function client({role='host',roomId='',clientToken=token(),ip='203.0.113.7',name='',viewerOrigin=origin}={}){
  const socket=new WebSocket(url,{origin:viewerOrigin,rejectUnauthorized:false,headers:{'x-forwarded-for':ip}}),queue=[],waiters=[];clients.push(socket);
  socket.on('error',()=>{});
  socket.on('message',data=>{const message=JSON.parse(data.toString()),index=waiters.findIndex(waiter=>waiter.type===message.type&&waiter.predicate(message));if(index>=0){const waiter=waiters.splice(index,1)[0];clearTimeout(waiter.timer);waiter.resolve(message);}else queue.push(message);});
  const next=(type,predicate=()=>true)=>{
   const index=queue.findIndex(message=>message.type===type&&predicate(message));if(index>=0)return Promise.resolve(queue.splice(index,1)[0]);
   return new Promise((resolve,reject)=>{const waiter={type,predicate,resolve,timer:setTimeout(()=>reject(new Error(`Timed out waiting for ${type}`)),2000)};waiters.push(waiter);});
  };
  await once(socket,'open');assert.equal((await next('hello')).ip,ip);
  const send=message=>socket.send(JSON.stringify(message));send({type:'join',role,roomId,clientToken,name});
  const joined=await next('joined');return {socket,next,send,joined,clientToken,queue};
 }
 return {client,relay,url,origin,clients};
}
async function play(host,...guests){host.send({type:'play'});await host.next('playing');await Promise.all(guests.map(guest=>guest.next('playing')));}
test('GitHub Pages hosts and guests exchange live state through a separate WSS relay',async t=>{
 const {client}=await fixture(t),viewerOrigin='https://ctzurcanu.github.io';
 const host=await client({viewerOrigin}),guest=await client({viewerOrigin,role:'guest',roomId:host.joined.roomId});
 await play(host,guest);const state=snapshot();host.send({type:'snapshot',snapshot:state});assert.deepEqual((await guest.next('snapshot')).snapshot,state);
 host.send({type:'pose',pose});assert.deepEqual((await guest.next('pose')).pose,pose);
 host.send({type:'stop'});await guest.next('ended');await host.next('stopped');
});
test('WSS defaults to IP names and only the host broadcasts complete viewer state and poses',async t=>{
 const {client,relay}=await fixture(t),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId,ip:'198.51.100.4'});
 assert.equal(host.joined.self.name,'203.0.113.7');assert.equal(guest.joined.self.name,'198.51.100.4');
 await play(host,guest);
 const state=snapshot();state.choice={id:'digestive',name:'Digestive system',elements:['ZA:Stomach'],group:true,children:[{id:'ZA:Stomach',name:'Stomach',elements:['ZA:Stomach'],group:false}]};state.ui.expanded=[{key:'digestive',expanded:true},{key:'details:After 1',expanded:true}];state.relationshipExpansion={anchorId:'ZA:Stomach',anchorIds:['ZA:Stomach'],selectedKey:'ZA:Stomach|ZA:Duodenum',frontier:['ZA:Duodenum'],depth:1,exhausted:false};assert.ok(validSnapshot(state));host.send({type:'snapshot',snapshot:state});assert.deepEqual((await guest.next('snapshot')).snapshot,state);
 const camera={...pose,camera:[.1,.2,.3,0,0,0,.1,-.2],rotation:[0,.2,0,.98]};host.send({type:'pose',pose:camera});assert.deepEqual((await guest.next('pose')).pose,camera);
 guest.send({type:'snapshot',snapshot:{...state,state:{...state.state,explode:1}}});guest.send({type:'rename',id:guest.joined.self.id,name:'Student'});await guest.next('self');
 assert.equal(relay.rooms.get(host.joined.roomId).snapshot.state.explode,0);
 assert.equal(host.queue.filter(message=>message.type==='snapshot').length,0);
});
test('sessions are isolated, names update, and guest Stop leaves the host running',async t=>{
 const {client,relay}=await fixture(t),host=await client({name:'Teacher'}),other=await client({name:'Other room'});
 const state=snapshot();state.state.selected=['ZA:Stomach'];state.state.contextOpacity=.18;state.state.sections=[{enabled:true,axis:'axial',position:.5,flip:true}];
 host.send({type:'rename',id:host.joined.self.id,name:'Instructor'});await host.next('self');
 const guest=await client({role:'guest',roomId:host.joined.roomId});assert.equal(guest.joined.host.name,'Instructor');await play(host,guest);host.send({type:'snapshot',snapshot:state});assert.deepEqual((await guest.next('snapshot')).snapshot,state);
 host.send({type:'rename',id:guest.joined.self.id,name:'Learner'});assert.equal((await guest.next('self')).self.name,'Learner');
 assert.equal(other.queue.some(message=>message.type==='snapshot'),false);
 guest.send({type:'stop'});await guest.next('ended');assert.ok(relay.rooms.has(host.joined.roomId));assert.equal(relay.rooms.get(host.joined.roomId).guests.size,0);
});
test('Host Stop ends every guest session and returns the host to an empty lobby',async t=>{
 const {client,relay}=await fixture(t),host=await client(),a=await client({role:'guest',roomId:host.joined.roomId}),b=await client({role:'guest',roomId:host.joined.roomId});
 await play(host,a,b);host.send({type:'stop'});assert.match((await a.next('ended')).reason,/host stopped/);await b.next('ended');await host.next('stopped');const room=relay.rooms.get(host.joined.roomId);assert.equal(room.playing,false);assert.equal(room.guests.size,0);assert.equal(room.participants.size,0);assert.equal(room.snapshot,null);assert.equal(room.pose,null);
});
test('only the host can remove a guest; other guests keep receiving slides and poses',async t=>{
 const {client,relay,url,origin,clients}=await fixture(t),host=await client(),a=await client({role:'guest',roomId:host.joined.roomId}),b=await client({role:'guest',roomId:host.joined.roomId});
 await play(host,a,b);
 const room=relay.rooms.get(host.joined.roomId);
 const roster=await host.next('roster',message=>message.guests.length===2&&message.guests.every(guest=>guest.inSession));assert.equal(roster.guests.length,2);
 a.send({type:'kick',id:b.joined.self.id});a.send({type:'rename',id:a.joined.self.id,name:'Still here'});await a.next('self');assert.equal(room.guests.size,2);
 host.send({type:'kick',id:a.joined.self.id});assert.match((await a.next('ended')).reason,/removed/);
 const remaining=await host.next('roster',message=>message.guests.length===1&&message.guests[0].id===b.joined.self.id);assert.equal(remaining.guests[0].inSession,true);
 assert.equal(room.playing,true);assert.equal(room.participants.has(a.clientToken),false);assert.equal(room.participants.has(b.clientToken),true);
 const state=snapshot();state.advanced={document:{schemaVersion:1,type:'slides',views:[{id:'slide',name:'Slide',url:'https://example.test/?model=male-detail',model:'male-detail',addedAt:'2026-09-28'}]},slideIndex:0};
 host.send({type:'snapshot',snapshot:state});assert.deepEqual((await b.next('snapshot')).snapshot,state);
 host.send({type:'pose',pose});assert.deepEqual((await b.next('pose')).pose,pose);
 const rejected=new WebSocket(url,{origin,rejectUnauthorized:false}),messages=[];clients.push(rejected);rejected.on('error',()=>{});
 rejected.on('message',raw=>messages.push(JSON.parse(raw.toString())));await once(rejected,'open');
 const closed=once(rejected,'close');rejected.send(JSON.stringify({type:'join',role:'guest',roomId:room.id,clientToken:a.clientToken}));
 assert.equal((await closed)[0],4003);assert.ok(messages.some(message=>message.type==='ended'&&/removed/.test(message.reason)));
 host.send({type:'stop'});await host.next('stopped');assert.equal(room.excluded.size,0);
 const next=await client({role:'guest',roomId:room.id,clientToken:a.clientToken});assert.equal(next.joined.playing,false);
});
test('a brief host transport interruption resumes the same room and preserves guest names',async t=>{
 const {client}=await fixture(t),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId,name:'Student'});
 await play(host,guest);
 host.socket.close();await guest.next('roster',message=>!message.host.online);
 const resumed=await client({roomId:host.joined.roomId,clientToken:host.clientToken});assert.equal(resumed.joined.roomId,host.joined.roomId);assert.equal(resumed.joined.playing,true);
 const roster=await guest.next('roster',message=>message.host.online);assert.equal(roster.guests[0].name,'Student');
});
test('an absent host ends the session after the reconnect grace period',async t=>{
 const {client,relay}=await fixture(t,30),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId});
 host.socket.terminate();assert.match((await guest.next('ended')).reason,/disconnected/);assert.equal(relay.rooms.size,0);
});
test('the relay rejects unrelated browser origins and malformed snapshots',async t=>{
 const {url,client,relay,clients}=await fixture(t);
 const forbidden=new WebSocket(url,{origin:'https://unrelated.example',rejectUnauthorized:false});clients.push(forbidden);
 const error=await new Promise(resolve=>forbidden.once('error',resolve));assert.match(error.message,/403/);
 const host=await client();await play(host);const bad=snapshot();bad.state.explode='invalid';host.send({type:'snapshot',snapshot:bad});assert.match((await host.next('error')).message,/could not be shared/);assert.equal(relay.rooms.get(host.joined.roomId).snapshot,null);
});
test('guests can wait and be counted before Play without receiving viewer actions',async t=>{
 const {client,relay}=await fixture(t),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId});
 assert.equal(host.joined.playing,false);assert.equal(guest.joined.playing,false);
 const roster=await host.next('roster',message=>message.guests.length===1);assert.equal(roster.guests.filter(peer=>peer.online).length,1);
 host.send({type:'snapshot',snapshot:snapshot()});host.send({type:'pose',pose});host.send({type:'rename',id:host.joined.self.id,name:'Ready host'});await host.next('self');
 const room=relay.rooms.get(host.joined.roomId);assert.equal(room.snapshot,null);assert.equal(room.pose,null);assert.equal(guest.queue.some(message=>['snapshot','pose'].includes(message.type)),false);
 guest.send({type:'stop'});await guest.next('ended');assert.equal((await host.next('roster',message=>message.guests.length===0)).guests.length,0);
});
test('Play captures online waiting guests, late guests wait, and captured guests can reconnect',async t=>{
 const {client,relay}=await fixture(t),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId});await play(host,guest);
 const state=snapshot();host.send({type:'snapshot',snapshot:state});await guest.next('snapshot');
 const late=await client({role:'guest',roomId:host.joined.roomId});assert.equal(late.joined.playing,false);assert.equal(late.queue.some(message=>['snapshot','pose','playing'].includes(message.type)),false);
 assert.deepEqual([...relay.rooms.get(host.joined.roomId).participants],[guest.clientToken]);
 guest.socket.close();await host.next('roster',message=>message.guests.some(peer=>peer.id===guest.joined.self.id&&!peer.online));
 const resumed=await client({role:'guest',roomId:host.joined.roomId,clientToken:guest.clientToken});assert.equal(resumed.joined.playing,true);assert.deepEqual((await resumed.next('snapshot')).snapshot,state);
 host.send({type:'stop'});await resumed.next('ended');await late.next('ended');await host.next('stopped');
 const next=await client({role:'guest',roomId:host.joined.roomId});assert.equal(next.joined.playing,false);await play(host,next);assert.equal(relay.rooms.get(host.joined.roomId).participants.has(next.clientToken),true);
});
test('normalized frames survive WSS snapshot, pose, and guest reconnection without viewport pixels',async t=>{
 const {client}=await fixture(t),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId});await play(host,guest);
 const framed={...pose,camera:[.5,.5,1,.5,.5,.5,.5,.5,.4,.5,1,.5]},state=snapshot();state.pose=framed;state.state.camera=framed.camera;
 assert.ok(validPose(framed));assert.ok(validSnapshot(state));host.send({type:'snapshot',snapshot:state});assert.deepEqual((await guest.next('snapshot')).snapshot,state);
 const moved={...framed,camera:[...framed.camera.slice(0,8),.35,...framed.camera.slice(9)]};host.send({type:'pose',pose:moved});assert.deepEqual((await guest.next('pose')).pose,moved);
 guest.socket.close();await host.next('roster',message=>message.guests.some(peer=>!peer.online));const resumed=await client({role:'guest',roomId:host.joined.roomId,clientToken:guest.clientToken});assert.deepEqual((await resumed.next('pose')).pose,moved);
 assert.equal(validPose({...framed,camera:[...framed.camera.slice(0,8),1,...framed.camera.slice(9)]}),false);
});

test('biological views, search and transcript branches are shared with session guests',async t=>{
 const {client}=await fixture(t),host=await client(),guest=await client({role:'guest',roomId:host.joined.roomId});await play(host,guest);
 const value=snapshot();value.hierarchy='guest:genes';Object.assign(value.state,{guestQuery:'HGNC:399',guestView:'types',guestExtensions:['HGNC:399']});assert.ok(validSnapshot(value));
 host.send({type:'snapshot',snapshot:value});assert.deepEqual((await guest.next('snapshot')).snapshot.state,value.state);
 assert.equal(validSnapshot({...value,state:{...value.state,guestQuery:123}}),false);
});
