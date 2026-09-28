import assert from 'node:assert/strict';
import {randomBytes} from 'node:crypto';
import {once} from 'node:events';
import WebSocket from 'ws';
import {connectionSocketUrl} from '../shared/connect-address.mjs';

const relay=new URL(process.argv[2]??'');
assert.equal(relay.protocol,'https:','Pass the deployed relay HTTPS address.');
const health=await fetch(new URL('/health',relay),{signal:AbortSignal.timeout(15_000)});
assert.equal(health.status,200);assert.equal((await health.json()).status,'ok');
const room=randomBytes(24).toString('base64url'),sockets=[];
async function client(role){
 const socket=new WebSocket(connectionSocketUrl(relay.origin,room),{origin:'https://ctzurcanu.github.io'}),queue=[],waiters=[];sockets.push(socket);
 socket.on('error',()=>{});
 socket.on('message',raw=>{
  const message=JSON.parse(raw.toString()),index=waiters.findIndex(waiter=>waiter.type===message.type&&waiter.predicate(message));
  if(index<0)queue.push(message);else{const waiter=waiters.splice(index,1)[0];clearTimeout(waiter.timer);waiter.resolve(message);}
 });
 const next=(type,predicate=()=>true)=>{
  const index=queue.findIndex(message=>message.type===type&&predicate(message));if(index>=0)return Promise.resolve(queue.splice(index,1)[0]);
  return new Promise((resolve,reject)=>{const waiter={type,predicate,resolve,timer:setTimeout(()=>reject(new Error(`Timed out waiting for ${type}`)),15_000)};waiters.push(waiter);});
 };
 await once(socket,'open');await next('hello');
 const send=message=>socket.send(JSON.stringify(message));send({type:'join',role,roomId:room,clientToken:randomBytes(24).toString('base64url'),name:`Relay verification ${role}`});
 const joined=await next('joined');assert.equal(joined.roomId,room);return {socket,next,send,joined};
}
try{
 const host=await client('host'),guest=await client('guest');
 host.send({type:'play'});await host.next('playing');await guest.next('playing');
 const pose={model:'male-detail',camera:[1,1,3,0,.8,0,0,0],up:[0,1,0],rotation:[0,0,0,1]};
 const snapshot={version:1,model:'male-detail',hierarchy:'systems',guestSources:[],state:{selected:['ZA:Stomach'],visible:['digestive'],explode:0,isolate:false,view:'three-quarter',rotate:false,reset:0},choice:null,ui:{panel:'advanced',frontPanel:'advanced',details:false,layersVisible:true,mobileLayersOpen:false,addSelection:false,query:''},covering:[],pose,advanced:null};
 host.send({type:'snapshot',snapshot});assert.deepEqual((await guest.next('snapshot')).snapshot,snapshot);
 host.send({type:'pose',pose});assert.deepEqual((await guest.next('pose')).pose,pose);
 const late=await client('guest');assert.equal(late.joined.playing,false);
 host.send({type:'admit',id:late.joined.self.id});await late.next('playing');
 assert.deepEqual((await late.next('snapshot')).snapshot,snapshot);assert.deepEqual((await late.next('pose')).pose,pose);await host.next('snapshot-request');
 await host.next('roster',message=>message.guests.length===2&&message.guests.every(peer=>peer.inSession));
 const waiting=await Promise.all([client('guest'),client('guest')]);assert.ok(waiting.every(peer=>!peer.joined.playing));
 host.send({type:'admit-all'});await Promise.all(waiting.map(async peer=>{await peer.next('playing');assert.deepEqual((await peer.next('snapshot')).snapshot,snapshot);assert.deepEqual((await peer.next('pose')).pose,pose);}));await host.next('snapshot-request');
 await host.next('roster',message=>message.guests.length===4&&message.guests.every(peer=>peer.inSession));
 const moved={...pose,camera:[2,1,3,0,.8,0,0,0]};host.send({type:'pose',pose:moved});await Promise.all([guest,late,...waiting].map(async peer=>assert.deepEqual((await peer.next('pose')).pose,moved)));
 host.send({type:'stop'});assert.match((await guest.next('ended')).reason,/host stopped/);await host.next('stopped');
 await late.next('ended');
 await Promise.all(waiting.map(peer=>peer.next('ended')));
 console.log('Live HTTPS/WSS relay verified: health, Host, Guest, Play, individual and bulk waiting-guest admission, anatomy state, camera, and Stop.');
}finally{for(const socket of sockets)socket.close(1000,'Verification complete');}
