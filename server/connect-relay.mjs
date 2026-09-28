import {randomBytes,randomUUID} from 'node:crypto';
import {isIP} from 'node:net';
import {WebSocketServer,WebSocket} from 'ws';
import {CONNECT_PATH,MAX_MESSAGE_BYTES,validSnapshot,validPose} from '../shared/connect-protocol.mjs';
import {hostAddresses} from './connect-address.mjs';
import {localOnlyHost} from '../shared/connect-address.mjs';

const code=()=>randomBytes(24).toString('base64url');
const validToken=value=>typeof value==='string'&&/^[\w-]{32}$/.test(value);
const cleanName=(name,fallback)=>typeof name==='string'&&name.trim()?name.trim().slice(0,80):fallback;
export function requestIp(request){
 const forwarded=String(request.headers['x-forwarded-for']??'').split(',')[0].trim();
 const address=isIP(forwarded)?forwarded:request.socket.remoteAddress??'Unknown IP';
 const ip=address.startsWith('::ffff:')?address.slice(7):address;
 return localOnlyHost(ip)?hostAddresses()[0]??'Unknown IP':ip;
}
export function attachConnectRelay(server,{path=CONNECT_PATH,graceMs=60_000}={}){
 const sockets=new WebSocketServer({noServer:true,maxPayload:MAX_MESSAGE_BYTES,perMessageDeflate:false});
 const rooms=new Map(),identities=new Map();
 const send=(socket,message)=>{if(socket?.readyState===WebSocket.OPEN){if(socket.bufferedAmount>MAX_MESSAGE_BYTES*2){socket.close(1013,'Connection is too slow');return;}socket.send(JSON.stringify(message));}};
 const peer=(identity,online,inSession=false)=>({id:identity.id,ip:identity.ip,name:identity.name,online,inSession});
 const roster=room=>{
  const message={type:'roster',roomId:room.id,host:peer(room.host,!!room.host.socket),guests:[...room.guests.values()].map(identity=>peer(identity,!!identity.socket,room.participants.has(identity.token)))};
  send(room.host.socket,message);for(const guest of room.guests.values())send(guest.socket,message);
 };
 const endRoom=(room,reason)=>{
  if(room.timer)clearTimeout(room.timer);
  rooms.delete(room.id);
  for(const identity of [room.host,...room.guests.values()]){send(identity.socket,{type:'ended',reason});identity.room=null;}
  room.guests.clear();room.snapshot=null;room.pose=null;
 };
 const stopPlay=room=>{
  room.playing=false;room.snapshot=null;room.pose=null;room.participants.clear();room.excluded.clear();
  for(const guest of room.guests.values()){send(guest.socket,{type:'ended',reason:'The host stopped the session.'});guest.room=null;}
  room.guests.clear();send(room.host.socket,{type:'stopped'});roster(room);
 };
 const leave=identity=>{
  const room=identity.room;if(!room)return;
  if(room.host===identity)stopPlay(room);
  else{room.guests.delete(identity.token);room.participants.delete(identity.token);identity.room=null;roster(room);}
 };
 const upgrade=(request,socket,head)=>{
  if((request.url??'').split('?')[0]!==path)return;
  // A room invitation is the access credential. Only the viewer's own origin
  // may use this relay; unrelated websites cannot open a room in a browser.
  let allowed=false;
  try{
   const origin=new URL(String(request.headers.origin));
   const host=String(request.headers['x-forwarded-host']??request.headers.host??'').split(',')[0].trim();
   const forwardedProtocol=String(request.headers['x-forwarded-proto']??'').split(',')[0].trim();
   const secure=request.socket.encrypted||forwardedProtocol==='https';
   const loopback=['localhost','127.0.0.1','[::1]'].includes(origin.hostname);
   allowed=origin.host===host&&(secure?origin.protocol==='https:':loopback&&origin.protocol==='http:');
  }catch{/* Reject missing or malformed Origin. */}
  if(!allowed){socket.write('HTTP/1.1 403 Forbidden\r\nConnection: close\r\n\r\n');socket.destroy();return;}
  sockets.handleUpgrade(request,socket,head,ws=>sockets.emit('connection',ws,request));
 };
 server.on('upgrade',upgrade);
 sockets.on('connection',(socket,request)=>{
  let identity=null,alive=true,windowStarted=Date.now(),messageCount=0;
  const ip=requestIp(request);
  send(socket,{type:'hello',ip});
  const joinTimer=setTimeout(()=>socket.close(1008,'Choose a role first'),10000);joinTimer.unref?.();
  socket.on('pong',()=>{alive=true;});socket.on('error',()=>{});
  socket.on('message',(raw,binary)=>{
   if(binary){socket.close(1003,'Use JSON messages');return;}
   if(Date.now()-windowStarted>1000){windowStarted=Date.now();messageCount=0;}
   if(++messageCount>90){socket.close(1008,'Too many messages');return;}
   let message;try{message=JSON.parse(raw.toString());}catch{send(socket,{type:'error',message:'Invalid connection message.'});return;}
   if(!message||typeof message!=='object')return;
   if(message.type==='join'){
    if(identity||!validToken(message.clientToken)||!['host','guest'].includes(message.role))return;
    const existing=identities.get(message.clientToken);
    identity=existing??{id:randomUUID(),token:message.clientToken,ip,name:cleanName(message.name,ip),socket:null,room:null};
    if(identity.socket&&identity.socket!==socket)identity.socket.close(4001,'Connected in another window');
    identity.socket=socket;identity.ip=ip;identity.name=cleanName(message.name,identity.name);identities.set(identity.token,identity);clearTimeout(joinTimer);
    let room;
    if(message.role==='host'){
     room=rooms.get(message.roomId);
     if(room&&room.host.token!==identity.token){send(socket,{type:'error',message:'This session belongs to another host.'});identity.socket=null;identity=null;socket.close(1008);return;}
     if(!room){room={id:code(),host:identity,guests:new Map(),participants:new Set(),excluded:new Set(),playing:false,snapshot:null,pose:null,timer:null};rooms.set(room.id,room);}
     if(room.timer){clearTimeout(room.timer);room.timer=null;}
    }else{
     room=rooms.get(message.roomId);
     if(!room){send(socket,{type:'error',message:'The host is not running this session. Ask for their current invitation.'});identity.socket=null;identity=null;socket.close(4004,'Host unavailable');return;}
     if(room.host.token===identity.token){send(socket,{type:'error',message:'Use another browser or tab to join your own session.'});identity.socket=null;identity=null;socket.close(1008);return;}
     if(room.excluded.has(identity.token)){send(socket,{type:'ended',reason:'You were removed from this session by the host.'});identity.socket=null;identity=null;socket.close(4003,'Removed by host');return;}
     room.guests.set(identity.token,identity);
    }
    identity.room=room;
    const playing=room.playing&&(message.role==='host'||room.participants.has(identity.token));
    send(socket,{type:'joined',role:message.role,roomId:room.id,playing,self:peer(identity,true,message.role==='guest'&&playing),host:peer(room.host,!!room.host.socket)});roster(room);
    if(message.role==='guest'&&playing&&room.snapshot)send(socket,{type:'snapshot',snapshot:room.snapshot});
    if(message.role==='guest'&&playing&&room.pose)send(socket,{type:'pose',pose:room.pose});
    return;
   }
   const room=identity?.room;if(!room||identity.socket!==socket)return;
   if(message.type==='stop'){const hosting=room.host===identity;leave(identity);if(!hosting)send(socket,{type:'ended',reason:'You left the session.'});return;}
   if(message.type==='rename'){
    const target=message.id===identity.id?identity:room.host===identity?[...room.guests.values()].find(guest=>guest.id===message.id):null;
    if(target){target.name=cleanName(message.name,target.ip);send(target.socket,{type:'self',self:peer(target,true)});roster(room);}return;
   }
   if(room.host!==identity)return; // Guests can never control the host or other guests.
   if(message.type==='kick'){
    const guest=[...room.guests.values()].find(guest=>guest.id===message.id);if(!guest)return;
    room.excluded.add(guest.token);leave(guest);send(guest.socket,{type:'ended',reason:'You were removed from this session by the host.'});guest.socket?.close(4003,'Removed by host');return;
   }
   if(message.type==='play'){
    if(!room.playing){room.playing=true;room.participants=new Set([...room.guests.values()].filter(guest=>guest.socket).map(guest=>guest.token));send(socket,{type:'playing'});for(const guest of room.guests.values())if(room.participants.has(guest.token))send(guest.socket,{type:'playing'});roster(room);}
    return;
   }
   if(!room.playing)return;
   if(message.type==='snapshot'&&validSnapshot(message.snapshot)){
    room.snapshot=message.snapshot;
    if(message.snapshot.pose)room.pose=message.snapshot.pose;
    for(const guest of room.guests.values())if(room.participants.has(guest.token))send(guest.socket,{type:'snapshot',snapshot:room.snapshot});
   }else if(message.type==='pose'&&validPose(message.pose)){
    room.pose=message.pose;for(const guest of room.guests.values())if(room.participants.has(guest.token))send(guest.socket,{type:'pose',pose:room.pose});
   }else send(socket,{type:'error',message:'The viewer state could not be shared.'});
  });
  socket.on('close',()=>{
   clearTimeout(joinTimer);
   if(!identity||identity.socket!==socket)return;
   identity.socket=null;
   const room=identity.room;
   if(room){roster(room);if(room.host===identity){room.timer=setTimeout(()=>endRoom(room,'The host disconnected.'),graceMs);room.timer.unref?.();}}
   setTimeout(()=>{if(!identity.socket&&!identity.room)identities.delete(identity.token);},graceMs).unref?.();
  });
  socket.isAlive=()=>alive;socket.markPing=()=>{alive=false;};
 });
 const heartbeat=setInterval(()=>{for(const socket of sockets.clients){if(!socket.isAlive()){socket.terminate();continue;}socket.markPing();socket.ping();}},20000);heartbeat.unref?.();
 const dispose=()=>{clearInterval(heartbeat);server.off('upgrade',upgrade);for(const room of rooms.values())if(room.timer)clearTimeout(room.timer);for(const socket of sockets.clients)socket.terminate();sockets.close();rooms.clear();identities.clear();};
 server.once('close',dispose);
 return {dispose,rooms};
}
