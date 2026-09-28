import {DurableObject} from 'cloudflare:workers';
import {CONNECT_PATH,MAX_MESSAGE_BYTES,validSnapshot,validPose} from '../shared/connect-protocol.mjs';

const validCode=value=>typeof value==='string'&&/^[\w-]{32}$/.test(value);
const cleanName=(name,ip)=>typeof name==='string'&&name.trim()?name.trim().slice(0,80):ip;
const GRACE_MS=60_000,JOIN_MS=10_000,MAX_GUESTS=200;
const allowedOrigin=(request,env)=>{
 const value=request.headers.get('Origin');if(!value)return false;
 try{
  const origin=new URL(value),relay=new URL(request.url);
  if(value!==origin.origin)return false;
  return origin.origin===relay.origin||(origin.protocol==='https:'&&(env.ALLOWED_VIEWER_ORIGINS??'').split(',').map(value=>value.trim()).includes(value));
 }catch{return false;}
};
const json=(data,headers={})=>Response.json(data,{headers:{'Cache-Control':'no-store',...headers}});

export default {
 async fetch(request,env){
  const url=new URL(request.url);
  if(url.pathname==='/'||url.pathname==='/health')return json({service:'human-atlas-connect',status:'ok',protocol:1});
  if(url.pathname!==CONNECT_PATH&&url.pathname!==`${CONNECT_PATH}/identity`)return new Response('Not found',{status:404});
  if(!allowedOrigin(request,env))return new Response('Forbidden',{status:403});
  const headers={'Access-Control-Allow-Origin':request.headers.get('Origin'),'Vary':'Origin'};
  if(request.method==='OPTIONS')return new Response(null,{status:204,headers:{...headers,'Access-Control-Allow-Methods':'GET','Access-Control-Allow-Headers':'Content-Type'}});
  if(request.method!=='GET')return new Response('Use GET',{status:405});
  if(url.pathname.endsWith('/identity'))return json({ip:request.headers.get('CF-Connecting-IP')??'Unknown IP',origin:url.origin},headers);
  if(request.headers.get('Upgrade')?.toLowerCase()!=='websocket')return new Response('Use a WebSocket connection',{status:426});
  const room=url.searchParams.get('room');
  if(!validCode(room))return new Response('Invalid session code',{status:400});
  return env.ROOMS.getByName(room).fetch(request);
 }
};

// One SQLite-backed object per invitation. Idle WebSockets can hibernate:
// no interval runs, and camera updates use socket attachments, not SQL writes.
export class AtlasRoom extends DurableObject {
 constructor(ctx,env){
  super(ctx,env);
  this.sql=ctx.storage.sql;
  this.sql.exec('CREATE TABLE IF NOT EXISTS atlas_state (key TEXT PRIMARY KEY, value TEXT NOT NULL)');
  const read=key=>{const row=this.sql.exec('SELECT value FROM atlas_state WHERE key = ?',key).toArray()[0];return row?JSON.parse(row.value):null;};
  this.meta=read('room');this.snapshot=read('snapshot');this.pose=read('pose');
  this.sockets=new Map();
  for(const socket of ctx.getWebSockets()){
   const attachment=socket.deserializeAttachment();
   if(attachment?.joined&&socket.readyState===WebSocket.OPEN){
    this.sockets.set(attachment.token,socket);
    if(attachment.role==='host'&&attachment.pose)this.pose=attachment.pose;
   }
  }
  ctx.setWebSocketAutoResponse(new WebSocketRequestResponsePair('ping','pong'));
  // A deployment/runtime restart can lose transports before a close event.
  if(this.meta&&!this.sockets.has(this.meta.host.token)&&!this.meta.hostMissingUntil){
   this.meta.hostMissingUntil=Date.now()+GRACE_MS;this.save();
   ctx.blockConcurrencyWhile(()=>this.schedule());
  }
 }
 put(key,value){this.sql.exec('INSERT INTO atlas_state (key,value) VALUES (?,?) ON CONFLICT(key) DO UPDATE SET value=excluded.value',key,JSON.stringify(value));}
 save(){this.put('room',this.meta);}
 send(socket,message){if(socket?.readyState===WebSocket.OPEN)try{socket.send(JSON.stringify(message));}catch{socket.close(1013,'Connection is too slow');}}
 peer(identity){return {id:identity.id,ip:identity.ip,name:identity.name,online:this.sockets.has(identity.token),inSession:identity.token!==this.meta.host.token&&this.meta.participants.includes(identity.token)};}
 roster(){
  if(!this.meta)return;
  const message={type:'roster',roomId:this.meta.id,host:this.peer(this.meta.host),guests:this.meta.guests.map(guest=>this.peer(guest))};
  for(const socket of this.sockets.values())this.send(socket,message);
 }
 broadcast(message){for(const token of this.meta.participants)this.send(this.sockets.get(token),message);}
 clearView(){
  this.snapshot=this.pose=null;this.sql.exec("DELETE FROM atlas_state WHERE key IN ('snapshot','pose')");
  const host=this.sockets.get(this.meta.host.token);
  if(host){const attachment=host.deserializeAttachment();host.serializeAttachment({...attachment,pose:null});}
 }
 detach(token,reason,code=1000){
  const socket=this.sockets.get(token);this.sockets.delete(token);
  if(socket){this.send(socket,{type:'ended',reason});const attachment=socket.deserializeAttachment();socket.serializeAttachment({...attachment,joined:false,pose:null});socket.close(code,reason);}
 }
 async schedule(){
  const deadlines=this.ctx.getWebSockets().flatMap(socket=>{const attachment=socket.deserializeAttachment();return socket.readyState===WebSocket.OPEN&&!attachment?.joined&&attachment?.joinUntil?[attachment.joinUntil]:[];});
  if(this.meta?.hostMissingUntil)deadlines.push(this.meta.hostMissingUntil);
  if(deadlines.length)await this.ctx.storage.setAlarm(Math.max(Date.now()+1,Math.min(...deadlines)));
  else await this.ctx.storage.deleteAlarm();
 }
 async expire(){
  if(this.meta?.hostMissingUntil&&this.meta.hostMissingUntil<=Date.now()){
   for(const token of [...this.sockets.keys()])this.detach(token,'The host disconnected.');
   this.meta=this.snapshot=this.pose=null;await this.ctx.storage.deleteAll();return true;
  }
  return false;
 }
 async fetch(request){
  await this.expire();
  const pair=new WebSocketPair(),[client,socket]=Object.values(pair),ip=request.headers.get('CF-Connecting-IP')??'Unknown IP';
  this.ctx.acceptWebSocket(socket);
  socket.serializeAttachment({ip,roomId:new URL(request.url).searchParams.get('room'),joined:false,joinUntil:Date.now()+JOIN_MS});
  this.send(socket,{type:'hello',ip});await this.schedule();
  return new Response(null,{status:101,webSocket:client});
 }
 async join(socket,message,attachment){
  if(!validCode(message.clientToken)||!['host','guest'].includes(message.role)||message.roomId!==attachment.roomId){socket.close(1008,'Invalid session');return;}
  await this.expire();
  const token=message.clientToken,hosting=message.role==='host';
  if(hosting&&this.meta&&this.meta.host.token!==token){this.send(socket,{type:'error',message:'This session belongs to another host.'});socket.close(1008,'Wrong host');return;}
  if(!hosting&&(!this.meta||this.meta.host.token===token)){
   this.send(socket,{type:'error',message:this.meta?'Use another browser or tab to join your own session.':'The host is not running this session. Ask for their current invitation.'});socket.close(4004,'Host unavailable');return;
  }
  if(!hosting&&this.meta.excluded.includes(token)){this.send(socket,{type:'ended',reason:'You were removed from this session by the host.'});socket.close(4003,'Removed by host');return;}
  let identity=hosting?this.meta?.host:this.meta.guests.find(guest=>guest.token===token);
  if(!hosting&&!identity&&this.meta.guests.length>=MAX_GUESTS){this.send(socket,{type:'error',message:'This session has reached its guest limit.'});socket.close(1008,'Session full');return;}
  identity=identity??{token,id:crypto.randomUUID(),ip:attachment.ip,name:attachment.ip};
  identity.ip=attachment.ip;identity.name=cleanName(message.name,identity.name);
  const old=this.sockets.get(token);if(old&&old!==socket){old.serializeAttachment({...old.deserializeAttachment(),joined:false});old.close(4001,'Connected in another window');}
  if(hosting){
   this.meta??={id:attachment.roomId,host:identity,guests:[],participants:[],excluded:[],playing:false,hostMissingUntil:null};
   this.meta.host=identity;this.meta.hostMissingUntil=null;
  }else if(!this.meta.guests.some(guest=>guest.token===token))this.meta.guests.push(identity);
  this.sockets.set(token,socket);socket.serializeAttachment({...attachment,joined:true,token,role:message.role,pose:hosting?this.pose:null,windowStarted:Date.now(),count:0});
  this.save();await this.schedule();
  const playing=this.meta.playing&&(hosting||this.meta.participants.includes(token));
  this.send(socket,{type:'joined',role:message.role,roomId:this.meta.id,playing,self:this.peer(identity),host:this.peer(this.meta.host)});this.roster();
  if(!hosting&&playing){if(this.snapshot)this.send(socket,{type:'snapshot',snapshot:this.snapshot});if(this.pose)this.send(socket,{type:'pose',pose:this.pose});}
 }
 async webSocketMessage(socket,raw){
  if(typeof raw!=='string'){socket.close(1003,'Use JSON messages');return;}
  if(new TextEncoder().encode(raw).byteLength>MAX_MESSAGE_BYTES){socket.close(1009,'Message too large');return;}
  let attachment=socket.deserializeAttachment();if(!attachment)return;
  const now=Date.now();if(!attachment.windowStarted||now-attachment.windowStarted>1000){attachment.windowStarted=now;attachment.count=0;}
  if(++attachment.count>90){socket.close(1008,'Too many messages');return;}
  socket.serializeAttachment(attachment);
  let message;try{message=JSON.parse(raw);}catch{this.send(socket,{type:'error',message:'Invalid connection message.'});return;}
  if(!message||typeof message!=='object')return;
  if(message.type==='join'){if(!attachment.joined)await this.join(socket,message,attachment);return;}
  if(!attachment.joined||!this.meta||this.sockets.get(attachment.token)!==socket||await this.expire())return;
  const token=attachment.token,hosting=this.meta.host.token===token;
  const remove=target=>{this.meta.guests=this.meta.guests.filter(guest=>guest.token!==target);this.meta.participants=this.meta.participants.filter(value=>value!==target);};
  if(message.type==='stop'){
   if(hosting){
    for(const guest of this.meta.guests)this.detach(guest.token,'The host stopped the session.');
    this.meta.playing=false;this.meta.guests=[];this.meta.participants=[];this.meta.excluded=[];this.clearView();this.send(socket,{type:'stopped'});
   }else{remove(token);this.detach(token,'You left the session.');}
   this.save();this.roster();return;
  }
  if(message.type==='rename'){
   const identity=message.id===this.meta.host.id?this.meta.host:this.meta.guests.find(guest=>guest.id===message.id);
   if(identity&&(identity.token===token||hosting&&identity!==this.meta.host)){identity.name=cleanName(message.name,identity.ip);this.save();this.send(this.sockets.get(identity.token),{type:'self',self:this.peer(identity)});this.roster();}return;
  }
  if(!hosting)return;
  if(message.type==='kick'){
   const guest=this.meta.guests.find(guest=>guest.id===message.id);if(!guest)return;
   this.meta.excluded.push(guest.token);remove(guest.token);this.detach(guest.token,'You were removed from this session by the host.',4003);this.save();this.roster();return;
  }
  if(message.type==='play'){
   if(!this.meta.playing){this.meta.playing=true;this.meta.participants=this.meta.guests.filter(guest=>this.sockets.has(guest.token)).map(guest=>guest.token);this.save();this.send(socket,{type:'playing'});this.broadcast({type:'playing'});this.roster();}return;
  }
  if(!this.meta.playing)return;
  if(message.type==='snapshot'&&validSnapshot(message.snapshot)){
   this.snapshot=message.snapshot;this.put('snapshot',this.snapshot);
   if(this.snapshot.pose){this.pose=this.snapshot.pose;socket.serializeAttachment({...attachment,pose:this.pose});}
   this.broadcast({type:'snapshot',snapshot:this.snapshot});
  }else if(message.type==='pose'&&validPose(message.pose)){
   this.pose=message.pose;socket.serializeAttachment({...attachment,pose:this.pose});this.broadcast({type:'pose',pose:this.pose});
  }else this.send(socket,{type:'error',message:'The viewer state could not be shared.'});
 }
 async webSocketClose(socket,code,reason){
  const attachment=socket.deserializeAttachment();
  if(attachment?.joined&&this.sockets.get(attachment.token)===socket){
   this.sockets.delete(attachment.token);
   if(this.meta?.host.token===attachment.token){this.meta.hostMissingUntil=Date.now()+GRACE_MS;if(this.pose)this.put('pose',this.pose);this.save();}
   this.roster();
  }
  if(socket.readyState===WebSocket.OPEN)socket.close(code===1005?1000:code,reason);
  await this.schedule();
 }
 async webSocketError(socket){socket.close(1011,'Connection failed');await this.webSocketClose(socket,1011,'Connection failed');}
 async alarm(){
  await this.expire();
  for(const socket of this.ctx.getWebSockets()){const attachment=socket.deserializeAttachment();if(!attachment?.joined&&attachment?.joinUntil<=Date.now())socket.close(1008,'Choose a role first');}
  await this.schedule();
 }
}
