import {useCallback,useEffect,useRef,useState} from 'react';
import {CONNECT_PATH,validSnapshot,validPose} from '../shared/connect-protocol.mjs';
import {guestInvitation,invitationUrl,localOnlyHost} from '../shared/connect-address.mjs';
import type {ConnectPeer,ConnectRole,ConnectSnapshot,ScenePose} from './connect-state';

type PreviousHost={id:string;ip:string;name:string;alias?:string;roomId:string;url?:string;lastConnected:string};
type Status='idle'|'connecting'|'ready'|'running'|'reconnecting';
const NAME_KEY='human-atlas-connect-name',HISTORY_KEY='human-atlas-connect-hosts';
const read=(key:string)=>{try{return localStorage.getItem(key);}catch{return null;}};
const write=(key:string,value:string)=>{try{localStorage.setItem(key,value);}catch{/* Keep the session usable without storage. */}};
const sessionToken=()=>{
 try{const saved=sessionStorage.getItem('human-atlas-connect-client');if(saved)return saved;}catch{/* Generate an in-memory identity. */}
 const token=btoa(String.fromCharCode(...crypto.getRandomValues(new Uint8Array(24)))).replace(/\+/g,'-').replace(/\//g,'_').replace(/=/g,'');
 try{sessionStorage.setItem('human-atlas-connect-client',token);}catch{/* Storage is optional. */}return token;
};
const previousHosts=():PreviousHost[]=>{
 try{const values:unknown=JSON.parse(read(HISTORY_KEY)??'[]');return Array.isArray(values)?values.filter((host):host is PreviousHost=>!!host&&typeof host==='object'&&['id','ip','name','roomId','lastConnected'].every(key=>typeof host[key]==='string')&&!localOnlyHost(host.ip)&&/^[\w-]{32}$/.test(host.roomId)).slice(0,20):[];}catch{return [];}
};
export function hostCode(input:string){
 return guestInvitation(input,location.href).code;
}
export function useAtlasConnection(snapshot:()=>ConnectSnapshot){
 const invite=new URLSearchParams(location.search).get('connect')??'';
 const [role,updateRole]=useState<ConnectRole>(invite?'guest':'host');
 const [playing,setPlaying]=useState(false);
 const [name,setName]=useState(()=>read(NAME_KEY)??'');
 const [ip,setIp]=useState('');
 const [hostInput,setHostInput]=useState(invite?location.href:'');
 const [hostAddress,updateHostAddress]=useState(()=>localOnlyHost(location.hostname)?'':location.host);
 const hostAddressEdited=useRef(false);
 const [status,setStatus]=useState<Status>('idle');
 const [error,setError]=useState(''),[notice,setNotice]=useState('');
 const [roomId,setRoomId]=useState('');
 const [self,setSelf]=useState<ConnectPeer|null>(null),[host,setHost]=useState<ConnectPeer|null>(null),[guests,setGuests]=useState<ConnectPeer[]>([]);
 const [hosts,setHosts]=useState(previousHosts);
 const [remoteSnapshot,setRemoteSnapshot]=useState<ConnectSnapshot|null>(null);
 const poseRef=useRef<ScenePose|null>(null),localPoseRef=useRef<ScenePose|null>(null);
 const socketRef=useRef<WebSocket|null>(null),snapshotRef=useRef(snapshot);snapshotRef.current=snapshot;
 const selfRef=useRef(self);selfRef.current=self;
 const nameRef=useRef(name);nameRef.current=name;
 const customNameRef=useRef(read(NAME_KEY)??'');
 const tokenRef=useRef<string>('');if(!tokenRef.current)tokenRef.current=sessionToken();
 const desired=useRef<{role:ConnectRole;roomId:string;playing:boolean}|null>(null);
 const retryRef=useRef<ReturnType<typeof setTimeout>|null>(null),snapshotTimer=useRef<ReturnType<typeof setTimeout>|null>(null),poseTimer=useRef<ReturnType<typeof setTimeout>|null>(null);
 const joinedRef=useRef(false),disposedRef=useRef(false),pendingStopRef=useRef(false),lastSnapshot=useRef('');
 const send=useCallback((message:unknown)=>{const socket=socketRef.current;if(joinedRef.current&&socket?.readyState===WebSocket.OPEN&&socket.bufferedAmount<1_048_576)socket.send(JSON.stringify(message));},[]);
 const rememberHost=useCallback((peer:ConnectPeer,id:string)=>setHosts(current=>{
  const old=current.find(item=>item.id===peer.id),url=new URL(location.href);url.searchParams.set('connect',id);
  const next=[{...old,id:peer.id,ip:peer.ip,name:peer.name,roomId:id,url:url.href,lastConnected:new Date().toISOString()},...current.filter(item=>item.id!==peer.id)].slice(0,20);
  write(HISTORY_KEY,JSON.stringify(next));return next;
 }),[]);
 const publishSnapshot=useCallback((force=false)=>{
  if(desired.current?.role!=='host'||!desired.current.playing||!joinedRef.current)return;
  if(snapshotTimer.current){if(!force)return;clearTimeout(snapshotTimer.current);}
  snapshotTimer.current=setTimeout(()=>{
   snapshotTimer.current=null;if(desired.current?.role!=='host'||!desired.current.playing||!joinedRef.current)return;
   const current=snapshotRef.current(),message={...current,pose:localPoseRef.current?.model===current.model?localPoseRef.current:current.pose};
   // Camera-only motion travels separately, without rebuilding anatomy state.
   const encoded=JSON.stringify({...message,pose:null,state:{...message.state,camera:undefined}});
   if(force||encoded!==lastSnapshot.current){lastSnapshot.current=encoded;send({type:'snapshot',snapshot:message});}
  },force?0:50);
 },[send]);
 const publishPose=useCallback((pose:ScenePose)=>{
  localPoseRef.current=pose;
  if(desired.current?.role!=='host'||!desired.current.playing||!joinedRef.current||poseTimer.current)return;
  poseTimer.current=setTimeout(()=>{poseTimer.current=null;if(desired.current?.role==='host'&&desired.current.playing&&localPoseRef.current)send({type:'pose',pose:localPoseRef.current});},50);
 },[send]);
 const finish=useCallback((message='')=>{
  desired.current=null;joinedRef.current=false;pendingStopRef.current=false;
  if(retryRef.current)clearTimeout(retryRef.current);
  if(snapshotTimer.current)clearTimeout(snapshotTimer.current);
  if(poseTimer.current)clearTimeout(poseTimer.current);
  retryRef.current=snapshotTimer.current=poseTimer.current=null;
  const socket=socketRef.current;socketRef.current=null;socket?.close();
  setStatus('idle');setPlaying(false);setRoomId('');setGuests([]);setHost(null);setNotice(message);
 },[]);
 const openSocketRef=useRef<(attempt?:number)=>void>(()=>{});
 openSocketRef.current=(attempt=0)=>{
  const session=desired.current;if(!session||disposedRef.current)return;
  const loopback=['localhost','127.0.0.1','[::1]'].includes(location.hostname);
  if(location.protocol!=='https:'&&!loopback){finish();setError('Open the viewer over HTTPS to use Connect.');return;}
  const url=new URL(CONNECT_PATH,location.origin);url.protocol=location.protocol==='https:'?'wss:':'ws:';
  const socket=new WebSocket(url);socketRef.current=socket;joinedRef.current=false;
  setStatus(attempt?'reconnecting':'connecting');
  const timeout=setTimeout(()=>{if(socketRef.current===socket&&!joinedRef.current){setError('Connect needs the live viewer server. Open your HTTPS/ngrok viewer or run npm run dev.');finish();}},10000);
  socket.onopen=()=>{if(socketRef.current===socket)socket.send(JSON.stringify({type:'join',clientToken:tokenRef.current,name:customNameRef.current,...session}));};
  socket.onmessage=event=>{
   if(socketRef.current!==socket)return;
   let message;try{message=JSON.parse(String(event.data));}catch{return;}
   switch(message.type){
    case 'hello':setIp(message.ip);if(!customNameRef.current)setName(message.ip);break;
    case 'joined':
     clearTimeout(timeout);joinedRef.current=true;session.roomId=message.roomId;setRoomId(message.roomId);setSelf(message.self);selfRef.current=message.self;setIp(message.self.ip);setName(message.self.name);setHost(message.host);setStatus(message.playing?'running':'ready');setError('');setNotice('');
     if(session.role==='host'){
      try{sessionStorage.setItem('human-atlas-connect-room',message.roomId);}catch{/* Optional resume. */}
      if(pendingStopRef.current){pendingStopRef.current=false;session.playing=false;setPlaying(false);setStatus('ready');send({type:'stop'});}
      else if(session.playing&&!message.playing)send({type:'play'});
      else{session.playing=!!message.playing;setPlaying(session.playing);if(session.playing)publishSnapshot(true);}
     }else{session.playing=!!message.playing;setPlaying(session.playing);rememberHost(message.host,message.roomId);}
     break;
    case 'playing':session.playing=true;setPlaying(true);setStatus('running');if(session.role==='host')publishSnapshot(true);break;
    case 'stopped':session.playing=false;setPlaying(false);setStatus('ready');setGuests([]);lastSnapshot.current='';break;
    case 'roster':setHost(message.host);setGuests(message.guests);if(session.role==='guest')rememberHost(message.host,message.roomId);break;
    case 'self':setSelf(message.self);selfRef.current=message.self;setName(message.self.name);customNameRef.current=message.self.name===message.self.ip?'':message.self.name;write(NAME_KEY,customNameRef.current);break;
    case 'snapshot':if(session.role==='guest'&&validSnapshot(message.snapshot)){poseRef.current=message.snapshot.pose;setRemoteSnapshot(message.snapshot);}break;
    case 'pose':if(session.role==='guest'&&validPose(message.pose))poseRef.current=message.pose;break;
    case 'ended':finish(message.reason);break;
    case 'error':setError(message.message);break;
   }
  };
  socket.onerror=()=>{if(socketRef.current===socket)setError('The live connection is unavailable. Retrying…');};
  socket.onclose=event=>{
   clearTimeout(timeout);if(socketRef.current!==socket||!desired.current||disposedRef.current)return;
   joinedRef.current=false;
   if([1008,4001,4004].includes(event.code)){finish();return;}
   setStatus('reconnecting');retryRef.current=setTimeout(()=>openSocketRef.current(attempt+1),Math.min(5000,500*2**Math.min(attempt,4)));
  };
 };
 const prepare=useCallback(async()=>{
  if(role==='host'&&!desired.current){let room='';try{room=sessionStorage.getItem('human-atlas-connect-room')??'';}catch{/* Storage is optional. */}desired.current={role:'host',roomId:room,playing:false};openSocketRef.current();}
  if(!ip)try{const response=await fetch(`${CONNECT_PATH}/identity`,{cache:'no-store'});if(!response.ok)return;const data=await response.json() as {ip?:unknown;origin?:unknown};if(typeof data.ip==='string'){setIp(data.ip);if(!customNameRef.current)setName(data.ip);}if(typeof data.origin==='string'&&data.origin&&!hostAddressEdited.current)updateHostAddress(new URL(data.origin).protocol===location.protocol?new URL(data.origin).host:data.origin);}catch{/* Show an error when Play tries to connect. */}
 },[ip,role]);
 const setHostAddress=useCallback((value:string)=>{hostAddressEdited.current=true;updateHostAddress(value);setError('');},[]);
 const start=useCallback(()=>{
  try{
   if(role==='host')invitationUrl(hostAddress,location.href,roomId);
   else{
    const target=guestInvitation(hostInput,location.href);
    if(new URL(target.url).origin!==location.origin){const url=new URL(target.url);url.searchParams.set('connectJoin','1');location.assign(url.href);return;}
   }
  }catch(cause){setError(cause instanceof Error?cause.message:'Invalid host invitation.');return;}
  if(desired.current){if(role==='host'&&desired.current.role==='host'){pendingStopRef.current=false;desired.current.playing=true;setPlaying(true);setError('');setNotice('');if(joinedRef.current)send({type:'play'});}return;}
  let id='';
  try{if(role==='guest')id=hostCode(hostInput);else id=sessionStorage.getItem('human-atlas-connect-room')??'';}catch(cause){setError(cause instanceof Error?cause.message:'Invalid host invitation.');return;}
  setError('');setNotice('');setRemoteSnapshot(null);poseRef.current=null;lastSnapshot.current='';setPlaying(role==='host');desired.current={role,roomId:id,playing:role==='host'};openSocketRef.current();
 },[role,hostInput,hostAddress,roomId,send]);
 const autoJoin=useRef(false);
 useEffect(()=>{if(autoJoin.current||role!=='guest'||new URLSearchParams(location.search).get('connectJoin')!=='1')return;autoJoin.current=true;const url=new URL(location.href);url.searchParams.delete('connectJoin');history.replaceState(null,'',url);start();},[role,start]);
 const stop=useCallback(()=>{
  if(role==='host'&&desired.current?.role==='host'){desired.current.playing=false;setPlaying(false);setGuests([]);lastSnapshot.current='';pendingStopRef.current=!joinedRef.current||socketRef.current?.readyState!==WebSocket.OPEN;if(!pendingStopRef.current)socketRef.current?.send(JSON.stringify({type:'stop'}));return;}
  if(socketRef.current?.readyState===WebSocket.OPEN)socketRef.current.send(JSON.stringify({type:'stop'}));finish('You left the session.');
 },[finish,role]);
 const setRole=useCallback((next:ConnectRole)=>{if(next===role)return;finish();updateRole(next);},[finish,role]);
 const renameSelf=useCallback((value:string)=>{
  const next=value.trim().slice(0,80);customNameRef.current=next;write(NAME_KEY,next);setName(next||ip);nameRef.current=next||ip;
  if(selfRef.current)send({type:'rename',id:selfRef.current.id,name:next});
 },[ip,send]);
 const renameGuest=useCallback((id:string,value:string)=>send({type:'rename',id,name:value}),[send]);
 const renameHost=useCallback((id:string,alias:string)=>setHosts(current=>{const next=current.map(item=>item.id===id?{...item,alias:alias.trim().slice(0,80)}:item);write(HISTORY_KEY,JSON.stringify(next));return next;}),[]);
 useEffect(()=>{
  disposedRef.current=false;
  return()=>{disposedRef.current=true;desired.current=null;joinedRef.current=false;for(const ref of [retryRef,snapshotTimer,poseTimer])if(ref.current)clearTimeout(ref.current);socketRef.current?.close();};
 },[]);
 const active=role==='host'?playing:status!=='idle',following=role==='guest'&&playing&&active;
 let invitation='';try{if(roomId)invitation=invitationUrl(hostAddress,location.href,roomId);}catch{/* Validate the address on Copy or Play. */}
 const copyInvitation=async()=>{try{if(!roomId)throw new Error('The host invitation is still being prepared.');await navigator.clipboard.writeText(invitationUrl(hostAddress,location.href,roomId));setError('');return true;}catch(cause){setError(cause instanceof Error?cause.message:'Could not copy the invitation.');return false;}};
 return {role,setRole,name,ip,renameSelf,hostAddress,setHostAddress,hostInput,setHostInput,status,error,notice,roomId,self,host,guests,hosts,renameGuest,renameHost,start,stop,prepare,active,following,invitation,copyInvitation,publishSnapshot,publishPose,remoteSnapshot,poseRef};
}
export type AtlasConnection=ReturnType<typeof useAtlasConnection>;
