import {networkInterfaces} from 'node:os';
import {localOnlyHost} from '../shared/connect-address.mjs';

export function hostAddresses(interfaces=networkInterfaces()){
 return [...new Set(Object.values(interfaces).flatMap(entries=>(entries??[]).filter(entry=>entry.family==='IPv4'&&!entry.internal&&!localOnlyHost(entry.address)).map(entry=>entry.address)))];
}
export async function connectionIdentity(request){
 const addresses=hostAddresses();
 const host=String(request.headers['x-forwarded-host']??request.headers.host??'').split(',')[0].trim();
 const protocol=request.socket.encrypted||String(request.headers['x-forwarded-proto']??'').split(',')[0].trim()==='https'?'https:':'http:';
 let origin='';
 try{
  const viewer=new URL(`${protocol}//${host}`);
  if(!localOnlyHost(viewer.hostname))origin=viewer.origin;
  else{
   // Prefer the running tunnel for this port when ngrok's local agent is available.
   try{
    const response=await fetch('http://127.0.0.1:4040/api/tunnels',{signal:AbortSignal.timeout(750)});
    const data=await response.json();
    for(const tunnel of data.tunnels??[]){
     const target=new URL(String(tunnel.config?.addr).includes('://')?tunnel.config.addr:`http://${tunnel.config?.addr}`),publicUrl=new URL(tunnel.public_url);
     if(target.port===viewer.port&&localOnlyHost(target.hostname)&&publicUrl.protocol==='https:'&&!localOnlyHost(publicUrl.hostname)){origin=publicUrl.origin;break;}
    }
   }catch{/* A LAN address can be entered when there is no local tunnel agent. */}
   if(!origin&&addresses[0])origin=`${protocol}//${addresses[0]}${viewer.port?`:${viewer.port}`:''}`;
  }
 }catch{/* Let the user enter the reachable host address. */}
 return {addresses,origin};
}
