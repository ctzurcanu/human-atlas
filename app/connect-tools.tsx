import {useEffect,useState} from 'react';
import {Check,Copy,Play} from 'lucide-react';
import {Button} from '@/components/ui/button';
import type {AtlasConnection} from './use-atlas-connection';
import type {ConnectRole} from './connect-state';

export function ConnectStopButton({role,onStop}:{role:ConnectRole;onStop:()=>void}){
 return <Button variant="ghost" className="connect-play-stop" aria-label="Stop Connect session" title={`Stop ${role==='host'?'Host':'Guest'} session`} onClick={onStop}><span className="connect-stop-symbol" aria-hidden="true">{role==='host'?'H':'G'}</span></Button>;
}

function NameField({value,label,onSave}:{value:string;label:string;onSave:(value:string)=>void}){
 const [draft,setDraft]=useState(value);useEffect(()=>setDraft(value),[value]);
 return <input aria-label={label} value={draft} maxLength={80} onChange={event=>setDraft(event.target.value)} onBlur={()=>{if(draft!==value)onSave(draft);}} onKeyDown={event=>{if(event.key==='Enter'){event.preventDefault();event.stopPropagation();event.currentTarget.blur();}else if(event.key==='Escape'){event.preventDefault();event.stopPropagation();setDraft(value);}}}/>;
}
export default function ConnectTools({connection:c}:{connection:AtlasConnection}){
 const [copied,setCopied]=useState(false);
 useEffect(()=>{void c.prepare();},[c.prepare]);
 useEffect(()=>{if(!copied)return;const timer=setTimeout(()=>setCopied(false),2000);return()=>clearTimeout(timer);},[copied]);
 if(c.active)return <div className="connect-tools is-running"><ConnectStopButton role={c.role} onStop={c.stop}/></div>;
 const previous=[...c.hosts].sort((a,b)=>Number(b.id===c.host?.id)-Number(a.id===c.host?.id));
 return <div className="connect-tools">
  <div className="connect-identity">
   <div className="connect-role-row">
    <select aria-label="Connect role" value={c.role} onChange={event=>c.setRole(event.target.value as 'host'|'guest')}><option value="host">Host</option><option value="guest">Guest</option></select>
    {c.role==='host'&&<Button variant="ghost" className="connect-copy" aria-label="Copy invitation URL" title={copied?'Invitation copied':'Copy invitation URL'} disabled={!c.roomId} onClick={async()=>setCopied(await c.copyInvitation())}>{copied?<Check size={16}/>:<Copy size={16}/>}</Button>}
   </div>
   <div className="connect-start-row">
    {c.role==='host'?<input aria-label="Host IP or HTTPS address" placeholder="Host IP or HTTPS address" title="Reachable host IP or HTTPS/ngrok address" value={c.hostAddress} onChange={event=>c.setHostAddress(event.target.value)} onKeyDown={event=>{if(event.key==='Enter'){event.preventDefault();c.start();}}}/>:<input aria-label="Host invitation URL" placeholder="Invitation URL" title="Paste the host invitation URL" value={c.hostInput} onChange={event=>c.setHostInput(event.target.value)} onKeyDown={event=>{if(event.key==='Enter'){event.preventDefault();c.start();}}}/>}
    <Button variant="ghost" className="connect-play-stop" aria-label={c.role==='host'?'Start Connect session':'Connect to host'} title={c.role==='host'?'Play: start hosting':'Play: join host'} onClick={c.start}><Play size={18}/></Button>
    {c.role==='host'&&<output className="connect-guest-count" aria-label="Connected guest count" title="Connected guests" aria-live="polite">{c.guests.filter(guest=>guest.online).length}</output>}
   </div>
  </div>
  <div className="connect-session-content">
   {c.role==='host'?<>
    <div className="connect-peer-list" aria-label="Connected guests">{c.guests.map(guest=><div className="connect-peer" key={guest.id} data-active={guest.online}><span className="connect-presence" aria-label={guest.online?'Connected':'Disconnected'}/><NameField value={guest.name} label={`Guest name ${guest.ip} (${guest.id})`} onSave={value=>c.renameGuest(guest.id,value)}/><span className="connect-ip">{guest.ip}</span></div>)}</div>
   </>:<>
    <div className="connect-peer-list" aria-label="Previous hosts">{previous.map(host=><div className="connect-peer" key={host.id}><span className="connect-presence" aria-label="Previous host"/><NameField value={host.alias||host.name} label={`Host name ${host.ip} (${host.id})`} onSave={value=>c.renameHost(host.id,value)}/><span className="connect-ip">{host.ip}</span><Button variant="ghost" className="connect-rejoin" aria-label={`Choose host ${host.alias||host.name}`} title="Choose this host" onClick={()=>c.setHostInput(host.url??new URL(`?connect=${host.roomId}`,location.href).href)}><Play size={14}/></Button></div>)}{!previous.length&&<span className="connect-empty">Previous hosts will appear here.</span>}</div>
   </>}
   {c.error&&<p className="connect-message" role="alert">{c.error}</p>}
  </div>
 </div>;
}
