// URL normalizes alternate IPv4 spellings, including 127.1 and integer forms.
export function localOnlyHost(hostname){
 const host=hostname.toLowerCase().replace(/^\[|\]$/g,'').replace(/\.$/,'');
 if(host==='localhost'||host.endsWith('.localhost')||host==='::'||host==='::1')return true;
 if(/^(127|0)\./.test(host))return true;
 const mapped=host.match(/^::ffff:([\da-f]+):([\da-f]+)$/);
 return !!mapped&&[0,127].includes(parseInt(mapped[1],16)>>>8);
}
export function invitationBase(input,currentOrigin){
 const value=input.trim();if(!value)throw new Error('Enter the host IP or HTTPS address.');
 const current=new URL(currentOrigin),explicit=/^[a-z][\w+.-]*:\/\//i.test(value);
 const url=new URL(explicit?value:`${current.protocol}//${value}`);
 if(!['http:','https:'].includes(url.protocol)||url.username||url.password)throw new Error('Enter the host IP or HTTPS address.');
 if(localOnlyHost(url.hostname))throw new Error('Use a reachable host address; localhost and loopback IPs cannot be shared.');
 // A bare IP refers to this viewer's port. A hostname can point at an HTTPS tunnel.
 if(!explicit&&!url.port&&(/^[\d.]+$/.test(url.hostname)||url.hostname.startsWith('[')))url.port=current.port;
 return url;
}
export function staticConnectHost(hostname){return hostname.toLowerCase().endsWith('.github.io');}
export function connectionOrigin(input,currentUrl){
 const current=new URL(currentUrl);
 const url=input.trim()?invitationBase(input,current.origin):current;
 if(staticConnectHost(url.hostname))throw new Error('GitHub Pages needs a live relay. Enter its HTTPS address in the Host field.');
 // Local development serves its own relay; the address is its public invitation.
 if(localOnlyHost(current.hostname))return current.origin;
 if(url.protocol!=='https:'){
  throw new Error('Use an HTTPS address for the live relay.');
 }
 return url.origin;
}
export function connectionSocketUrl(origin,roomId=''){
 const url=new URL('/atlas-connect',origin);url.protocol=url.protocol==='https:'?'wss:':'ws:';if(roomId)url.searchParams.set('room',roomId);return url.href;
}
export function invitationUrl(input,currentUrl,roomId,relayOrigin=''){
 const current=new URL(currentUrl),url=invitationBase(input,current.origin);
 // A static viewer stays on its own origin; its invitation carries the relay.
 const separate=relayOrigin&&relayOrigin!==current.origin&&!localOnlyHost(current.hostname);
 const viewer=separate?new URL(current):url;
 viewer.pathname=current.pathname;viewer.search='';viewer.hash='';viewer.searchParams.set('connect',roomId);
 if(separate)viewer.searchParams.set('relay',connectionOrigin(relayOrigin,currentUrl));
 if(current.searchParams.get('embed')==='1'){viewer.searchParams.set('embed','1');viewer.searchParams.set('ui','study,systems,details,model,open');}
 return viewer.href;
}
export function guestInvitation(input,currentUrl){
 const value=input.trim();if(!value)throw new Error('Paste the host’s invitation URL.');
 const url=new URL(value,currentUrl),code=url.searchParams.get('connect');
 if(!['http:','https:'].includes(url.protocol)||url.username||url.password||!code||!/^[\w-]{32}$/.test(code))throw new Error('Paste the host’s invitation URL.');
 if(localOnlyHost(url.hostname))throw new Error('Ask the host for an invitation using their IP or HTTPS address, not localhost.');
 if(url.searchParams.has('relay'))connectionOrigin(url.searchParams.get('relay')??'',url.href);
 return {url:url.href,code};
}
