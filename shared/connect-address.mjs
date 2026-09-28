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
export function invitationUrl(input,currentUrl,roomId){
 const current=new URL(currentUrl),url=invitationBase(input,current.origin);
 url.pathname=current.pathname;url.search='';url.hash='';url.searchParams.set('connect',roomId);
 if(current.searchParams.get('embed')==='1'){url.searchParams.set('embed','1');url.searchParams.set('ui','study,systems,details,model,open');}
 return url.href;
}
export function guestInvitation(input,currentUrl){
 const value=input.trim();if(!value)throw new Error('Paste the host’s invitation URL.');
 const url=new URL(value,currentUrl),code=url.searchParams.get('connect');
 if(!['http:','https:'].includes(url.protocol)||url.username||url.password||!code||!/^[\w-]{32}$/.test(code))throw new Error('Paste the host’s invitation URL.');
 if(localOnlyHost(url.hostname))throw new Error('Ask the host for an invitation using their IP or HTTPS address, not localhost.');
 return {url:url.href,code};
}
