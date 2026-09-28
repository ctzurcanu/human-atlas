import {test} from 'node:test';
import assert from 'node:assert/strict';
import {connectionOrigin,connectionSocketUrl,guestInvitation,invitationBase,invitationUrl} from '../shared/connect-address.mjs';
import {connectionIdentity,hostAddresses} from './connect-address.mjs';
import {requestIp} from './connect-relay.mjs';
import {allowedConnectOrigin,connectIdentityCors} from './connect-origin.mjs';

const room='abcdefghijklmnopqrstuvwxyz123456';
test('invitations reject loopback aliases and unspecified addresses',()=>{
 for(const address of ['localhost','LOCALHOST.','viewer.localhost','127.0.0.1','127.0.01','127.1','127.44.55.66','2130706433','0x7f000001','0.0.0.0','[::]','[::1]','[0:0:0:0:0:0:0:1]','[::ffff:127.0.0.1]']){
  assert.throws(()=>invitationBase(address,'http://localhost:3016'),/loopback/,address);
  assert.throws(()=>guestInvitation(`https://${address}/?connect=${room}`,'https://atlas.example'),/not localhost/,address);
 }
 for(const address of ['ftp://atlas.example','https://user:secret@atlas.example'])assert.throws(()=>invitationBase(address,'https://atlas.example'));
});
test('reachable IPs and HTTPS tunnel invitations preserve the viewer path and embedding',()=>{
 assert.equal(invitationUrl('192.168.1.25','http://localhost:3016/?model=male-detail',room),`http://192.168.1.25:3016/?connect=${room}`);
 assert.equal(invitationUrl('https://atlas.ngrok-free.app','http://localhost:3016/viewer/?embed=1&model=male-detail',room),`https://atlas.ngrok-free.app/viewer/?connect=${room}&embed=1&ui=study%2Csystems%2Cdetails%2Cmodel%2Copen`);
 assert.equal(invitationBase('https://192.168.1.25:8443','http://localhost:3016').origin,'https://192.168.1.25:8443');
});
test('Guest accepts another viewer origin and validates the session code',()=>{
 const url=`https://atlas.ngrok-free.app/?connect=${room}`;
 assert.deepEqual(guestInvitation(url,'http://localhost:3016'),{url,code:room});
 assert.throws(()=>guestInvitation('https://atlas.example/?connect=bad','http://localhost:3016'),/invitation URL/);
});
test('host discovery excludes loopback interfaces and honors the HTTPS proxy address',async()=>{
 assert.deepEqual(hostAddresses({lo:[{family:'IPv4',internal:true,address:'127.0.0.1'}],en0:[{family:'IPv4',internal:false,address:'192.168.1.25'},{family:'IPv6',internal:false,address:'fe80::1'}]}),['192.168.1.25']);
 const identity=await connectionIdentity({headers:{host:'localhost:3016','x-forwarded-host':'atlas.ngrok-free.app','x-forwarded-proto':'https'},socket:{}});
 assert.equal(identity.origin,'https://atlas.ngrok-free.app');
 assert.notEqual(requestIp({headers:{},socket:{remoteAddress:'::ffff:127.0.0.1'}}),'127.0.0.1');
});
test('GitHub Pages uses an HTTPS relay and keeps invitations on the static viewer',()=>{
 const viewer='https://ctzurcanu.github.io/human-atlas/?model=male-detail&hide=long-state',relay='https://atlas.example';
 assert.throws(()=>connectionOrigin('',viewer),/GitHub Pages needs a live relay/);
 assert.throws(()=>connectionOrigin('ctzurcanu.github.io',viewer),/GitHub Pages needs a live relay/);
 assert.throws(()=>connectionOrigin('http://atlas.example',viewer),/HTTPS address/);
 assert.throws(()=>connectionOrigin('https://127.0.0.1',viewer),/loopback/);
 assert.equal(connectionOrigin(relay+'/ignored-path',viewer),relay);
 assert.equal(connectionSocketUrl(connectionOrigin(relay,viewer)),'wss://atlas.example/atlas-connect');
 assert.equal(connectionSocketUrl(relay,room),`wss://atlas.example/atlas-connect?room=${room}`);
 const invitation=invitationUrl(relay,viewer,room,relay),url=new URL(invitation);
 assert.equal(url.origin,'https://ctzurcanu.github.io');assert.equal(url.pathname,'/human-atlas/');
 assert.equal(url.searchParams.get('connect'),room);assert.equal(url.searchParams.get('relay'),relay);
 assert.equal(url.searchParams.has('hide'),false);
 assert.deepEqual(guestInvitation(invitation,viewer),{url:invitation,code:room});
 assert.throws(()=>guestInvitation(invitation.replace(encodeURIComponent(relay),encodeURIComponent('http://atlas.example')),viewer),/HTTPS/);
 const embed=new URL(invitationUrl(relay,viewer+'&embed=1',room,relay));assert.equal(embed.searchParams.get('embed'),'1');assert.ok(embed.searchParams.has('ui'));
});
test('local development and same-origin HTTPS retain their existing relay and invitations',()=>{
 assert.equal(connectionOrigin('https://atlas.ngrok-free.app','http://localhost:3016'),'http://localhost:3016');
 assert.equal(connectionOrigin('192.168.1.25','http://localhost:3016'),'http://localhost:3016');
 assert.equal(connectionSocketUrl('http://localhost:3016'),'ws://localhost:3016/atlas-connect');
 assert.equal(connectionOrigin('atlas.ngrok-free.app','https://atlas.ngrok-free.app'),'https://atlas.ngrok-free.app');
 assert.equal(invitationUrl('https://atlas.ngrok-free.app','http://localhost:3016/',room,'http://localhost:3016'),`https://atlas.ngrok-free.app/?connect=${room}`);
});
test('identity CORS and WebSockets allow our GitHub viewer only over a secure relay',()=>{
 const request={headers:{origin:'https://ctzurcanu.github.io',host:'relay.example','x-forwarded-proto':'https'},socket:{}};
 const headers=new Map(),response={setHeader:(key,value)=>headers.set(key,value)};
 assert.equal(allowedConnectOrigin(request),true);assert.equal(connectIdentityCors(request,response),true);
 assert.equal(headers.get('Access-Control-Allow-Origin'),'https://ctzurcanu.github.io');assert.equal(headers.get('Vary'),'Origin');
 for(const origin of ['https://unrelated.example','https://ctzurcanu.github.io.attacker.example','http://ctzurcanu.github.io','https://ctzurcanu.github.io/path','null']){
  const rejected={...request,headers:{...request.headers,origin}};
  assert.equal(allowedConnectOrigin(rejected),false,origin);assert.equal(connectIdentityCors(rejected,response),false,origin);
 }
 assert.equal(allowedConnectOrigin({...request,headers:{origin:request.headers.origin,host:'relay.example'}}),false);
});
