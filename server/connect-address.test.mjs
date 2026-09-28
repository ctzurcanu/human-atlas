import {test} from 'node:test';
import assert from 'node:assert/strict';
import {guestInvitation,invitationBase,invitationUrl} from '../shared/connect-address.mjs';
import {connectionIdentity,hostAddresses} from './connect-address.mjs';
import {requestIp} from './connect-relay.mjs';

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
