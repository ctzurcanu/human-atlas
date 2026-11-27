import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';

const built=await build({entryPoints:['app/cardiac-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {cardiacRelations}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const stage='.local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1';
const read=path=>JSON.parse(fs.readFileSync(path));
const contracts=a=>{
 const parts=new Map(a.parts.map(p=>[p.id,p])),result=new Map();
 for(const p of a.parts)for(const c of p.cardiacCoupledAttachment?.interfaces??[]){
  if(!c.continuous||!['papillary-wall-insertion','leaflet-host-attachment','vessel-host-attachment','annular-leaflet-insertion','annular-myocardial-insertion'].includes(c.kind))continue;
  const q=parts.get(c.peer);
  if(q?.cardiacCoupledAttachment?.bankSha256!==p.cardiacCoupledAttachment.bankSha256||!q.cardiacCoupledAttachment.interfaces.some(r=>r.peer===p.id&&r.kind===c.kind&&r.continuous))continue;
  result.set([p.id,q.id].sort().join('|'),[p.id,q.id]);
 }
 return result;
};

test('heart-only fit invalidates fixed cuff continuity and preserves internal shared interfaces',()=>{
 const a=read(`${stage}/atlas.json`),g=cardiacRelations(a.parts),pending=a.finiteCardiacIntegration.pendingRootInterfaces;
 assert.equal(pending.length,16);
 const pendingKeys=new Set(pending.map(r=>[r.first,r.second].sort().join('|')));
 for(const r of pending)for(const [from,to] of [[r.first,r.second],[r.second,r.first]])assert.ok(!g.get(from)?.some(x=>x.target.id===to&&x.kind==='continuous'),`${from} / ${to}`);
 let downgraded=0,retained=0;
 for(const [key,[from,to]] of contracts(a)){
  const links=g.get(from)??[];
  if(pendingKeys.has(key)){assert.ok(links.some(x=>x.target.id===to&&x.kind==='adjacent'));downgraded++;}
  else{assert.ok(links.some(x=>x.target.id===to&&x.kind==='continuous'),key);retained++;}
 }
 assert.equal(downgraded,3);assert.ok(retained>20);
});

for(const sex of ['male','female'])test(`${sex}: common-pose coupled tissues retain continuity, including equivalent quaternion signs`,()=>{
 const a=read(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`);
 // Use a co-moving internal pair; this also remains valid after the staged fit.
 const pair=[...contracts(a).values()].find(([f,t])=>[f,t].every(id=>id.startsWith(`TA98HEART:${sex}:`)));
 assert.ok(pair);const [from,to]=pair,selected=a.parts.filter(p=>pair.includes(p.id));
 assert.ok(cardiacRelations(selected).get(from)?.some(r=>r.target.id===to&&r.kind==='continuous'));
 const flipped=structuredClone(selected),p=flipped.find(p=>p.id===to);
 p.geometryInstance={rotation:(p.geometryInstance?.rotation??[0,0,0,1]).map(v=>-v),scale:p.geometryInstance?.scale??1};
 assert.ok(cardiacRelations(flipped).get(from)?.some(r=>r.target.id===to&&r.kind==='continuous'));
 const displaced=structuredClone(selected),q=displaced.find(p=>p.id===to);q.geometryOrigin=[...(q.geometryOrigin??[0,0,0])];q.geometryOrigin[0]+=.001;
 assert.ok(!cardiacRelations(displaced).get(from)?.some(r=>r.target.id===to&&r.kind==='continuous'));
 assert.ok(cardiacRelations(displaced).get(from)?.some(r=>r.target.id===to&&r.kind==='adjacent'));
});
