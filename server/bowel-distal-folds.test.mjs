import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import crypto from 'node:crypto';
const sha=b=>crypto.createHash('sha256').update(b).digest('hex'),buffer=async p=>Uint8Array.from(await fs.readFile(p)).buffer;
for(const sex of ['female','male'])test(`${sex}: ileal distal direction replays same-source caecal evidence and has no distal physical folds`,async()=>{
 const dir='.local-models/ta98-review/small-intestine-layer-v1'+(sex==='male'?'-male':''),meta=JSON.parse(await fs.readFile(dir+'/manifest.json')),raw=await buffer(dir+'/wall.bin'),heightRaw=await buffer(dir+'/fold-height.bin'),organ=meta.organs.find(o=>o.id==='ileum'),n=organ.vertices;
 const chart=new Float32Array(raw,organ.chart,n*2),outer=new Float32Array(raw,organ.outer,n*3),height=new Float32Array(heightRaw,meta.foldRelief.organs.find(o=>o.id==='ileum').offset,n);
 const report=JSON.parse(await fs.readFile(`reports/ta98-${sex}/small-intestine-physical-folds.json`)),entry=report.evidence.find(o=>o.organ==='ileum'),proof=entry.regionalFoldOrientation;assert.equal(proof.terminalBoundaryDefined,false);
 const bytes=await fs.readFile(`.local-models/ta98-${sex}.json`);assert.equal(sha(bytes),proof.catalogueSha256);const atlas=JSON.parse(bytes),part=atlas.parts.find(p=>p.id===proof.caecumSource);assert.ok(part.id.startsWith(`LOCAL:${sex}:`));
 const url=atlas.chunks[part.chunk].url,path=url.startsWith('/local-models/')?'.local-models/'+url.slice('/local-models/'.length):'public'+url,chunk=await buffer(path),cv=new Float32Array(chunk,part.positions,part.vertexCount*3);
 for(const port of proof.ports)for(const [j,id] of port.controlIds.entries()){
  let best=Infinity;for(let i=0;i<part.vertexCount;i++){let distance=0;for(let d=0;d<3;d++)distance+=(outer[id*3+d]+meta.renderOrigin[d]-cv[i*3+d]-(part.geometryOrigin?.[d]??0))**2;best=Math.min(best,distance);}
  assert.ok(Math.abs(Math.sqrt(best)*1000-port.nearestCaecumVertexMm[j])<1e-6);
 }
 const s=Array.from({length:n},(_,i)=>chart[i*2]),lo=Math.min(...s),hi=Math.max(...s),reverse=Math.abs(proof.distalChartMeters-lo)<1e-5;let proximalRaised=0,distal=0;
 for(let i=0;i<n;i++){const q=reverse?1-(s[i]-lo)/(hi-lo):(s[i]-lo)/(hi-lo);if(q>=.65){assert.equal(height[i],0);distal++;}if(q<.35&&height[i]>0)proximalRaised++;}
 assert.ok(distal>n*.2);assert.ok(proximalRaised>0);assert.equal(entry.newInvertedPrisms,0);assert.equal(entry.outerFiveBoundariesExact,true);assert.equal(entry.allOpeningControlsExact,true);
});
