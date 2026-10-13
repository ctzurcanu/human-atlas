import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs';import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/gallbladder-wall';export * from './app/bowel-wall-peeling';export * from './app/viewer-state';export * from './app/requested-selection-choice';export * from './app/guest-hierarchy';export * from './app/digestive-duct-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));const read=p=>JSON.parse(fs.readFileSync(p)),a=read('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'),ids=a.gallbladderWallPreview.partIds,byId=new Map(a.parts.map(p=>[p.id,p])),keys=['covering','subserosa','muscular','mucosa'],codes=['002','004','003','005'];
test('ordinary male atlas selects one regional/coat assembly and peels four complete intervals',()=>{
 assert.equal(ids.length,16);const selected=api.resolveSelection(a,['ZA:Gallbladder']);assert.deepEqual(new Set(selected),new Set(ids));
 const choice=api.requestedSelectionChoice(a,['ZA:Gallbladder'],selected);const steps=api.bowelWallPeelSteps(a,choice);assert.deepEqual(steps.map(s=>s.ids.length),[4,4,4,4]);
 for(const code of codes){const source='TA98PART:ZA:Gallbladder:A05.8.02.'+code;assert.ok(byId.get(source).suppressed);const selected=api.resolveSelection(a,[source]);assert.equal(selected.length,4);assert.ok(selected.every(id=>id.endsWith('A05.8.02.'+code)));const choice=api.requestedSelectionChoice(a,[source],selected);assert.equal(api.bowelWallPeelSteps(a,choice).length,4);}
 for(const key of keys){assert.ok(byId.get('TA98PART:ZA:Gallbladder:GALLWALL:'+key).suppressed);assert.equal(a.concepts.find(c=>c.id==='GALLWALL:male:'+key).elements.length,4);}
 assert.ok(ids.every(id=>!byId.get(id).suppressed));assert.ok(byId.get('ZA:Bile duct').suppressed);assert.ok(byId.get('FJ3080').pairedTissue&&byId.get('FJ3080').sectionSolid);assert.match(api.digestiveDuctDescription(byId.get('FJ3080')).text,/finite wall|open gallbladder neck/);
});
test('TA98 region and coat controls select their finite subsets, without claiming unallocated serosa',()=>{
 const nodes=new Map(),walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};api.resolveGuestHierarchy(a,api.parseGuestHierarchy(read('public/assets/ta98.json'))).forEach(walk);
 for(const code of codes){const p=nodes.get('A05.8.02.'+code).parts;assert.equal(p.length,4);assert.ok(p.every(p=>ids.includes(p.id)));}
 for(const code of ['007','008','009'])assert.equal(nodes.get('A05.8.02.'+code).directParts.length,4);
 assert.equal(nodes.get('GALLWALL:male:wall-assembly').children.length,4);assert.ok(!nodes.get('A05.8.02.006').parts.some(p=>ids.includes(p.id)));
});
test('every internal interface is reciprocal with identical triangle coordinates and reversed winding',()=>{
 const p0=byId.get(ids[0]);const chunk=a.chunks[p0.chunk];const raw=fs.readFileSync('.local-models/'+chunk.url.replace('/local-models/',''));
 const coords=new Float32Array(raw.buffer,raw.byteOffset+p0.positions,p0.vertexCount*3);const signature=(p,g,reverse=false)=>{const ix=new Uint32Array(raw.buffer,raw.byteOffset+p.indices,p.indexCount),rows=[];
  for(let j=g.start;j<g.start+g.count;j+=3){let pts=[0,1,2].map(k=>[...coords.slice(ix[j+k]*3,ix[j+k]*3+3)].join(','));if(reverse)pts=[pts[0],pts[2],pts[1]];rows.push([0,1,2].map(k=>pts.slice(k).concat(pts.slice(0,k)).join('|')).sort()[0]);}return rows.sort();};
 for(const id of ids){const p=byId.get(id);assert.equal(p.chunk,p0.chunk);assert.equal(p.positions,p0.positions);assert.equal(p.normals,p0.normals);assert.equal(p.uvs,undefined);
  for(const g of p.tissueDrawGroups.filter(g=>g.neighbor)){const target=byId.get(g.neighbor);assert.ok(ids.includes(target.id)||target.id==='FJ3080'&&api.gallbladderWallInterface(p,target));const reverse=target.tissueDrawGroups.find(g=>g.neighbor===id);assert.ok(reverse);assert.equal(reverse.count,g.count);assert.deepEqual(signature(p,g),signature(target,reverse,true));}
 }
});
