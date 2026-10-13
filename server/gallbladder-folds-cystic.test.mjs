import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/gallbladder-wall';export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const read=p=>JSON.parse(fs.readFileSync(p));
for(const sex of ['male','female'])test(`${sex} folds and neck use the installed finite tissue, not overlapping feature solids`,()=>{
 const a=read(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`),byId=new Map(a.parts.map(p=>[p.id,p]));
 const fold=a.concepts.find(c=>c.id===`GALLWALL:${sex}:folds`);
 assert.equal(fold.ta98Term,'A05.8.02.010');assert.equal(fold.partOf,`GALLWALL:${sex}:mucosa`);
 assert.deepEqual(fold.elements,a.gallbladderWallPreview.partIds.filter(id=>byId.get(id).tissueWallLayer==='mucosa'));
 const tree=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(read('public/assets/ta98.json'))),nodes=new Map();
 const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};tree.forEach(walk);
 assert.deepEqual(new Set(nodes.get('A05.8.02.010').directParts.map(p=>p.id)),new Set(fold.elements));
 const femaleCystic=sex==='female'&&a.biliaryJunctionIntegration?.continuousConfluence;
 const duct=byId.get(sex==='male'?'FJ3080':femaleCystic?'LOCAL:female:cystic_duct:authored':'LOCAL:female:bile_duct:1031'),contacts=a.gallbladderWallPreview.partIds.map(id=>byId.get(id)).filter(p=>p.tissueDrawGroups.some(g=>g.neighbor===duct.id));
 assert.equal(contacts.length,4);
 const raw=fs.readFileSync('.local-models/'+a.chunks[duct.chunk].url.replace('/local-models/',''));
 const signature=(p,target,reverse=false)=>{
  const coords=new Float32Array(raw.buffer,raw.byteOffset+p.positions,p.vertexCount*3),ix=new Uint32Array(raw.buffer,raw.byteOffset+p.indices,p.indexCount),rows=[];
  for(const g of p.tissueDrawGroups.filter(g=>g.neighbor===target))for(let i=g.start;i<g.start+g.count;i+=3){let pts=[0,1,2].map(k=>[...coords.slice(ix[i+k]*3,ix[i+k]*3+3)].join(','));if(reverse)pts=[pts[0],pts[2],pts[1]];rows.push([0,1,2].map(k=>pts.slice(k).concat(pts.slice(0,k)).join('|')).sort()[0]);}
  return rows.sort();
 };
 const graph=api.gallbladderWallRelations(a.parts);
 for(const p of contacts){assert.ok(api.gallbladderWallInterface(p,duct));assert.equal(p.chunk,duct.chunk);assert.deepEqual(signature(p,duct.id),signature(duct,p.id,true));assert.ok(graph.get(p.id).some(r=>r.kind==='continuous'&&r.target.id===duct.id));assert.ok(graph.get(duct.id).some(r=>r.kind==='continuous'&&r.target.id===p.id));assert.equal(api.gallbladderWallInterface({...p,sourceOffset:[.01,0,0]},duct),false);}
 const report=read(`reports/ta98-${sex}/${femaleCystic?'cystic-local-connection-installed':'gallbladder-folds-cystic-installed'}.json`);assert.ok(report.installed&&report.allUnrelatedPartsExact&&report.defaultHiddenExact);
 assert.ok(report.finiteCells.every(p=>p.components===1&&p.openEdges===0&&p.selfIntersections===0));
 assert.equal(report.construction.folds.newUVBytes,0);assert.equal(report.construction.folds.newImageBytes,0);
 assert.ok(report.construction.folds.peakReliefMetres>0);
 if(femaleCystic){assert.equal(report.construction.folds.seamUnchanged,false);assert.ok(report.construction.folds.newNeckSharedFieldsExact&&report.construction.folds.originalTerminalSupportReplaced);}
 else assert.ok(report.construction.folds.seamUnchanged);
});
