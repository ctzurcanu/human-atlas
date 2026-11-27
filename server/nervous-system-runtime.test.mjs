import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createRequire} from 'node:module';
import {build} from 'esbuild';

const bundle=await build({stdin:{contents:`
 export * from './app/shared-tissue-geometry';export * from './app/bowel-wall-bounds';
 export * from './app/hierarchy-visibility';export * from './app/hierarchy-navigation';
 export * from './app/preferred-source-representations';export * from './app/anatomy-hierarchy';
 export * from './app/guest-hierarchy';export * from './app/anatomical-item-presentation';
 export * from './app/viewer-state';export {BufferAttribute}from'three';
 export * from './app/organ-models';
 `,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const uiBundle=await build({stdin:{contents:`
 import React from 'react';import{renderToStaticMarkup}from'react-dom/server';
 import SystemTree,{buildSystemTree}from'./app/system-tree';import DepthTree from'./app/depth-tree';import GuestTree from'./app/guest-tree';
 export {buildSystemTree};export function render(type,props){return renderToStaticMarkup(React.createElement(type==='depth'?DepthTree:type==='guest'?GuestTree:SystemTree,props));}
 `,resolveDir:process.cwd(),sourcefile:'nervous-system-runtime-ui.tsx',loader:'tsx'},bundle:true,platform:'node',format:'cjs',write:false,define:{'import.meta.env.DEV':'false'}});
const ui={exports:{}};
new Function('require','module','exports',uiBundle.outputFiles[0].text)(createRequire(import.meta.url),ui,ui.exports);
const read=sex=>api.withHiddenReviewAlternatives(JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`)));
const atlas=read('female'),byId=new Map(atlas.parts.map(p=>[p.id,p])),marker=atlas.finiteSpinalIntegration;
const entries=api.hierarchyEntries(atlas),tree=ui.exports.buildSystemTree(entries,'systems',atlas),nervous=tree.find(b=>b.id==='nervous');
const guest=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
const base={focus:0,contextOpacity:1,skinOpacity:1,region:'all',peel:0,depthHidden:[],hidden:atlas.defaultHidden,labels:false,explode:0,visible:['skeletal'],selected:[],isolate:false,view:'right',rotate:false,reset:0,camera:[.5,.6,.7,.4,.3,.2,.5,.45,.6,.5,1,.5]};
const noop=()=>{};

test('the real nervous-system checkbox reveals its assembly and preserves the camera',t=>{
 assert.ok(nervous.parts.length>3400,'Exercise the whole system that caused the freeze');
 const off={...base,hidden:[...new Set([...base.hidden,...nervous.parts.map(p=>p.id)])]};
 const began=performance.now(),shown=api.toggleSystemParts(off,nervous.parts,atlas.parts.filter(p=>!p.suppressed));
 assert.ok(performance.now()-began<500);
 assert.deepEqual(shown.camera,base.camera);assert.equal(shown.view,base.view);assert.equal(shown.focus,base.focus);
 const selected=new Set(nervous.parts.map(p=>p.id));
 for(const p of nervous.parts)if(!p.spinalBoundaryReference)assert.ok(api.partVisible(p,shown,atlas.parts),p.name);
 for(const p of atlas.parts)if(!selected.has(p.id))assert.equal(api.partVisible(p,shown,atlas.parts),api.partVisible(p,off,atlas.parts),'The checkbox cannot move/hide the skeleton or enable unrelated anatomy');
 const hidden=api.toggleSystemParts(shown,nervous.parts,atlas.parts.filter(p=>!p.suppressed));
 assert.ok(nervous.parts.every(p=>!api.partVisible(p,hidden,atlas.parts)));assert.deepEqual(hidden.camera,base.camera);
 t.diagnostic(`${nervous.parts.length} physical pieces controlled by ${nervous.nodes.length} immediate named groups`);
});

test('all shared spinal cells allocate one exact selection ID and retain their original bank',t=>{
 const cells=[...marker.partIds,...marker.coveringPartIds,...marker.referencePartIds].map(id=>byId.get(id)),first=cells[0],raw=fs.readFileSync('.local-models'+atlas.chunks[first.chunk].url.replace('/local-models',''));
 const position=new api.BufferAttribute(new Float32Array(raw.buffer,raw.byteOffset+first.positions,first.vertexCount*3),3);
 let bytes=0,previousBytes=0,checked=0;const began=performance.now();
 for(const p of cells){
  const index=atlas.parts.indexOf(p),g=api.anatomyPartGeometry(p,index,atlas.parts.length),id=g.getAttribute('partIndex');
  assert.ok(g.isInstancedBufferGeometry);assert.equal(g.instanceCount,1);assert.ok(id.isInstancedBufferAttribute);assert.equal(id.count,1);assert.equal(id.getX(0),index);
  bytes+=id.array.byteLength;previousBytes+=p.vertexCount*2;
  g.setAttribute('position',position);g.setIndex(new api.BufferAttribute(new Uint32Array(raw.buffer,raw.byteOffset+p.indices,p.indexCount),1));
  const sphere=api.bowelWallBoundingSphere(p);assert.ok(sphere,'Every shared core cell must avoid a full position-bank bounds scan');
  g.boundingSphere=sphere;
  const display=api.anatomyGeometryView(g);assert.equal(display.instanceCount,1);assert.equal(display.getAttribute('position'),position);assert.equal(display.getIndex(),g.getIndex());assert.equal(display.getAttribute('partIndex'),id);
  const ids=g.getIndex().array,positions=position.array,centre=sphere.center,limit=sphere.radius*sphere.radius;
  for(let i=0;i<ids.length;i++){const j=ids[i]*3,dx=positions[j]-centre.x,dy=positions[j+1]-centre.y,dz=positions[j+2]-centre.z;if(dx*dx+dy*dy+dz*dz>limit+1e-16)assert.fail(`${p.id}: stored bounds exclude indexed vertex ${ids[i]}`);checked++;}
  display.dispose();g.dispose();
 }
 assert.equal(bytes,cells.length*2);assert.ok(previousBytes>8_000_000_000);assert.ok(bytes<10_000);
 t.diagnostic(`${previousBytes.toLocaleString()} → ${bytes.toLocaleString()} selection-ID bytes; ${checked.toLocaleString()} actual indexed vertices checked in ${Math.round(performance.now()-began)} ms`);
});

test('ordinary merged geometry retains per-vertex IDs and exact large-atlas IDs',()=>{
 const g=api.anatomyPartGeometry({vertexCount:8},5,100);assert.ok(!g.isInstancedBufferGeometry);assert.deepEqual([...g.getAttribute('partIndex').array],Array(8).fill(5));g.dispose();
 const large=api.anatomyPartGeometry({vertexCount:1_396_358,geometryOrigin:[0,0,0]},70000,80000);assert.equal(large.getAttribute('partIndex').getX(0),70000);assert.equal(large.getAttribute('partIndex').array.byteLength,4);large.dispose();
});

test('whole systems, regions, depth and TA98 badges sum their own mesh membership',()=>{
 for(const sex of ['female','male']){
  const a=api.withHiddenReviewAlternatives(JSON.parse(fs.readFileSync(`.local-models/organ-models/spinal-cord/${sex}-body.json`)));
  for(const mode of ['systems','regions','depth','guest:ta98']){
   const root=api.hierarchyNavigation(a,mode,mode==='guest:ta98'?guest:undefined),seen=new Set();
   const visit=c=>{if(seen.has(c))return;seen.add(c);assert.ok((c.children?.length??0)<1000,`${sex} ${mode}: ${c.name} has excessive siblings`);c.children?.forEach(visit);};visit(root);
   const props={atlas:a,state:{...base,hidden:a.defaultHidden},mode,onChoose:noop,setState:noop,hierarchy:guest};
   const html=ui.exports.render(mode==='guest:ta98'?'guest':mode==='depth'?'depth':'system',props);
   const counts=[...html.matchAll(/class="tree-count"[^>]*>([\d,]+)</g)].map(m=>Number(m[1].replaceAll(',','')));
   assert.ok(counts.length>0);
   assert.equal(counts[0],api.meshCount(root.elements.map(id=>a.parts.find(p=>p.id===id))));
  }
 }
});

test('generic member and relation lists group all organs by anatomical name',()=>{
 const cells=marker.partIds.map(id=>byId.get(id)),members=api.namedAnatomyMembers(cells,atlas.concepts),relations=api.namedAnatomyRelations(cells.map(target=>({kind:'contains',target})),atlas.concepts);
 assert.equal(members.length,96);assert.equal(relations.length,96);
 assert.deepEqual(new Set(members.flatMap(c=>c.elements)),new Set(marker.partIds));
 assert.deepEqual(new Set(relations.flatMap(c=>c.targetConcept.elements)),new Set(marker.partIds));
});
