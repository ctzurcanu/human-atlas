import assert from 'node:assert/strict';
import {readFile,readdir} from 'node:fs/promises';
import {build} from 'esbuild';
const result=await build({stdin:{contents:"export * from './app/guest-hierarchy';export * from './app/hierarchy-navigation';export * from './app/anatomy-hierarchy';export * from './app/hierarchical-explosion';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import(`data:text/javascript;base64,${Buffer.from(result.outputFiles[0].contents).toString('base64')}`);
const raw=await Promise.all(['genes','cell-types'].map(async id=>JSON.parse(await readFile(`public/assets/${id}.json`))));
const hierarchies=raw.map(api.parseGuestHierarchy),provenance=JSON.parse(await readFile('public/assets/biology-provenance.json'));
const [genes,cells]=hierarchies,byGene=new Map(genes.nodes.map(node=>[node.id,node])),byCell=new Map(cells.nodes.map(node=>[node.id,node]));
assert.equal(genes.nodes.filter(node=>/^HGNC:\d+$/.test(node.id)).length,provenance.statistics.genes);
assert.equal(genes.nodes.filter(node=>node.id.startsWith('HGNC-GROUP:')).length,provenance.statistics.families);
assert.equal(cells.nodes.filter(node=>node.id.startsWith('CL:')&&!node.name.endsWith(' · ASCT+B')).length,provenance.statistics.humanCellTypes);
for(const h of hierarchies){
 const byId=new Map(h.nodes.map(node=>[node.id,node]));
 for(const node of h.nodes)for(const id of node.childIds??[])assert(byId.has(id));
 for(const node of h.nodes)for(const link of node.links??[])assert((link.hierarchy==='genes'?byGene:byCell).has(link.id),`Missing biology target ${link.id}`);
 const depth=new Map();const maxDepth=id=>{if(depth.has(id))return depth.get(id);const n=1+Math.max(0,...(byId.get(id).childIds??[]).map(maxDepth));depth.set(id,n);return n;};
 const max=Math.max(...h.roots.map(maxDepth));assert(max>=6,`${h.name} must preserve deep source branches`);console.log(`${h.name}: ${h.nodes.length} unique nodes, ${max} levels, source references intact`);
}
// Shared source nodes are resolved once, regardless of their number of parents.
const fixture={schema:'human-atlas-hierarchy/v2',id:'fixture',name:'Fixture',roots:['a','b'],nodes:[{id:'a',name:'A',childIds:['c']},{id:'b',name:'B',childIds:['c']},{id:'c',name:'C',matches:['Heart']}]};
assert.throws(()=>api.parseGuestHierarchy({...fixture,nodes:[...fixture.nodes.slice(0,2),{id:'c',name:'C',childIds:['a']}]}),/cycle/);
assert.throws(()=>api.parseGuestHierarchy({...fixture,roots:['missing']}),/Unknown/);
assert.throws(()=>api.parseGuestHierarchy({...fixture,nodes:[...fixture.nodes,{id:'unreachable',name:'Unreachable'}]}),/unreachable/);
const files=provenance.statistics.catalogues;
for(const file of files){
 const atlas=JSON.parse(await readFile(file)),available=atlas.parts.filter(part=>!part.suppressed),ids=new Set(available.map(part=>part.id));
 const shared=api.resolveGuestHierarchy(atlas,api.parseGuestHierarchy(fixture));assert.equal(shared[0].children[0],shared[1].children[0]);
 for(const h of hierarchies)for(const view of h.views??[{id:'',roots:h.roots}]){
  const active=view.id?{...h,roots:view.roots}:h,start=performance.now(),nodes=api.resolveGuestHierarchy(atlas,active);assert.equal(api.resolveGuestHierarchy(atlas,active),nodes,'Resolution cache must be reused');
  const seen=new Map(),visit=node=>{if(seen.has(node.id))return;seen.set(node.id,node);assert.equal(new Set(node.parts.map(p=>p.id)).size,node.parts.length);node.parts.forEach(p=>assert(ids.has(p.id)));node.children.forEach(visit);};nodes.forEach(visit);
  const nav=api.hierarchyNavigation(atlas,`guest:${h.id}`,active),first=nodes.find(node=>node.parts.length);if(first){const choice=nav.children.find(node=>node.id.endsWith(':'+first.id));assert(choice);assert(api.hierarchyAncestors(nav,choice).some(node=>node===nav));}
  const layout=api.createHierarchicalExplosionLayout(atlas,api.hierarchyEntries(atlas),ids,'guest',1.7,nodes);
  for(const stage of layout.stages){const members=stage.clusters.flat();assert.equal(members.length,available.length);assert.equal(new Set(members).size,available.length);assert(stage.positions.every(Number.isFinite));}
  console.log(`${file} · ${h.name} ${view.id}: ${Math.round(performance.now()-start)} ms, ${layout.steps} Explode levels, ${available.length} meshes assigned once`);
 }
}
const atlas=JSON.parse(await readFile('public/models/atlas-male-complete.json'));
for(const [h,id] of [[genes,'HGNC:399'],[cells,'CL:0000182']]){
 const seen=new Map(),visit=n=>{if(seen.has(n.id))return;seen.set(n.id,n);n.children.forEach(visit);};api.resolveGuestHierarchy(atlas,h).forEach(visit);
 const node=seen.get(id);assert(node.directParts.some(p=>/Liver/i.test(p.name)));assert(!node.directParts.some(p=>/Liver · Ligament/.test(p.name)),'Parenchymal mapping must not include bundled ligament materials');
}
const state={guestView:'lineage',guestQuery:'CL:0000182',guestExtensions:['HGNC:399'],selected:[],visible:[],explode:0,isolate:false,view:'front',rotate:false,reset:0};
const restored=api.readViewUrl(new URL(api.viewUrl('http://localhost:3016','male-detail',state)).search,atlas,state);assert.equal(restored.guestView,state.guestView);assert.equal(restored.guestQuery,state.guestQuery);assert.deepEqual(restored.guestExtensions,state.guestExtensions);
const extensions=genes.nodes.filter(node=>node.extension);
let sampleChunk;
for(const node of extensions){const chunk=api.parseGuestHierarchy(JSON.parse(await readFile('public'+node.extension)));assert.equal(chunk.roots[0],node.id);assert(chunk.nodes.length>1);assert(!chunk.nodes[0].extension);for(const n of chunk.nodes.slice(1))assert(!byGene.has(n.id));if(node.id==='HGNC:399')sampleChunk=chunk;}
if(sampleChunk){const added=new Map(sampleChunk.nodes.map(n=>[n.id,n])),existing=new Set(genes.nodes.map(n=>n.id));const merged=api.parseGuestHierarchy({...genes,nodes:genes.nodes.map(n=>added.get(n.id)??n).concat(sampleChunk.nodes.filter(n=>!existing.has(n.id)))});assert(merged.nodes.length>genes.nodes.length);}

if(extensions.length)console.log(`${extensions.length} transcript/protein chunks are valid and merge without duplicate or dangling IDs.`);
