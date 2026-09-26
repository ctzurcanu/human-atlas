import assert from 'node:assert/strict';
import {readFileSync,existsSync} from 'node:fs';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/anatomy-hierarchy.ts','app/depth-layers.ts','app/anatomy.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/hierarchy-taxonomy-validation'});
const load=async name=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(name)).text).toString('base64')}`);
const {REGION_ORDER,SKIN_DETAIL_NAMES,buildAnatomyNodes,depthPathFor,displaySystemFor,hierarchyEntries,regionPathFor,systemPathFor}=await load('anatomy-hierarchy.js');
const {DEPTH_LAYERS,depthLayerFor,isSkinPart}=await load('depth-layers.js');
const {SYSTEMS}=await load('anatomy.js');
const canonical=new Set([
 ...REGION_ORDER,...SYSTEMS.map(system=>system.name),'Skin','Body surface','Named skin regions','Left arm','Right arm','Left leg','Right leg',
 'Hand & wrist','Forearm & elbow','Shoulder','Upper arm','Foot','Thigh & hip','Lower leg',
 'Head','Neck','Brain','Skull','Face & senses','Scalp & epicranium','Back & spine','Pelvis','Abdomen','Chest','Other',
 'Mammary gland','Lacrimal glands','Subcutaneous tissue',
 'Cell boundary','Nucleus','Cytoplasm',...SKIN_DETAIL_NAMES,
]);
function leaves(nodes){return nodes.flatMap(node=>node.kind==='group'?leaves(node.nodes):node.kind==='bilateral'?node.entries.flatMap(entry=>entry.parts):node.entry.parts);}
function inspect(nodes){
 for(const node of nodes){
  if(node.kind!=='group')continue;
  assert(node.parts.length>0,`${node.name}: empty intermediate node`);
  assert(canonical.has(node.name),`${node.name}: noncanonical middle node`);
  assert.deepEqual(node.parts.map(part=>part.id).sort(),leaves(node.nodes).map(part=>part.id).sort(),`${node.name}: parent and leaf membership differ`);
  if(node.name==='Other')assert(node.nodes.every(child=>child.kind!=='group'),'Other should contain unmatched leaves directly');
  else assert(node.nodes.every(child=>child.kind==='group'),'Unassigned leaves must be under Other');
  inspect(node.nodes);
 }
}
function verify(atlas,name){
 const entries=hierarchyEntries(atlas),available=atlas.parts.filter(part=>!part.suppressed);
 const expected=available.map(part=>part.id).sort();
 assert.deepEqual(entries.flatMap(entry=>entry.parts).map(part=>part.id).sort(),expected,`${name}: entry coverage`);
 const modes={
  systems:()=>[...new Set(entries.map(entry=>displaySystemFor(entry.system)))].flatMap(system=>buildAnatomyNodes(entries.filter(entry=>displaySystemFor(entry.system)===system),entry=>systemPathFor(entry,atlas.scope))),
  regions:()=>[...new Set(entries.map(entry=>entry.region))].flatMap(region=>buildAnatomyNodes(entries.filter(entry=>entry.region===region),regionPathFor)),
  depth:()=>DEPTH_LAYERS.flatMap(layer=>{
   const members=entries.flatMap(entry=>{const parts=entry.parts.filter(part=>depthLayerFor(part)===layer.id);return parts.length?[{...entry,parts}]:[];});
   return members.length?buildAnatomyNodes(members,entry=>depthPathFor(entry,layer.name)):[];
  }),
 };
 for(const [mode,make] of Object.entries(modes)){
  const nodes=make();inspect(nodes);
  assert.deepEqual(leaves(nodes).map(part=>part.id).sort(),expected,`${name}: ${mode} has missing or duplicated items`);
 }
 if(name==='local-reference'){
  const skin=available.filter(part=>depthLayerFor(part)==='skin');
  assert.equal(skin.length,702,'All reference skin must be in the Skin depth layer');
  assert(skin.every(part=>displaySystemFor(part.system)==='integumentary'),'All reference skin must be under Body surface in Systems');
 }
 const expectedSkin={'male-detail':494,male:5,female:7,'local-male':2,'local-female':3,'local-reference':702};
 if(name in expectedSkin)assert.equal(available.filter(isSkinPart).length,expectedSkin[name],`${name}: skin classification changed`);
 for(const part of available.filter(part=>part.system==='integumentary'||part.system==='regions')){
  assert.equal(depthLayerFor(part)==='skin',isSkinPart(part),`${name}: ${part.name} has wrong skin depth`);
  const entry=entries.find(entry=>entry.parts.some(member=>member.id===part.id));
  assert(entry,`${name}: missing ${part.name}`);
  if(isSkinPart(part))assert(systemPathFor(entry,atlas.scope)[0]==='Skin',`${name}: ${part.name} is outside Skin`);
  else assert(systemPathFor(entry,atlas.scope)[0]!=='Skin',`${name}: ${part.name} is incorrectly under Skin`);
 }
 console.log(`${name}: ${available.length} pieces have one Systems, Regions, and Depth parent`);
}
const catalogues=[
 ['male-detail','public/models/atlas-male-complete.json'],['male','public/models/atlas.json'],
 ['female','public/models/atlas-hra-female.json'],['embryo','public/models/atlas-embryo.json'],
 ['cell','public/models/atlas-cell.json'],['local-male','.local-models/male.json'],
 ['local-female','.local-models/female.json'],['local-reference','.local-models/reference.json'],
];
const loaded=new Map();
for(const [name,file] of catalogues)if(existsSync(file)){
 const atlas=JSON.parse(readFileSync(file));loaded.set(name,atlas);verify(atlas,name);
}
if(loaded.has('male-detail')&&loaded.has('local-reference')){
 const detailed=loaded.get('male-detail'),reference=loaded.get('local-reference');
 const byName=new Map(hierarchyEntries(detailed).filter(entry=>systemPathFor(entry,detailed.scope)[0]==='Skin').map(entry=>[entry.name,entry]));
 let common=0;
 for(const entry of hierarchyEntries(reference).filter(entry=>systemPathFor(entry,reference.scope)[0]==='Skin')){
  const peer=byName.get(entry.name);if(!peer)continue;common++;
  assert.deepEqual(systemPathFor(entry,reference.scope),systemPathFor(peer,detailed.scope),`${entry.name}: system path differs between models`);
  assert.deepEqual(regionPathFor(entry),regionPathFor(peer),`${entry.name}: region path differs between models`);
 }
 assert(common>=250,`Only ${common} named skin structures share the canonical model paths`);
}
