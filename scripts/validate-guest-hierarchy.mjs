import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';

async function load(entryPoint){
 const result=await build({entryPoints:[entryPoint],bundle:true,platform:'node',format:'esm',write:false});
 return import(`data:text/javascript;base64,${Buffer.from(result.outputFiles[0].contents).toString('base64')}`);
}

const {parseGuestHierarchy,resolveGuestHierarchy}=await load('app/guest-hierarchy.ts');
const {hierarchyEntries}=await load('app/anatomy-hierarchy.ts');
const {createHierarchicalExplosionLayout}=await load('app/hierarchical-explosion.ts');
const raw=JSON.parse(await readFile(new URL('../public/assets/chakras.json',import.meta.url)));
const hierarchy=parseGuestHierarchy(raw);
assert.equal(hierarchy.id,'chakras');
assert.equal(hierarchy.nodes.length,8);
assert.deepEqual(hierarchy.nodes.map(node=>node.name.split(' · ')[0]),['Ajna','Vishuddha','Anahata','Manipura','Manipura','Swaddhishtana','Muladhara','Other']);

const atlas=JSON.parse(await readFile(new URL('../public/models/atlas-male-complete.json',import.meta.url)));
const nodes=resolveGuestHierarchy(atlas,hierarchy);
for(const name of ['Retinas','Pineal','Cerebellum','Thyroid','Heart','Stomach','Liver','Pancreas','Kidneys','Testes','Colon','Bones']){
 const find=(items)=>items.flatMap(node=>[node,...find(node.children)]);
 const item=find(nodes).find(node=>node.name===name);
 assert(item?.directParts.length>0,`${name} should resolve to model anatomy`);
}
const ajna=nodes.find(node=>node.name==='Ajna');
assert(ajna?.children.find(node=>node.name==='Glands')?.children.find(node=>node.name==='Pituitary'),'Unmodeled source entries should remain in the tree');
assert(nodes.every(node=>node.parts.length>0),'Every branch must have modeled anatomy');
const mapped=new Set(nodes.slice(0,-1).flatMap(node=>node.parts.map(part=>part.id)));
assert.deepEqual(new Set(nodes.at(-1).parts.map(part=>part.id)),new Set(atlas.parts.filter(part=>!mapped.has(part.id)).map(part=>part.id)),'Other must contain exactly the unmapped parts');
const visible=new Set(atlas.parts.map(part=>part.id));
const layout=createHierarchicalExplosionLayout(atlas,hierarchyEntries(atlas),visible,'guest',1.7,nodes);
assert(layout.steps>=4,'Guest Explode should follow several hierarchy levels');
for(const stage of layout.stages){
 const members=stage.clusters.flat();
 assert.equal(members.length,atlas.parts.length,'Every mesh should have a rotation cluster');
 assert.equal(new Set(members).size,atlas.parts.length,'Guest groups must not duplicate meshes');
 assert(stage.positions.every(Number.isFinite),'Guest positions must be finite');
}
assert(layout.stages.at(-1).clusters.every(members=>members.length===1),'Guest Explode must end at individual meshes');
console.log(`Chakras: ${nodes.length} root branches, ${layout.steps} Explode levels, ${atlas.parts.length} meshes assigned once`);

const input=await loadChakraInputs();
let nerves=0,innervated=0,pieces=0;
for(const {file,atlas} of input.catalogues){
 const roots=resolveGuestHierarchy(atlas,hierarchy).slice(0,-1);
 const sets=roots.map(root=>new Set(root.parts.map(part=>part.id)));
 const entries=input.hierarchyEntries(atlas),byConcept=new Map(entries.map(entry=>[entry.id,entry]));
 const canonicalParts=atlas.parts.filter(part=>!part.suppressed).map(part=>({...part,name:byConcept.get(part.conceptId)?.name??part.name}));
 const mapped=new Set(roots.flatMap(root=>root.parts.map(part=>part.id)));
 const other=resolveGuestHierarchy(atlas,hierarchy).at(-1);
 assert.deepEqual(new Set(other.parts.map(part=>part.id)),new Set(canonicalParts.filter(part=>!mapped.has(part.id)).map(part=>part.id)),`${file}: Other must be the exact complement`);
 for(const entry of entries){
  const path=input.metadata.byConceptPath[entry.id]??[];
  const spinal=path.some(code=>/^A14\.2\.0[2-7]\./.test(code))||path.includes('A14.2.00.027')||/spinal (?:nerve|ganglion)|cauda equina|musculocutan(?:eous|eus)|palmar digital|plantar digital|median nerve|ulnar nerve|radial nerve|sciatic nerve|tibial nerve|(?:fibular|peroneal) nerve|femoral nerve|obturator nerve|saphenous nerve|phrenic nerve/i.test(entry.name);
  const vagus=path.includes('A14.2.01.153')||/vagus|vagal|laryngeal nerve/i.test(entry.name);
  if((spinal||vagus)&&input.majorSystemFor(entry)==='nervous'){
   for(const part of entry.parts)assert(mapped.has(part.id),`${file}: ${entry.name} must belong to a chakra`);
   nerves++;
  }
  const expected={circulatory:2,respiratory:2,digestive:3,muscular:3,skeletal:6,lymphatic:5}[input.majorSystemFor(entry)];
  if(expected!==undefined)for(const part of entry.parts)assert(sets[expected].has(part.id),`${file}: ${entry.name} is missing its requested system chakra`);
  if(entry.parts.some(input.isSkinPart))for(const part of entry.parts)assert(sets[2].has(part.id),`${file}: skin must belong to Anahata`);
  if(path.some(code=>code.startsWith('A15.3.03.'))||/vestibulocochlear|cochlear nerve|vestibular nerve/i.test(entry.name))for(const part of entry.parts)assert(sets[1].has(part.id),`${file}: inner ear and hearing must belong to Vishuddha`);
  if(input.majorSystemFor(entry)!=='nervous'||!entry.parts.some(part=>mapped.has(part.id)))continue;
  const relations=entry.parts.flatMap(part=>input.anatomicalRelations({...part,name:entry.name},canonicalParts));
  const sourceRoots=sets.filter(set=>entry.parts.some(part=>set.has(part.id)));
  for(const relation of relations.filter(link=>link.kind==='innervates')){
   for(const set of sourceRoots)assert(set.has(relation.target.id),`${file}: ${entry.name} and innervated ${relation.target.name} must share every assigned chakra`);
   innervated++;
  }
 }
 pieces+=canonicalParts.length;
 console.log(`${file}: nerve coverage, organ membership and system assignments pass`);
}

// Every one of the 31 named spinal levels resolves even if absent in an import.
const levelNames=[...Array.from({length:8},(_,i)=>`Cervical spinal nerve c${i+1}`),...Array.from({length:12},(_,i)=>`Thoracic spinal nerve t${i+1}`),...Array.from({length:5},(_,i)=>`Lumbar spinal nerve l${i+1}`),...Array.from({length:5},(_,i)=>`Sacral spinal nerve s${i+1}`),'Coccygeal spinal nerve co1'];
const fixtures=levelNames.flatMap(name=>['left','right'].map(side=>({name:`${name} (${side})`,id:`fixture:${name}:${side}`,conceptId:`fixture:${name}:${side}`,system:'nervous',groups:[]})));
const fixtureAtlas={...atlas,parts:fixtures,concepts:fixtures.map(part=>({id:part.conceptId,name:part.name,elements:[part.id]}))};
const fixtureMapped=new Set(resolveGuestHierarchy(fixtureAtlas,hierarchy).slice(0,-1).flatMap(root=>root.parts.map(part=>part.id)));
assert.equal(fixtureMapped.size,62,'All 31 spinal levels must be mapped bilaterally');
assert.throws(()=>parseGuestHierarchy({...raw,nodes:[{name:'Invalid',unmapped:true,matchBases:['Heart']}]}));
assert.throws(()=>parseGuestHierarchy({...raw,nodes:[{name:'Invalid',matchBases:Array(2001).fill('Heart')}]}));
console.log(`${pieces} available catalogue pieces checked; ${nerves} spinal/vagus concepts covered; ${innervated} innervation links share a chakra; all 31 spinal levels covered bilaterally.`);
