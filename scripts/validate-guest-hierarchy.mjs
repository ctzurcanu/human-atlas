import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';

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
