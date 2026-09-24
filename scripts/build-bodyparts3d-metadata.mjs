/** Join official BodyParts3D IS-A and PART-OF tables into converter inputs.
 * The element-to-system and default-concept choices are Human Atlas curation,
 * kept separately in scripts/data/bodyparts3d-elements.json. No intermediary
 * atlas geometry, labels, or catalogue is read during a rebuild.
 * Usage: node scripts/build-bodyparts3d-metadata.mjs TABLE_DIRECTORY OUTPUT_DIRECTORY
 */
import fs from 'node:fs/promises';
import path from 'node:path';
import assert from 'node:assert/strict';
import {fileURLToPath} from 'node:url';

const tables=path.resolve(process.argv[2]),out=path.resolve(process.argv[3]);
await fs.mkdir(out,{recursive:true});
const curated=JSON.parse(await fs.readFile(new URL('./data/bodyparts3d-elements.json',import.meta.url)));
const names=new Map(),members=new Map();
for(const tree of ['isa','partof']){
 const parts=await fs.readFile(path.join(tables,`${tree}_parts_list_e.txt`),'utf8');
 for(const line of parts.trim().split(/\r?\n/).slice(1)){
  const [id,representation,name]=line.split('\t');
  if(id?.startsWith('FMA')&&name&&!names.has(id))names.set(id,name);
 }
 const elements=await fs.readFile(path.join(tables,`${tree}_element_parts.txt`),'utf8');
 for(const line of elements.trim().split(/\r?\n/).slice(1)){
  const [conceptId,name,elementId]=line.split('\t');
  if(!conceptId?.startsWith('FMA')||!elementId?.startsWith('FJ'))continue;
  if(!names.has(conceptId))names.set(conceptId,name);
  const set=members.get(conceptId)??new Set();set.add(elementId);members.set(conceptId,set);
 }
}
const ids=new Set(Object.keys(curated));
assert.equal(ids.size,2234,'Curated BodyParts3D element count changed');
for(const id of ids)assert.ok([...members.values()].some(set=>set.has(id)),`No publisher metadata for ${id}`);
const concepts=[...members].map(([id,set])=>({id,name:names.get(id),elements:[...set].filter(id=>ids.has(id)).sort()})).filter(c=>c.elements.length).sort((a,b)=>a.id.localeCompare(b.id));
assert.equal(concepts.length,3432,'Publisher concept count changed');
const elements=[...ids].sort().map(id=>({id,name:names.get(curated[id].conceptId)??id,conceptId:curated[id].conceptId}));
assert.ok(elements.every(e=>concepts.some(c=>c.id===e.conceptId&&c.elements.includes(e.id))));
await fs.writeFile(path.join(out,'concept-map.json'),JSON.stringify({elements,concepts}));
await fs.writeFile(path.join(out,'system-map.json'),JSON.stringify({systems:Object.fromEntries(elements.map(e=>[e.id,curated[e.id].system]))}));
console.log(`Joined ${elements.length} publisher OBJ identities and ${concepts.length} publisher concepts.`);
