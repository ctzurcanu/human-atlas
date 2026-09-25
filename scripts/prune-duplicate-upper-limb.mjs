/** Remove previously imported Open 3D Model anatomy represented in Z-Anatomy.
 * The current importer skips these; this also upgrades existing manifests.
 */
import assert from 'node:assert/strict';
import {readFile,writeFile} from 'node:fs/promises';
import path from 'node:path';

const filename=path.resolve(process.argv[2]??'public/models/atlas-z-anatomy.json');
const atlas=JSON.parse(await readFile(filename,'utf8'));
const bone=part=>part.system==='skeletal'&&!/Sesamoid bones of hand/i.test(part.name);
const superficialVein=part=>/^Arm superficial vein-(Basilic|Cephalic|Median antebrachial|Median cubital) vein/i.test(part.name);
const brachiocephalic=part=>/^Brachiocephalic artery/i.test(part.name);
const o3m=atlas.parts.filter(part=>part.id.startsWith('O3M:'));
for(const [name,count] of [[bone,52],[superficialVein,8],[brachiocephalic,2]]){
 const found=o3m.filter(name).length;
 assert.ok(found===count||found===0,`Unexpected Open 3D Model duplicate set: ${found} rather than ${count}`);
}
const redundant=o3m.filter(part=>bone(part)||superficialVein(part)||brachiocephalic(part));
assert.equal(atlas.parts.filter(part=>part.id.startsWith('O3M:')&&part.system==='skeletal'&&/Sesamoid bones of hand/i.test(part.name)).length,2);
if(redundant.length){
 const removed=new Set(redundant.map(part=>part.id));
 atlas.parts=atlas.parts.filter(part=>!removed.has(part.id));
 atlas.concepts=atlas.concepts.map(concept=>({...concept,elements:concept.elements.filter(id=>!removed.has(id))})).filter(concept=>concept.elements.length);
 atlas.triangles-=redundant.reduce((sum,part)=>sum+part.indexCount/3,0);
 await writeFile(filename,JSON.stringify(atlas));
}
console.log(`Removed ${redundant.length} redundant Open 3D Model surfaces from ${filename}.`);
