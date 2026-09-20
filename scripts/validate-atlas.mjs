import fs from 'node:fs';
import assert from 'node:assert/strict';
const filename=process.argv[2]??'atlas.json',female=filename.includes('female'),full=filename.includes('male-full'),detailed=filename.includes('male-detail')||full;
const base=new URL('../public/models/',import.meta.url),atlas=JSON.parse(fs.readFileSync(new URL(filename,base)));
if(!detailed){assert.equal(atlas.parts.length,female?1038:2234);assert.equal(atlas.concepts.length,female?1253:3432);}
const ids=new Set(atlas.parts.map(p=>p.id));assert.equal(ids.size,atlas.parts.length);
const files=atlas.chunks.map(c=>{const b=fs.readFileSync(new URL(c.url.split('/').pop(),base));assert.equal(b.length,c.bytes);return b;});
if(female){assert.equal(atlas.parts.filter(p=>p.system==='pregnancy').length,8);for(const label of ['uterus','ovary','vagina'])assert.ok(atlas.parts.some(p=>p.name.toLowerCase().includes(label)));assert.ok(!atlas.parts.some(p=>/prostate|testis|penis/i.test(p.name)));}
if(female){
 assert.equal(atlas.provenance.femaleOnly,true);
 assert.equal(atlas.parts.filter(p=>p.system==='muscular').length,90);
 assert.equal(atlas.parts.filter(p=>p.id.startsWith('VHF')).length,76);
 assert.ok(atlas.parts.every(p=>/^HRAF\d+$|^VHF\d+$/.test(p.id)));
 assert.ok(atlas.parts.every(p=>p.provenance?.label&&p.provenance?.url));
 assert.ok(!atlas.parts.some(p=>p.system==='borrowed'||p.id==='HRAF0394'||p.id==='HRAF0396'));
 for(const side of ['left','right'])assert.equal(atlas.parts.filter(p=>p.name===`Rectus femoris (${side})`).length,1);
 for(const name of ['Gluteus maximus (left)','Soleus (right)','Sartorius (left)'])assert.ok(atlas.parts.some(p=>p.name===name));
}
if(full){assert.equal(atlas.triangles,atlas.sourceTriangles);assert.equal(atlas.optimized.maximumRelativeError,0);assert.ok(atlas.parts.some(p=>p.sourceId==='Body surface (derived)'&&p.system==='integumentary'));}
if(detailed){
 assert.equal(atlas.sex,'male');assert.equal(atlas.parts.length,4391);assert.equal(atlas.concepts.length,4683);
 assert.ok(atlas.parts.filter(p=>p.schematic).every(p=>p.system==='schematic'));
 assert.ok(atlas.parts.filter(p=>p.tissue==='fascia').every(p=>p.system==='fascia'));
 assert.ok(atlas.parts.filter(p=>p.source==='o3m').length>650);
 assert.ok(atlas.parts.filter(p=>p.system==='attachments').length>500);
 assert.ok(atlas.parts.filter(p=>p.system==='muscular').length>500);
 assert.ok(atlas.parts.every(p=>p.provenance?.label&&p.license));
 assert.equal(atlas.parts.filter(p=>p.license==='CC BY-NC 4.0').length,2);
 assert.equal(atlas.parts.filter(p=>p.license==='CC BY-NC-SA 4.0').length,4);
 for(const side of ['l','r'])for(const name of ['Median nerve','Femur','Kidney'])assert.ok(atlas.parts.some(p=>p.sourceId===`${name}.${side}`),`${name}.${side} missing`);
 const systems=new Set(['skeletal','muscular','arterial','venous','nervous','digestive','respiratory','urinary','reproductive','lymphatic','endocrine','integumentary','connective','sensory','cardiac','attachments','regions','fascia','schematic']);
 assert.ok(atlas.parts.every(p=>systems.has(p.system)));
}
let tris=0;
for(const p of atlas.parts){assert.ok(p.name.trim()&&p.name!=='-'&&!p.name.includes('Bounds('));assert.ok(p.conceptId!=='-');assert.ok(p.system);const b=files[p.chunk];assert.ok(p.indices+p.indexCount*4<=b.length);const pos=new Float32Array(b.buffer,b.byteOffset+p.positions,p.vertexCount*3),indices=new Uint32Array(b.buffer,b.byteOffset+p.indices,p.indexCount);assert.ok(indices.length>=3);for(const i of indices)assert.ok(i<p.vertexCount,`${p.id}: invalid vertex`);for(const value of pos)assert.ok(Number.isFinite(value));tris+=p.indexCount/3;}
assert.equal(new Set(atlas.concepts.map(c=>c.id)).size,atlas.concepts.length);
for(const c of atlas.concepts){assert.ok(c.elements.length);for(const id of c.elements)assert.ok(ids.has(id),`${c.id}: missing ${id}`);}
assert.equal(tris,atlas.triangles);
console.log(`Verified ${ids.size} individually indexed meshes, ${atlas.concepts.length} complete concept mappings, ${tris.toLocaleString()} triangles, and every binary buffer.`);
