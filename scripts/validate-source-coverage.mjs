import fs from 'node:fs';
import assert from 'node:assert/strict';
const read=name=>JSON.parse(fs.readFileSync(new URL(name,import.meta.url)));
const catalogue=read('./data/male-source-catalog.json'),aliases=read('./data/male-source-aliases.json'),groups=read('./data/male-representation-groups.json'),synonyms=read('./data/male-representation-aliases.json');
const index=read('./data/male-source-index.json'),expected=new Set();
for(const region of index.regions){const b=fs.readFileSync(`/tmp/male-atlas-source/${region.region}.glb`),doc=JSON.parse(b.subarray(20,20+b.readUInt32LE(12)));
 function visit(i,owner){const node=doc.nodes[i];const candidate=aliases[`${region.region}|${node.name}`]??node.name;if(catalogue[candidate])owner=candidate;if(node.mesh!==undefined&&owner)for(const primitive of doc.meshes[node.mesh].primitives)expected.add(catalogue[owner].split?.[doc.materials?.[primitive.material]?.name]??owner);for(const child of node.children??[])visit(child,owner);}
 for(const root of doc.scenes[doc.scene??0].nodes)visit(root,null);
}
for(const file of ['atlas-male-detail.json','atlas-male-full.json']){
 const atlas=read('../public/models/'+file),parts=new Set(atlas.parts.map(p=>p.sourceId));assert.equal(atlas.physicalParts,4350);assert.equal(atlas.parts.filter(p=>p.derivedMaterial).length,254);const missing=[...expected].filter(id=>!parts.has(id));assert.deepEqual(missing,[],`${file}: source-owned geometry or material-defined patch missing`);
 for(const [name,members] of Object.entries(groups)){const available=members.filter(m=>parts.has(m));if(available.length){const c=atlas.concepts.find(c=>c.id==='DETAIL:'+name);assert.ok(c);assert.ok(available.every(n=>c.elements.includes('DETAIL:'+n)));}}
 for(const [name,target] of Object.entries(synonyms))if(atlas.concepts.some(c=>c.id==='DETAIL:'+target))assert.ok(atlas.concepts.some(c=>c.id==='DETAIL:'+name));
 console.log(`${file}: all ${expected.size} source-identifiable structures, including material-defined patches, whole-structure groups and aliases are selectable.`);
}
const missingCatalogue=Object.keys(catalogue).filter(n=>!expected.has(n));fs.writeFileSync(new URL('../public/SOURCE-COVERAGE.json',import.meta.url),JSON.stringify({source:'https://anatomy-atlas.brianp.chatgpt.site/',audit:'Source GLB ownership using the reference node-alias table and material-defined patches. Names without geometry are not fabricated.',geometryBackedRecords:expected.size,physicalRepresentations:4350,catalogueEntriesWithoutDirectGeometry:missingCatalogue},null,2));
