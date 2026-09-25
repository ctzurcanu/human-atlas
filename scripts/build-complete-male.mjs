/** Add missing geometry from the original BodyParts3D catalogue to the
 * original Z-Anatomy/Open 3D Model catalogue. Both inputs are built by this repo.
 */
import assert from 'node:assert/strict';
import {readFile,writeFile} from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const models=process.argv[2]?path.resolve(process.argv[2]):path.join(root,'public/models');
const detail=JSON.parse(await readFile(path.join(models,'atlas-z-anatomy.json'),'utf8'));
const standard=JSON.parse(await readFile(path.join(models,'atlas.json'),'utf8'));

// The Z-Anatomy jejunum is a detailed continuous small-bowel surface. The
// BodyParts3D ileum, ileocecal junction, and mesentery belong to a different
// bowel representation: combining them creates duplicate, crossing loops.
// The detailed model also contains its own pancreatic duct, deferent ducts,
// and dorsal penile veins. Supplement only structures missing from that model.
const digestive=/rectum|hepatic duct|cystic duct|biliary tree/i;
const added=standard.parts.filter(part=>part.system==='digestive'&&digestive.test(part.name));
const replacementStomach=standard.parts.find(part=>part.name==='Stomach');
const sourceStomach=detail.parts.find(part=>part.id==='ZA:Stomach');
assert(replacementStomach&&sourceStomach,'Missing stomach geometry for wall repair');
const stomachBuffer=await readFile(path.join(models,path.basename(standard.chunks[replacementStomach.chunk].url)));
const stomachIndices=new Uint32Array(stomachBuffer.buffer,stomachBuffer.byteOffset+replacementStomach.indices,replacementStomach.indexCount);
const stomachEdges=new Map();
for(let i=0;i<stomachIndices.length;i+=3)for(let side=0;side<3;side++){
 const a=stomachIndices[i+side],b=stomachIndices[i+(side+1)%3],key=a<b?`${a},${b}`:`${b},${a}`;
 stomachEdges.set(key,(stomachEdges.get(key)??0)+1);
}
assert([...stomachEdges.values()].filter(count=>count===1).length<=24,'Replacement stomach has a large open boundary; review its source mesh');
assert.equal(detail.parts.length,5192,'Unexpected detailed-source catalogue');
assert.equal(standard.version,'BodyParts3D 4.0','Unexpected supplement source');
assert.equal(added.length,17,'The BodyParts3D supplement changed; review the merge');
const existingParts=new Set(detail.parts.map(part=>part.id));
const existingConcepts=new Set(detail.concepts.map(concept=>concept.id));
for(const part of added)assert(!existingParts.has(part.id),`Duplicate part ${part.id}`);

const selectedIds=new Set(added.map(part=>part.id));
const usedChunks=[...new Set([...added.map(part=>part.chunk),replacementStomach.chunk])].sort((a,b)=>a-b);
const remap=new Map(usedChunks.map((chunk,index)=>[chunk,detail.chunks.length+index]));
const provenance={
 label:'BodyParts3D 4.0 original publisher OBJ',
 url:standard.provenance.sourceUrl,
 detail:'CC BY 4.0; original publisher geometry, optimized by the Human Atlas import script.',
};
const parts=added.map(part=>({...part,chunk:remap.get(part.chunk),sourceId:part.id,license:'CC BY 4.0',provenance}));
// The Z-Anatomy stomach has a 39-edge opening in its anterior wall. The
// aligned BodyParts3D stomach has only its two small physiological ports.
// Keep the established ZA identifier and hierarchy while replacing its mesh.
const detailedParts=detail.parts.map(part=>part.id==='ZA:Stomach'?{
 ...part,chunk:remap.get(replacementStomach.chunk),positions:replacementStomach.positions,normals:replacementStomach.normals,
 indices:replacementStomach.indices,vertexCount:replacementStomach.vertexCount,indexCount:replacementStomach.indexCount,
 bounds:replacementStomach.bounds,sourceId:replacementStomach.id,license:'CC BY 4.0',
 provenance:{...provenance,detail:'CC BY 4.0; BodyParts3D publisher stomach surface replaces the Z-Anatomy surface with a large anterior-wall opening.'},
}:part);
const concepts=standard.concepts
 .map(concept=>({...concept,elements:concept.elements.filter(id=>selectedIds.has(id))}))
 .filter(concept=>concept.elements.length&&(digestive.test(concept.name)||parts.some(part=>part.conceptId===concept.id))&&!existingConcepts.has(concept.id));
const small={id:'HA:small-intestine',name:'Small intestine',elements:[...detail.concepts.find(concept=>concept.id==='ZA:Duodenum').elements,...detail.concepts.find(concept=>concept.id==='ZA:Jejunum').elements]};
concepts.push(small);
const large=[...['ZA:Ascending colon','ZA:Transverse colon','ZA:Descending colon','ZA:Sigmoid colon','ZA:Vermiform appendix'].flatMap(id=>detail.concepts.find(concept=>concept.id===id).elements),...parts.filter(part=>part.name==='Rectum').map(part=>part.id)];
concepts.push({id:'HA:large-intestine',name:'Large intestine',elements:large});
concepts.push({id:'HA:intestines',name:'Intestines',elements:[...small.elements,...large]});
concepts.push({id:'HA:ureters',name:'Ureters',elements:['ZA:Ureter.l','ZA:Ureter.r']});
const ducts=detail.parts.filter(part=>/ductus deferens/i.test(part.name)).map(part=>part.id);
const penileVessels=detail.parts.filter(part=>/dorsal arter(y|ies) of penis|dorsal veins? of penis/i.test(part.name)).map(part=>part.id);
assert.equal(ducts.length,2);
assert.equal(penileVessels.length,4);
concepts.push({id:'HA:ductus-deferens',name:'Ductus deferens',elements:ducts});
concepts.push({id:'HA:penile-blood-vessels',name:'Blood vessels of the penis',elements:penileVessels});
const conceptIds=new Set(concepts.map(concept=>concept.id));
for(const part of parts)assert(conceptIds.has(part.conceptId),`Missing concept for ${part.id}`);
for(const name of ['jejunum','small intestine','rectum','oesophagus','ureters','ductus deferens','blood vessels of the penis'])
 assert([...detail.concepts,...concepts].some(concept=>concept.name.toLowerCase()===name),`Missing ${name} concept`);

const atlas={
 ...detail,
 version:`${detail.version} + BodyParts3D 4.0 targeted supplement`,
 source:`${detail.source} + BodyParts3D 4.0`,
 scope:`${detail.scope} Missing digestive structures are supplemented from the original BodyParts3D 4.0 publisher geometry.`,
 parts:[...detailedParts,...parts],
 concepts:[...detail.concepts,...concepts],
 chunks:[...detail.chunks,...usedChunks.map(index=>standard.chunks[index])],
 triangles:detail.triangles+(replacementStomach.indexCount-sourceStomach.indexCount)/3+parts.reduce((total,part)=>total+part.indexCount/3,0),
 provenance:{...detail.provenance,additionalSources:[...(detail.provenance.additionalSources??[]),{...standard.provenance,license:'CC BY 4.0',supplement:'Missing digestive and pelvic structures'}]},
};
const output=path.join(models,'atlas-male-complete.json');
await writeFile(output,JSON.stringify(atlas));
console.log(`Built ${output}: ${atlas.parts.length} surfaces, ${atlas.concepts.length} concepts; ${parts.length} BodyParts3D meshes added.`);
