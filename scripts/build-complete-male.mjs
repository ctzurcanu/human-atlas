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

// Restrict the supplement to named gaps so matching geometry in the two
// sources is not duplicated.
const digestive=/ileu|small intestine|rectum|pancreatic duct tree|hepatic duct|cystic duct|biliary tree|mesentery of small intestine|ileocecal junction/i;
const pelvic=/deferent duct|superficial dorsal vein of penis/i;
const added=standard.parts.filter(part=>part.system==='digestive'&&digestive.test(part.name)||['reproductive','arterial','venous'].includes(part.system)&&pelvic.test(part.name));
assert.equal(detail.parts.length,5254,'Unexpected detailed-source catalogue');
assert.equal(standard.version,'BodyParts3D 4.0','Unexpected supplement source');
assert.equal(added.length,56,'The BodyParts3D supplement changed; review the merge');
const existingParts=new Set(detail.parts.map(part=>part.id));
const existingConcepts=new Set(detail.concepts.map(concept=>concept.id));
for(const part of added)assert(!existingParts.has(part.id),`Duplicate part ${part.id}`);

const selectedIds=new Set(added.map(part=>part.id));
const usedChunks=[...new Set(added.map(part=>part.chunk))].sort((a,b)=>a-b);
const remap=new Map(usedChunks.map((chunk,index)=>[chunk,detail.chunks.length+index]));
const provenance={
 label:'BodyParts3D 4.0 original publisher OBJ',
 url:standard.provenance.sourceUrl,
 detail:'CC BY 4.0; original publisher geometry, optimized by the Human Atlas import script.',
};
const center=part=>part.bounds[0].map((min,axis)=>(min+part.bounds[1][axis])/2);
const alignment=['Prostate','Corpus cavernosum of penis'].map(name=>{
 const source=standard.parts.find(part=>part.name===name),target=detail.parts.find(part=>part.name===name);
 assert(source&&target,`Missing alignment landmark ${name}`);
 return center(target).map((value,axis)=>value-center(source)[axis]);
});
const pelvicOffset=[0,1,2].map(axis=>alignment.reduce((sum,offset)=>sum+offset[axis],0)/alignment.length);
const parts=added.map(part=>{
 const offset=part.system==='digestive'?null:pelvicOffset;
 return {...part,chunk:remap.get(part.chunk),sourceId:part.id,license:'CC BY 4.0',provenance,
  ...(offset?{sourceOffset:offset,bounds:part.bounds.map(point=>point.map((value,axis)=>value+offset[axis]))}:{})};
});
const concepts=standard.concepts
 .map(concept=>({...concept,elements:concept.elements.filter(id=>selectedIds.has(id))}))
 .filter(concept=>concept.elements.length&&(digestive.test(concept.name)||pelvic.test(concept.name)||parts.some(part=>part.conceptId===concept.id))&&!existingConcepts.has(concept.id));
const small=concepts.find(concept=>concept.name.toLowerCase()==='small intestine');
assert(small,'Missing small intestine concept');
small.elements=[...detail.concepts.find(concept=>concept.id==='ZA:Duodenum').elements,...detail.concepts.find(concept=>concept.id==='ZA:Jejunum').elements,...small.elements];
const large=[...['ZA:Ascending colon','ZA:Transverse colon','ZA:Descending colon','ZA:Sigmoid colon','ZA:Vermiform appendix'].flatMap(id=>detail.concepts.find(concept=>concept.id===id).elements),...parts.filter(part=>part.name==='Rectum').map(part=>part.id)];
concepts.push({id:'HA:large-intestine',name:'Large intestine',elements:large});
concepts.push({id:'HA:intestines',name:'Intestines',elements:[...small.elements,...large]});
concepts.push({id:'HA:ureters',name:'Ureters',elements:['ZA:Ureter.l','ZA:Ureter.r']});
const ducts=parts.filter(part=>/deferent duct/i.test(part.name)).map(part=>part.id);
const penileVessels=[...detail.parts,...parts].filter(part=>/dorsal (artery|vein) of penis/i.test(part.name)).map(part=>part.id);
assert.equal(ducts.length,2);
assert.equal(penileVessels.length,6);
concepts.push({id:'HA:ductus-deferens',name:'Ductus deferens',elements:ducts});
concepts.push({id:'HA:penile-blood-vessels',name:'Blood vessels of the penis',elements:penileVessels});
const conceptIds=new Set(concepts.map(concept=>concept.id));
for(const part of parts)assert(conceptIds.has(part.conceptId),`Missing concept for ${part.id}`);
for(const name of ['jejunum','ileum','small intestine','rectum','oesophagus','ureters','ductus deferens','blood vessels of the penis'])
 assert([...detail.concepts,...concepts].some(concept=>concept.name.toLowerCase()===name),`Missing ${name} concept`);

const atlas={
 ...detail,
 version:`${detail.version} + BodyParts3D 4.0 targeted supplement`,
 source:`${detail.source} + BodyParts3D 4.0`,
 scope:`${detail.scope} Missing digestive and pelvic structures are supplemented from the original BodyParts3D 4.0 publisher geometry.`,
 parts:[...detail.parts,...parts],
 concepts:[...detail.concepts,...concepts],
 chunks:[...detail.chunks,...usedChunks.map(index=>standard.chunks[index])],
 triangles:detail.triangles+parts.reduce((total,part)=>total+part.indexCount/3,0),
 provenance:{...detail.provenance,additionalSources:[...(detail.provenance.additionalSources??[]),{...standard.provenance,license:'CC BY 4.0',supplement:'Missing digestive and pelvic structures'}]},
};
const output=path.join(models,'atlas-male-complete.json');
await writeFile(output,JSON.stringify(atlas));
console.log(`Built ${output}: ${atlas.parts.length} surfaces, ${atlas.concepts.length} concepts; ${parts.length} BodyParts3D meshes added.`);
