/** Separate the HRA placenta branch from the female body, including its meshes. */
import fs from 'node:fs/promises';
import path from 'node:path';
import {gzipSync} from 'node:zlib';
import {fileURLToPath} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const directory=path.resolve(process.argv[2]??path.join(root,'public/models'));
const femaleFile=path.join(directory,'atlas-hra-female.json');
const embryoFile=path.join(directory,'atlas-embryo.json');
const atlas=JSON.parse(await fs.readFile(femaleFile,'utf8'));
const placenta=atlas.concepts.find(concept=>concept.id==='HRA:VH_F_placenta');
if(!placenta||placenta.elements.length!==8)throw new Error('Expected the eight original HRA placenta meshes before splitting.');
const embryoIds=new Set(placenta.elements);
const embryoParts=atlas.parts.filter(part=>embryoIds.has(part.id));
if(embryoParts.length!==8)throw new Error('Placenta concept and meshes disagree.');
const bodyParts=atlas.parts.filter(part=>!embryoIds.has(part.id));
const inputChunks=atlas.chunks.map(chunk=>path.join(directory,path.basename(chunk.url)));
const buffers=await Promise.all(inputChunks.map(file=>fs.readFile(file)));

function pack(parts,prefix){
 const outputChunks=[],outputParts=[],files=[];
 let segments=[],bytes=0;
 const append=data=>{const padding=(4-bytes%4)%4;if(padding){segments.push(Buffer.alloc(padding));bytes+=padding;}const offset=bytes;segments.push(data);bytes+=data.length;return offset;};
 const flush=()=>{if(!bytes)return;const name=`${prefix}-${outputChunks.length}.bin`,raw=Buffer.concat(segments),gz=gzipSync(raw,{level:9,mtime:0});files.push([name,raw],[name+'.gz',gz]);outputChunks.push({url:'/models/'+name,bytes:raw.length,gzip:'/models/'+name+'.gz',gzipBytes:gz.length});segments=[];bytes=0;};
 for(const part of parts){
  if(bytes>6_000_000)flush();
  const source=buffers[part.chunk];
  const positions=append(source.subarray(part.positions,part.positions+part.vertexCount*12));
  const normals=append(source.subarray(part.normals,part.normals+part.vertexCount*6));
  const indices=append(source.subarray(part.indices,part.indices+part.indexCount*4));
  outputParts.push({...part,chunk:outputChunks.length,positions,normals,indices,...(prefix==='embryo'?{system:'pregnancy'}:{})});
 }
 flush();
 return {parts:outputParts,chunks:outputChunks,files};
}
const body=pack(bodyParts,'female-body'),embryo=pack(embryoParts,'embryo');
const conceptsFor=(ids)=>atlas.concepts.map(concept=>({...concept,elements:concept.elements.filter(id=>ids.has(id))})).filter(concept=>concept.elements.length);
const bodyIds=new Set(bodyParts.map(part=>part.id));
const bodyAtlas={...atlas,parts:body.parts,concepts:conceptsFor(bodyIds),chunks:body.chunks,triangles:body.parts.reduce((sum,part)=>sum+part.indexCount/3,0),scope:'Female body anatomy; the separate Embryo model contains the HRA placenta branch.'};
const embryoConcepts=conceptsFor(embryoIds).filter(concept=>concept.id!=='HRA:VH_F'&&concept.id!=='HRA:VH_F_reproductive_system');
embryoConcepts.push({id:'HRA:embryo',name:'Embryo',elements:[...embryoIds]});
const embryoAtlas={...atlas,version:'Human Reference Atlas united-female v1.10 placenta branch',scope:'embryo',parts:embryo.parts,concepts:embryoConcepts,chunks:embryo.chunks,materials:Object.fromEntries(Object.entries(atlas.materials).filter(([id])=>embryo.parts.some(part=>part.material===id))),triangles:embryo.parts.reduce((sum,part)=>sum+part.indexCount/3,0),provenance:{sourceUrl:atlas.provenance.sourceUrl,sourceSha256:atlas.provenance.sourceSha256,sourceFile:atlas.provenance.sourceFile,license:atlas.provenance.license,notes:'Original HRA placenta branch: basal plate, chorionic plate, placental vessels, amnion, and umbilical cord with its vessels. The source does not contain a fetal body.'}};
for(const [name,data] of [...body.files,...embryo.files])await fs.writeFile(path.join(directory,name),data);
await fs.writeFile(embryoFile,JSON.stringify(embryoAtlas));
await fs.writeFile(femaleFile,JSON.stringify(bodyAtlas));
for(const file of inputChunks){
 if(path.basename(file).startsWith('hra-female-')){
  await fs.rm(file);
  await fs.rm(file+'.gz');
 }
}
console.log(JSON.stringify({female:body.parts.length,embryo:embryo.parts.length,femaleChunks:body.chunks.length,embryoChunks:embryo.chunks.length}));
