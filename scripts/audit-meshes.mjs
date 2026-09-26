/** Diagnose source geometry without altering publisher meshes. Run after imports. */
import {readFile,writeFile,mkdir} from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {build} from 'esbuild';
import * as T from 'three';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const models={male:'public/models/atlas.json','male-detail':'public/models/atlas-z-anatomy.json','male-full':'public/models/atlas-male-complete.json',female:'public/models/atlas-hra-female.json',embryo:'public/models/atlas-embryo.json',cell:'public/models/atlas-cell.json','local-reference':'.local-models/reference.json','local-male':'.local-models/male.json','local-female':'.local-models/female.json'};
const bundle=await build({entryPoints:[path.join(root,'app/mesh-topology.ts')],bundle:true,platform:'node',format:'esm',write:false});
const {analyzeMeshTopology}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);
const selected=process.argv.slice(2).filter(value=>!value.startsWith('--'));
for(const id of selected.length?selected:Object.keys(models))if(!models[id])throw Error(`Unknown model ${id}`);
const summaries={};
let missingMappedUvs=0;
for(const id of selected.length?selected:Object.keys(models)){
 const atlas=JSON.parse(await readFile(path.join(root,models[id]),'utf8'));
 const byChunk=new Map();atlas.parts.forEach((part,index)=>{const group=byChunk.get(part.chunk)??[];group.push({part,index});byChunk.set(part.chunk,group);});
 const flags=[],totals={parts:atlas.parts.length,triangles:0,open:0,nonManifold:0,reversedNormals:0,duplicateFaces:0,missingMappedUvs:0,outsideBodySurface:0,coincidentBounds:0};
 const shell=atlas.parts.filter(part=>/^(skin(?: of body)?|body surface)$/i.test(part.name)).sort((a,b)=>b.vertexCount-a.vertexCount)[0];
 const body=shell?.bounds,coincident=new Map();
 for(const [chunkIndex,entries] of byChunk){
  const url=atlas.chunks[chunkIndex].url;
  const file=url.startsWith('/models/')?path.join(root,'public',url):url.startsWith('/local-models/')?path.join(root,'.local-models',url.slice('/local-models/'.length)):path.join(root,path.dirname(models[id]),path.basename(url));
  const bytes=await readFile(file);
  for(const {part,index} of entries){
   const geometry=new T.BufferGeometry();
   geometry.setAttribute('position',new T.BufferAttribute(new Float32Array(bytes.buffer,bytes.byteOffset+part.positions,part.vertexCount*3),3));
   if(part.normals!==undefined)geometry.setAttribute('normal',new T.BufferAttribute(new Int16Array(bytes.buffer,bytes.byteOffset+part.normals,part.vertexCount*3),3,true));
   if(part.uvs!==undefined)geometry.setAttribute('uv',new T.BufferAttribute(new Float32Array(bytes.buffer,bytes.byteOffset+part.uvs,part.vertexCount*2),2));
   geometry.setIndex(new T.BufferAttribute(new Uint32Array(bytes.buffer,bytes.byteOffset+part.indices,part.indexCount),1));
   const t=analyzeMeshTopology(geometry);geometry.dispose();totals.triangles+=t.triangles;
   const sourceMaterial=atlas.materials?.[part.material];
   const mapped=!!(sourceMaterial?.map||sourceMaterial?.normalMap);
   const outside=body&&part!==shell&&part.bounds.some((corner,side)=>corner.some((v,axis)=>side===0?v<body[0][axis]-.06:v>body[1][axis]+.06));
   const bKey=part.bounds.flat().map(v=>Math.round(v*1000)).join(',')+':'+part.indexCount;
   const earlier=coincident.get(bKey);if(earlier!==undefined&&earlier!==part.conceptId){totals.coincidentBounds++;flags.push({index,id:part.id,name:part.name,issue:'coincident-bounds',peer:earlier});}else coincident.set(bKey,part.conceptId);
   if(t.boundaryEdges)totals.open++;
   if(t.nonManifoldEdges||t.windingConflicts)totals.nonManifold++;
   if(t.reversedNormals)totals.reversedNormals++;
   if(t.duplicateTriangles)totals.duplicateFaces++;
   if(mapped&&!t.hasUv)totals.missingMappedUvs++;
   if(outside)totals.outsideBodySurface++;
   if(t.boundaryEdges||t.nonManifoldEdges||t.windingConflicts||t.duplicateTriangles||t.reversedNormals||mapped&&!t.hasUv||outside)flags.push({index,id:part.id,name:part.name,system:part.system,topology:t,...(mapped&&!t.hasUv?{missingMappedUvs:true}:{}),...(outside?{outsideBodySurface:true}:{})});
  }
 }
 const report={model:id,source:models[id],version:atlas.version,totals,flags};
 const reportDir=path.join(root,'reports','mesh-audit');await mkdir(reportDir,{recursive:true});await writeFile(path.join(reportDir,`${id}.json`),JSON.stringify(report));
 summaries[id]=totals;console.log(`${id}: ${totals.parts} parts; ${totals.open} open, ${totals.nonManifold} nonmanifold, ${totals.reversedNormals} with reversed normals, ${totals.missingMappedUvs} mapped without UVs`);
 missingMappedUvs+=totals.missingMappedUvs;
}
console.log(JSON.stringify(summaries,null,2));
if(missingMappedUvs){console.error(`${missingMappedUvs} textured parts have no UVs; inspect reports/mesh-audit before deploying.`);process.exitCode=1;}
