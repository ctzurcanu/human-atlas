import {readFile,writeFile} from 'node:fs/promises';
import path from 'node:path';
import {gzipSync} from 'node:zlib';
import {MeshoptSimplifier} from 'meshoptimizer/simplifier';

// Separate assets keep the original anatomy and desktop rendering untouched.
await MeshoptSimplifier.ready;
const files=process.argv.slice(2);
if(!files.length)files.push(...['atlas-male-complete','atlas','atlas-hra-female','atlas-embryo','atlas-cell'].map(name=>`public/models/${name}.json`));
for(const filename of files){
 const atlas=JSON.parse(await readFile(filename,'utf8')),directory=path.dirname(filename),name=path.basename(filename,'.json');
 const sourceTriangles=atlas.parts.reduce((sum,part)=>sum+part.indexCount/3,0),budget=Math.min(350_000,sourceTriangles);
 const minimum=atlas.parts.reduce((sum,part)=>sum+Math.min(12,part.indexCount/3),0);
 const weights=atlas.parts.map(part=>Math.pow(part.indexCount/3,.8)),weightSum=weights.reduce((a,b)=>a+b,0);
 const parts=new Array(atlas.parts.length),segments=[];let bytes=0,triangles=0,maxError=0;
 const append=values=>{const padding=(4-bytes%4)%4;if(padding){segments.push(Buffer.alloc(padding));bytes+=padding;}const offset=bytes,buffer=Buffer.from(values.buffer,values.byteOffset,values.byteLength);segments.push(buffer);bytes+=buffer.byteLength;return offset;};
 for(let chunk=0;chunk<atlas.chunks.length;chunk++){
  const buffer=await readFile(path.join(directory,path.basename(atlas.chunks[chunk].url)));
  for(let i=0;i<atlas.parts.length;i++){
   const part=atlas.parts[i];if(part.chunk!==chunk)continue;
   const sourcePositions=new Float32Array(buffer.buffer,buffer.byteOffset+part.positions,part.vertexCount*3),sourceNormals=new Int16Array(buffer.buffer,buffer.byteOffset+part.normals,part.vertexCount*3);
   const original=new Uint32Array(buffer.buffer,buffer.byteOffset+part.indices,part.indexCount);
   // Weld split triangle vertices so disconnected GLB exports can simplify.
   const remap=MeshoptSimplifier.generatePositionRemap(sourcePositions,3),weldMap=new Map(),vertices=[],normalSums=[];
   const welded=new Uint32Array(part.vertexCount);
   for(let v=0;v<part.vertexCount;v++){
    let index=weldMap.get(remap[v]);
    if(index===undefined){index=vertices.length/3;weldMap.set(remap[v],index);vertices.push(...sourcePositions.subarray(v*3,v*3+3));normalSums.push(0,0,0);}
    welded[v]=index;for(let axis=0;axis<3;axis++)normalSums[index*3+axis]+=sourceNormals[v*3+axis];
   }
   const positions=new Float32Array(vertices),normals=new Int16Array(normalSums.length),indices=Uint32Array.from(original,index=>welded[index]);
   for(let v=0;v<normals.length;v+=3){const length=Math.hypot(...normalSums.slice(v,v+3))||1;for(let axis=0;axis<3;axis++)normals[v+axis]=Math.round(normalSums[v+axis]/length*32767);}
   const target=Math.min(part.indexCount/3,Math.max(12,Math.floor((budget-minimum)*weights[i]/weightSum)+Math.min(12,part.indexCount/3)))*3;
   const [simplified,error]=MeshoptSimplifier.simplify(indices,positions,3,target,.03,['Permissive']);
   if(!simplified.length)throw new Error(`VR simplification removed ${part.id}`);
   maxError=Math.max(maxError,error);
   const [compact,count]=MeshoptSimplifier.compactMesh(simplified),compactPositions=new Float32Array(count*3),compactNormals=new Int16Array(count*3);
   for(let v=0;v<compact.length;v++){const index=compact[v];if(index===0xffffffff)continue;for(let axis=0;axis<3;axis++){compactPositions[index*3+axis]=positions[v*3+axis]+(part.sourceOffset?.[axis]??0);compactNormals[index*3+axis]=normals[v*3+axis];}}
   parts[i]={id:part.id,sourceVertexCount:part.vertexCount,sourceIndexCount:part.indexCount,positions:append(compactPositions),normals:append(compactNormals),indices:append(simplified),vertexCount:count,indexCount:simplified.length};
   triangles+=simplified.length/3;
  }
 }
 const binary=Buffer.concat(segments),compressed=gzipSync(binary,{level:9}),url=atlas.chunks[0].url.replace(/[^/]+$/,`${name}.vr.bin.gz`);
 await writeFile(path.join(directory,`${name}.vr.bin.gz`),compressed);
 await writeFile(path.join(directory,`${name}.vr.json`),JSON.stringify({version:1,sourceTriangles,triangles,maxRelativeError:maxError,parts,chunk:{url,bytes,gzip:true}}));
 console.log(JSON.stringify({model:name,parts:parts.length,sourceTriangles,triangles,maxError,bytes,gzipBytes:compressed.length}));
}
