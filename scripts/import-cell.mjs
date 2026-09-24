#!/usr/bin/env node
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {gzipSync} from 'node:zlib';
import {Matrix3,Matrix4,Quaternion,Vector3} from 'three';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const source=path.join(root,'assets/source/human_eukaryotic_cell_and_its_main_components.glb');
const destination=path.join(root,'public/models');
const sourceHash='826a9b440515176cfa1a5e8d567d73a215dd3ce48c5602a35fdc18dcb18dc2ea';
const sourceUrl='https://sketchfab.com/3d-models/human-eukaryotic-cell-and-its-main-components-e769645f40fc47c38b89bcfc6b8b4d44';
const components=[
 ['Nuclear envelope','cell-nucleus','Nucleus'],
 ['Nuclear pores','cell-nucleus','Nucleus'],
 ['Nucleolus','cell-nucleus','Nucleus'],
 ['Heterochromatin','cell-nucleus','Nucleus'],
 ['Microfilaments (set 2)','cell-cytoskeleton','Cytoplasm'],
 ['Microfilaments','cell-cytoskeleton','Cytoplasm'],
 ['Golgi apparatus','cell-endomembrane','Cytoplasm'],
 ['Rough endoplasmic reticulum','cell-endomembrane','Cytoplasm'],
 ['Ribosomes','cell-protein-synthesis','Cytoplasm'],
 ['Proteins','cell-protein-synthesis','Cytoplasm'],
 ['Mitochondrial outer membranes','cell-mitochondria','Cytoplasm'],
 ['Endosome contents','cell-vesicles','Cytoplasm'],
 ['Endosome outer membrane','cell-vesicles','Cytoplasm'],
 ['Endosome inner membrane','cell-vesicles','Cytoplasm'],
 ['Centrioles','cell-division','Cytoplasm'],
 ['Plasma membrane','cell-boundary','Cell boundary'],
 ['Lysosomal membrane','cell-vesicles','Cytoplasm'],
 ['Lysosomal contents','cell-vesicles','Cytoplasm'],
 ['Mitochondrial interior','cell-mitochondria','Cytoplasm'],
 ['RNA','cell-protein-synthesis','Cytoplasm'],
];

const file=await fs.readFile(source);
assert.equal(createHash('sha256').update(file).digest('hex'),sourceHash,'The source GLB has changed');
assert.equal(file.toString('ascii',0,4),'glTF');
assert.equal(file.readUInt32LE(4),2);
assert.equal(file.readUInt32LE(8),file.length);
let cursor=12,json,bin;
while(cursor<file.length){
 const size=file.readUInt32LE(cursor),kind=file.readUInt32LE(cursor+4);
 if(kind===0x4e4f534a)json=JSON.parse(file.toString('utf8',cursor+8,cursor+8+size));
 if(kind===0x004e4942)bin=file.subarray(cursor+8,cursor+8+size);
 cursor+=8+size;
}
assert.ok(json&&bin);
assert.equal(json.meshes.length,components.length);
assert.equal(json.asset.extras.source,sourceUrl);
assert.equal(json.asset.extras.license,'CC-BY-NC-SA-4.0 (http://creativecommons.org/licenses/by-nc-sa/4.0/)');

const componentBytes={5120:1,5121:1,5122:2,5123:2,5125:4,5126:4};
const componentsIn={SCALAR:1,VEC2:2,VEC3:3,VEC4:4};
const readComponent=(view,offset,type)=>type===5120?view.getInt8(offset):type===5121?view.getUint8(offset):type===5122?view.getInt16(offset,true):type===5123?view.getUint16(offset,true):type===5125?view.getUint32(offset,true):view.getFloat32(offset,true);
function accessor(index){
 const a=json.accessors[index],b=json.bufferViews[a.bufferView],width=componentsIn[a.type],size=componentBytes[a.componentType],stride=b.byteStride??width*size,base=(b.byteOffset??0)+(a.byteOffset??0),view=new DataView(bin.buffer,bin.byteOffset,bin.byteLength);
 assert.ok(!a.sparse&&base+(a.count-1)*stride+width*size<=bin.length);
 const values=new Array(a.count*width);
 for(let i=0;i<a.count;i++)for(let j=0;j<width;j++)values[i*width+j]=readComponent(view,base+i*stride+j*size,a.componentType);
 return {values,width,count:a.count};
}
const instances=[];
function walk(index,parent){
 const node=json.nodes[index],local=node.matrix?new Matrix4().fromArray(node.matrix):new Matrix4().compose(new Vector3().fromArray(node.translation??[0,0,0]),new Quaternion().fromArray(node.rotation??[0,0,0,1]),new Vector3().fromArray(node.scale??[1,1,1]));
 const world=parent.clone().multiply(local);
 if(node.mesh!==undefined)instances.push({mesh:node.mesh,node:node.name,world});
 for(const child of node.children??[])walk(child,world);
}
for(const index of json.scenes[json.scene??0].nodes)walk(index,new Matrix4());
assert.equal(instances.length,components.length);
const raw=new Map(),min=new Vector3(Infinity,Infinity,Infinity),max=new Vector3(-Infinity,-Infinity,-Infinity);
for(const {mesh,node,world} of instances){
 const primitive=json.meshes[mesh].primitives[0];
 assert.equal(json.meshes[mesh].primitives.length,1);
 assert.equal(primitive.mode??4,4);
 const position=accessor(primitive.attributes.POSITION),normal=accessor(primitive.attributes.NORMAL),sourceColor=accessor(primitive.attributes.COLOR_0),indices=accessor(primitive.indices);
 assert.equal(position.width,3);assert.equal(normal.width,3);assert.equal(sourceColor.width,4);assert.equal(indices.width,1);
 const points=new Float32Array(position.values.length),normals=new Int16Array(normal.values.length),colors=new Uint8Array(sourceColor.count*3),normalMatrix=new Matrix3().getNormalMatrix(world),v=new Vector3();
 for(let i=0;i<position.count;i++){
  v.fromArray(position.values,i*3).applyMatrix4(world);points.set(v.toArray(),i*3);min.min(v);max.max(v);
  v.fromArray(normal.values,i*3).applyMatrix3(normalMatrix).normalize();for(let j=0;j<3;j++)normals[i*3+j]=Math.round(v.getComponent(j)*32767);
  for(let j=0;j<3;j++)colors[i*3+j]=Math.round(Math.max(0,Math.min(1,sourceColor.values[i*4+j]))*255);
 }
 assert.ok(indices.values.every(value=>value<position.count));
 raw.set(mesh,{node,primitive,points,normals,colors,indices:Uint32Array.from(indices.values)});
}
const center=min.clone().add(max).multiplyScalar(.5),scale=1.85/Math.max(...max.clone().sub(min).toArray());
const segments=[];let byteLength=0,triangles=0;
function append(data){
 const padding=(4-byteLength%4)%4;if(padding){segments.push(Buffer.alloc(padding));byteLength+=padding;}
 const offset=byteLength,bytes=Buffer.from(data.buffer,data.byteOffset,data.byteLength);segments.push(bytes);byteLength+=bytes.length;return offset;
}
const linearToSrgb=value=>value<=.0031308?value*12.92:1.055*Math.pow(value,1/2.4)-.055;
const parts=[],concepts=[],materials={};
for(let mesh=0;mesh<components.length;mesh++){
 const [name,system,compartment]=components[mesh],item=raw.get(mesh);
 assert.ok(item,`Missing source mesh ${mesh}`);
 const boundsMin=new Vector3(Infinity,Infinity,Infinity),boundsMax=new Vector3(-Infinity,-Infinity,-Infinity),point=new Vector3();
 for(let i=0;i<item.points.length;i+=3){
  point.fromArray(item.points,i).sub(center).multiplyScalar(scale);point.y+=.85;
  item.points.set(point.toArray(),i);boundsMin.min(point);boundsMax.max(point);
 }
 const material=json.materials[item.primitive.material],diffuse=material.extensions?.KHR_materials_pbrSpecularGlossiness?.diffuseFactor??[1,1,1,1];
 const materialId=`cell-material-${mesh}`;
 materials[materialId]={color:diffuse.slice(0,3).map(value=>Math.round(linearToSrgb(Math.max(0,Math.min(1,value)))*255)),opacity:diffuse[3],vertexColors:true};
 const id=`CELL:${mesh}`,positions=append(item.points),normals=append(item.normals),colors=append(item.colors),indices=append(item.indices);
 parts.push({id,name,conceptId:id,system,chunk:0,positions,normals,colors,indices,vertexCount:item.points.length/3,indexCount:item.indices.length,bounds:[boundsMin.toArray(),boundsMax.toArray()],material:materialId,sourceId:item.node,groups:[compartment]});
 concepts.push({id,name,elements:[id]});triangles+=item.indices.length/3;
}
for(const [id,name,meshes] of [
 ['cell','Human eukaryotic cell',components.map((_,i)=>i)],
 ['nucleus','Nucleus',[0,1,2,3]],
 ['microfilaments','Microfilaments (all)',[4,5]],
 ['mitochondria','Mitochondria',[10,18]],
 ['endosome','Endosome',[11,12,13]],
 ['lysosome','Lysosome',[16,17]],
])concepts.push({id:`CELL:${id}`,name,elements:meshes.map(i=>`CELL:${i}`)});
const binary=Buffer.concat(segments,byteLength),compressed=gzipSync(binary,{level:9});
await fs.mkdir(destination,{recursive:true});
await fs.writeFile(path.join(destination,'cell-0.bin'),binary);
await fs.writeFile(path.join(destination,'cell-0.bin.gz'),compressed);
const manifest={version:'E-learning UMCG human eukaryotic cell',scope:'cell',source:'E-learning UMCG via Sketchfab',sourceUrl,sourceSha256:sourceHash,license:'CC BY-NC-SA 4.0',licenseUrl:'https://creativecommons.org/licenses/by-nc-sa/4.0/',author:'E-learning UMCG',parts,concepts,materials,chunks:[{url:'/models/cell-0.bin',bytes:binary.length,gzip:'/models/cell-0.bin.gz',gzipBytes:compressed.length}],triangles};
await fs.writeFile(path.join(destination,'atlas-cell.json'),JSON.stringify(manifest));
console.log(`Imported ${parts.length} cell components, ${concepts.length} concepts, ${triangles.toLocaleString()} triangles; ${(binary.length/1e6).toFixed(1)} MB binary, ${(compressed.length/1e6).toFixed(1)} MB gzip`);
