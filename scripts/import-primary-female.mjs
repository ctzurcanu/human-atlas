/** Convert official HRA united-female v1.10 with missing v1.5 pelvic pieces.
 * Usage: node scripts/import-primary-female.mjs [--source SOURCE.glb] [--supplement V1.5.glb] [--out DIRECTORY]
 * Both sources are downloaded from cdn.humanatlas.io when omitted.
 */
import fs from 'node:fs/promises';
import fssync from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {Matrix3,Matrix4,Quaternion,Vector3} from 'three';
import {MeshoptSimplifier} from 'meshoptimizer';

const SOURCE_URL='https://cdn.humanatlas.io/digital-objects/ref-organ/united-female/v1.10/assets/3d-vh-f-united.glb';
const SOURCE_SHA256='95f0c3d2f918582608692ca1139e8bdb18c147a16470e9ee9af8b276bd77c422';
const SUPPLEMENT_URL='https://cdn.humanatlas.io/digital-objects/ref-organ/united-female/v1.5/assets/3d-vh-f-united.glb';
const SUPPLEMENT_SHA256='472567a56896b9b7890508da6501fbf858e56aaa30745365f7a71ade782b529c';
const ROOT=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const option=(key,defaultValue)=>{const index=process.argv.indexOf(key);return index<0?defaultValue:path.resolve(process.argv[index+1]);};
const sourceFile=option('--source',path.join(os.tmpdir(),'human-atlas-primary-models','hra-united-female-v1.10.glb'));
const supplementFile=option('--supplement',path.join(os.tmpdir(),'human-atlas-primary-models','hra-united-female-v1.5.glb'));
const out=option('--out',path.join(ROOT,'public','models'));
await fs.mkdir(out,{recursive:true});
async function verifiedSource(filename,url,expected){
 await fs.mkdir(path.dirname(filename),{recursive:true});
 if(!fssync.existsSync(filename)){
  const temporary=filename+'.download';
  const child=spawn('curl',['--fail','--location','--silent','--show-error','--output',temporary,url],{stdio:'inherit'});
  const code=await new Promise((resolve,reject)=>{child.once('error',reject);child.once('close',resolve)});
  if(code!==0)throw new Error(`HRA source download failed (curl exit ${code})`);
  await fs.rename(temporary,filename);
 }
 const data=await fs.readFile(filename);
 assert.equal(createHash('sha256').update(data).digest('hex'),expected,`HRA source hash differs: ${filename}`);
 return data;
}
const glb=await verifiedSource(sourceFile,SOURCE_URL,SOURCE_SHA256);
const supplementGlb=await verifiedSource(supplementFile,SUPPLEMENT_URL,SUPPLEMENT_SHA256);
assert.equal(glb.readUInt32LE(0),0x46546c67);
assert.equal(glb.readUInt32LE(4),2);
const jsonLength=glb.readUInt32LE(12),doc=JSON.parse(glb.subarray(20,20+jsonLength).toString());
const binary=glb.subarray(28+jsonLength);
assert.equal(supplementGlb.readUInt32LE(0),0x46546c67);
const supplementJsonLength=supplementGlb.readUInt32LE(12);
const supplementDoc=JSON.parse(supplementGlb.subarray(20,20+supplementJsonLength).toString());
const supplementBinary=supplementGlb.subarray(28+supplementJsonLength);
assert.equal(doc.images?.length??0,0,'Image textures require an explicit conversion path');
assert.ok(!doc.extensionsUsed?.length,'Unsupported GLB extension');
const shapes={SCALAR:1,VEC2:2,VEC3:3,VEC4:4};
const types={5121:[1,'getUint8'],5123:[2,'getUint16'],5125:[4,'getUint32'],5126:[4,'getFloat32']};
function accessor(index,sourceDoc=doc,sourceBinary=binary){
 const a=sourceDoc.accessors[index],v=sourceDoc.bufferViews[a.bufferView],width=shapes[a.type],format=types[a.componentType];
 assert.ok(a&&!a.sparse&&v.buffer===0&&format);
 const [size,method]=format,stride=v.byteStride??size*width,offset=(v.byteOffset??0)+(a.byteOffset??0);
 const view=new DataView(sourceBinary.buffer,sourceBinary.byteOffset,sourceBinary.byteLength);
 const values=a.componentType===5126?new Float32Array(a.count*width):new Uint32Array(a.count*width);
 for(let row=0;row<a.count;row++)for(let col=0;col<width;col++)values[row*width+col]=view[method](offset+row*stride+col*size,true);
 return values;
}
const topSystems={VH_F_integumentary_system:'integumentary',VH_F_nervous_system:'nervous',VH_F_muscular_system:'muscular',VH_F_reproductive_system:'reproductive',VH_F_digestive_system:'digestive',VH_F_urinary_system:'urinary',VH_F_circulatory_system:'cardiac',VH_F_respiratory_system:'respiratory',VH_F_lymphatic_system:'lymphatic',VH_F_skeletal_system:'skeletal',VH_F_lower_limb:'skeletal'};
const palette={integumentary:[196,155,132],nervous:[218,184,104],sensory:[187,165,143],muscular:[165,84,79],reproductive:[173,116,117],digestive:[169,108,98],urinary:[177,124,104],cardiac:[170,77,77],arterial:[190,79,77],venous:[105,116,167],respiratory:[167,139,149],lymphatic:[136,161,128],skeletal:[217,207,183],connective:[181,188,178],endocrine:[183,130,142],pregnancy:[184,131,130]};
function systemFor(node,branch){
 const name=node.name.toLowerCase(),label=(node.extras?.label??'').toLowerCase(),text=name+' '+label;
 if(/placenta|umbilical.cord|chorion|decidua/.test(text))return 'pregnancy';
 if(branch==='cardiac')return /vein|venous|vena|sinus/.test(text)?'venous':/arter|aorta/.test(text)?'arterial':'cardiac';
 if(branch==='nervous'&&/eye|ear|retina|lens|optic/.test(text))return 'sensory';
 if(branch==='skeletal'&&/cartilage|ligament|meniscus|enthesis|perichondr/.test(text))return 'connective';
 if(branch==='digestive'&&/gland|pancreas/.test(text))return /pancreas/.test(text)?'digestive':'endocrine';
 return branch;
}
const localMatrix=node=>node.matrix?new Matrix4().fromArray(node.matrix):new Matrix4().compose(new Vector3().fromArray(node.translation??[0,0,0]),new Quaternion().fromArray(node.rotation??[0,0,0,1]),new Vector3().fromArray(node.scale??[1,1,1]));
const displayName=node=>{const label=node.extras?.label;const raw=label&&label!=='-'?label.replace(/\s+/g,' ').trim():node.name.replace(/^VH_F_/,'').replaceAll('_',' ');const side=/^(left|right) (.+)$/i.exec(raw);const name=side?`${side[2]} (${side[1].toLowerCase()})`:raw;return name.replace(/^./,letter=>letter.toUpperCase());};
const chunks=[],parts=[],materials={},meshParts=new Map(),nodeParents=new Map();
let segments=[],bytes=0,sourceTriangles=0;
function append(values){const padding=(4-bytes%4)%4;if(padding){segments.push(Buffer.alloc(padding));bytes+=padding;}const offset=bytes;segments.push(Buffer.from(values.buffer,values.byteOffset,values.byteLength));bytes+=values.byteLength;return offset;}
async function flush(){if(!bytes)return;const filename=`hra-female-${chunks.length}.bin`,buffer=Buffer.concat(segments),gz=gzipSync(buffer,{level:9,mtime:0});await fs.writeFile(path.join(out,filename),buffer);await fs.writeFile(path.join(out,filename+'.gz'),gz);chunks.push({url:'/models/'+filename,bytes:buffer.length,gzip:'/models/'+filename+'.gz',gzipBytes:gz.length});segments=[];bytes=0;}
const mpos=new Vector3(),mnorm=new Vector3();
async function visit(index,parent,branch,parentIndex){
 const node=doc.nodes[index],world=parent.clone().multiply(localMatrix(node));nodeParents.set(index,parentIndex);
 const active=index===0?null:parentIndex===0?topSystems[node.name]:branch;
 if(index!==0)assert.ok(active,`Unknown source branch: ${node.name}`);
 if(node.mesh!==undefined){
  const id='HRA:'+node.name,system=systemFor(node,active),normalMatrix=new Matrix3().getNormalMatrix(world),mirror=world.determinant()<0;
  const name=displayName(node);
  const mesh=doc.meshes[node.mesh];
  for(let pi=0;pi<mesh.primitives.length;pi++){
   const primitive=mesh.primitives[pi];assert.equal(primitive.mode??4,4);
   const sourcePos=accessor(primitive.attributes.POSITION),sourceNorm=accessor(primitive.attributes.NORMAL),sourceIndices=accessor(primitive.indices);
   assert.equal(sourcePos.length,sourceNorm.length);
   const positions=new Float32Array(sourcePos.length),normalsFloat=new Float32Array(sourceNorm.length);
   for(let v=0;v<sourcePos.length;v+=3){
    mpos.fromArray(sourcePos,v).applyMatrix4(world);positions[v]=mpos.x;positions[v+1]=mpos.y+.794760942;positions[v+2]=mpos.z;
    mnorm.fromArray(sourceNorm,v).applyMatrix3(normalMatrix).normalize();normalsFloat[v]=mnorm.x;normalsFloat[v+1]=mnorm.y;normalsFloat[v+2]=mnorm.z;
   }
   const original=new Uint32Array(sourceIndices);
   if(mirror)for(let i=0;i<original.length;i+=3){const temp=original[i+1];original[i+1]=original[i+2];original[i+2]=temp;}
   sourceTriangles+=original.length/3;
   const target=Math.max(96,Math.floor(original.length*.35/3)*3);
   const optimized=target<original.length?MeshoptSimplifier.simplify(original,positions,3,target,.002)[0]:original;
   const [remap,count]=MeshoptSimplifier.compactMesh(optimized);
   const compactPositions=new Float32Array(count*3),compactNormals=new Int16Array(count*3),bounds=[[Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity]];
   for(let v=0;v<remap.length;v++)if(remap[v]!==0xffffffff){const to=remap[v];for(let axis=0;axis<3;axis++){const value=positions[v*3+axis];compactPositions[to*3+axis]=value;compactNormals[to*3+axis]=Math.round(Math.max(-1,Math.min(1,normalsFloat[v*3+axis]))*32767);bounds[0][axis]=Math.min(bounds[0][axis],value);bounds[1][axis]=Math.max(bounds[1][axis],value);}}
   if(bytes>6_000_000)await flush();
   const materialId=`HRA:${system}:${primitive.material??'none'}`;
   const factor=doc.materials?.[primitive.material]?.pbrMetallicRoughness?.baseColorFactor??[.5,.5,.5,1];
   const srgb=value=>Math.round(255*(value<=.0031308?12.92*value:1.055*value**(1/2.4)-.055));
   materials[materialId]={color:palette[system],sourceColor:factor.slice(0,3).map(srgb)};
   const partId=mesh.primitives.length===1?id:id+':'+pi;
   const part={id:partId,conceptId:id,name:mesh.primitives.length===1?name:`${name} · surface ${pi+1}`,system,chunk:chunks.length,positions:append(compactPositions),normals:append(compactNormals),indices:append(optimized),vertexCount:count,indexCount:optimized.length,bounds,material:materialId,sourceId:node.name,provenance:{label:'Human Reference Atlas united-female v1.10',url:'https://purl.humanatlas.io/ref-organ/united-female/v1.10',detail:`Original GLB node ${index}; ${node.extras?.ontologyid??'no ontology ID'}; CC BY 4.0.`}};
   parts.push(part);const entries=meshParts.get(index)??[];entries.push(partId);meshParts.set(index,entries);
  }
 }
 for(const child of node.children??[])await visit(child,world,active,index);
}
await visit(doc.scenes[doc.scene??0].nodes[0],new Matrix4(),null,-1);
await flush();
const concepts=[];
function descendants(index){const node=doc.nodes[index],members=[...(meshParts.get(index)??[])];for(const child of node.children??[])members.push(...descendants(child));if(members.length)concepts.push({id:'HRA:'+node.name,name:displayName(node),elements:members});return members;}
descendants(doc.scenes[doc.scene??0].nodes[0]);
// v1.10 omits four pubis and four ischium bone surfaces that v1.5 supplied.
// The two editions share identical ilium world bounds, so these source meshes
// use the same world transform and stage offset without a fitted alignment.
const missingPelvis=/^VH_F_(pubis|ischium)_(spongy|compact)_bone_[LR]$/;
const supplementParts=[];
let supplementSourceTriangles=0;
async function visitSupplement(index,parent){
 const node=supplementDoc.nodes[index],world=parent.clone().multiply(localMatrix(node));
 const match=missingPelvis.exec(node.name);
 if(match&&node.mesh!==undefined){
  assert.ok(!parts.some(p=>p.sourceId===node.name),`Already present in v1.10: ${node.name}`);
  const mesh=supplementDoc.meshes[node.mesh];
  assert.equal(mesh.primitives.length,1,`Unexpected surfaces: ${node.name}`);
  const primitive=mesh.primitives[0];
  const sourcePos=accessor(primitive.attributes.POSITION,supplementDoc,supplementBinary);
  const sourceNorm=accessor(primitive.attributes.NORMAL,supplementDoc,supplementBinary);
  const sourceIndices=accessor(primitive.indices,supplementDoc,supplementBinary);
  const normalMatrix=new Matrix3().getNormalMatrix(world),mirrored=world.determinant()<0;
  const positions=new Float32Array(sourcePos.length),normalsFloat=new Float32Array(sourceNorm.length);
  for(let v=0;v<sourcePos.length;v+=3){
   mpos.fromArray(sourcePos,v).applyMatrix4(world);positions[v]=mpos.x;positions[v+1]=mpos.y+.794760942;positions[v+2]=mpos.z;
   mnorm.fromArray(sourceNorm,v).applyMatrix3(normalMatrix).normalize();normalsFloat[v]=mnorm.x;normalsFloat[v+1]=mnorm.y;normalsFloat[v+2]=mnorm.z;
  }
  const original=new Uint32Array(sourceIndices);
  if(mirrored)for(let i=0;i<original.length;i+=3){const temp=original[i+1];original[i+1]=original[i+2];original[i+2]=temp;}
  supplementSourceTriangles+=original.length/3;
  const target=Math.max(96,Math.floor(original.length*.35/3)*3);
  const optimized=target<original.length?MeshoptSimplifier.simplify(original,positions,3,target,.002)[0]:original;
  const [remap,count]=MeshoptSimplifier.compactMesh(optimized);
  const compactPositions=new Float32Array(count*3),compactNormals=new Int16Array(count*3),bounds=[[Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity]];
  for(let v=0;v<remap.length;v++)if(remap[v]!==0xffffffff){const to=remap[v];for(let axis=0;axis<3;axis++){const value=positions[v*3+axis];compactPositions[to*3+axis]=value;compactNormals[to*3+axis]=Math.round(Math.max(-1,Math.min(1,normalsFloat[v*3+axis]))*32767);bounds[0][axis]=Math.min(bounds[0][axis],value);bounds[1][axis]=Math.max(bounds[1][axis],value);}}
  if(bytes>6_000_000)await flush();
  const materialId=`HRA15:skeletal:${primitive.material??'none'}`;
  const factor=supplementDoc.materials?.[primitive.material]?.pbrMetallicRoughness?.baseColorFactor??[.5,.5,.5,1];
  const srgb=value=>Math.round(255*(value<=.0031308?12.92*value:1.055*value**(1/2.4)-.055));
  materials[materialId]={color:palette.skeletal,sourceColor:factor.slice(0,3).map(srgb)};
  const id='HRA15:'+node.name;
  const name=`${match[1]==='pubis'?'Pubis':'Ischium'}, ${match[2]==='spongy'?'spongy':'compact'} bone (${match[3]==='L'?'left':'right'})`;
  const part={id,conceptId:id,name,system:'skeletal',chunk:chunks.length,positions:append(compactPositions),normals:append(compactNormals),indices:append(optimized),vertexCount:count,indexCount:optimized.length,bounds,material:materialId,sourceId:node.name,license:'CC BY 4.0',provenance:{label:'Human Reference Atlas united-female v1.5 pelvic supplement',url:'https://purl.humanatlas.io/ref-organ/united-female/v1.5',detail:`Original GLB node ${index}; ${node.extras?.ontologyid??'no ontology ID'}; CC BY 4.0. Omitted from v1.10.`}};
  parts.push(part);supplementParts.push(part);concepts.push({id,name,elements:[id]});
 }
 for(const child of node.children??[])await visitSupplement(child,world);
}
for(const rootIndex of supplementDoc.scenes[supplementDoc.scene??0].nodes)await visitSupplement(rootIndex,new Matrix4());
assert.equal(supplementParts.length,8,'Expected four pubis and four ischium surfaces from v1.5');
await flush();
for(const bone of ['pubis','ischium']){
 const members=supplementParts.filter(p=>p.sourceId.startsWith(`VH_F_${bone}_`)).map(p=>p.id);
 concepts.push({id:`HRA15:VH_F_${bone}`,name:bone==='pubis'?'Pubis':'Ischium',elements:members});
}
for(const id of ['HRA:VH_F_pelvis','HRA:VH_F_skeletal_system','HRA:VH_F']){
 const concept=concepts.find(c=>c.id===id);
 assert.ok(concept,`Missing v1.10 pelvis ancestor ${id}`);
 concept.elements.push(...supplementParts.map(p=>p.id));
}
assert.equal(new Set(parts.map(p=>p.id)).size,parts.length);
assert.equal(new Set(concepts.map(c=>c.id)).size,concepts.length);
const atlas={version:'Human Reference Atlas united-female v1.10 + v1.5 pelvic supplement',sex:'female',source:'Human Reference Atlas',scope:'Original united-female references; partial muscle coverage.',parts,concepts,chunks,materials,triangles:parts.reduce((sum,p)=>sum+p.indexCount/3,0),sourceTriangles,quality:'optimized',provenance:{sourceUrl:SOURCE_URL,sourceSha256:SOURCE_SHA256,sourceFile:'3d-vh-f-united.glb',supplementUrl:SUPPLEMENT_URL,supplementSha256:SUPPLEMENT_SHA256,supplementParts:supplementParts.map(p=>p.id),supplementSourceTriangles,license:'CC BY 4.0',notes:'Direct HRA GLB imports; v1.5 supplies eight pelvic bone surfaces omitted from v1.10. Original node hierarchy and material colors retained as sourceColor. Neither GLB has image textures.'}};
await fs.writeFile(path.join(out,'atlas-hra-female.json'),JSON.stringify(atlas));
const split=spawn(process.execPath,[path.join(ROOT,'scripts','split-female-embryo.mjs'),out],{stdio:'inherit'});
const splitCode=await new Promise((resolve,reject)=>{split.once('error',reject);split.once('close',resolve);});
if(splitCode!==0)throw new Error(`Female/embryo split failed (exit ${splitCode})`);
console.log(JSON.stringify({parts:parts.length,concepts:concepts.length,sourceTriangles,triangles:atlas.triangles,chunks:chunks.length},null,2));
