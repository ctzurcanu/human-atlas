/** Convert licensed anatomical geometry, never the source viewer software.
 * node scripts/import-male-detail.mjs [CACHE_DIRECTORY]
 * Each upstream GLB is pinned by SHA-256 in data/male-source-index.json.
 */
import fs from 'node:fs/promises';
import {finalizeMaleCatalogue} from './male-catalogue.mjs';
import path from 'node:path';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';
import {MeshoptDecoder,MeshoptSimplifier} from 'meshoptimizer';
import {Matrix4,Vector3,Matrix3,Quaternion} from 'three';

await Promise.all([MeshoptDecoder.ready,MeshoptSimplifier.ready]);
const readJson=async file=>JSON.parse(await fs.readFile(new URL(file,import.meta.url),'utf8'));
const index=await readJson('./data/male-source-index.json');
const catalogue=await readJson('./data/male-source-catalog.json');
const aliases=await readJson('./data/male-source-aliases.json');
const representations=await readJson('./data/male-representation-groups.json');
const synonyms=await readJson('./data/male-representation-aliases.json');
const resolve=(name,region)=>aliases[`${region}|${name}`]??(catalogue[name]?name:normalizedNames.get(normalize(name))); 
const cache=process.argv[2]??'/tmp/male-atlas-source';
const full=process.argv.includes('--full');
const prefix=full?'male-full':'male-detail';
const out=new URL('../public/models/',import.meta.url);
await fs.mkdir(cache,{recursive:true});
function displayName(record){
 const label=record.label.replace(/\.l$/,' (left)').replace(/\.r$/,' (right)');
 return record.side&&!/\b(left|right)\b/i.test(label)?`${label} (${record.side})`:label;
}
const tissueSystem={bone:'skeletal',cartilage:'connective',ligament:'connective',capsule:'connective',synovium:'connective',fascia:'fascia',muscle:'muscular',tendon:'muscular',attachment:'attachments',region:'regions',surface:'integumentary',artery:'arterial',vein:'venous',nerve:'nervous',cns:'nervous',lymphoid:'lymphatic',sense:'sensory',heart:'cardiac'};
function classify(name,record){
 if(record.schematic)return 'schematic';
 if(tissueSystem[record.tissue])return tissueSystem[record.tissue];
 const groups=record.groups.join(' ');
 if(/Urinary/.test(groups))return 'urinary';
 if(/Genital/.test(groups))return 'reproductive';
 if(/Endocrine/.test(groups))return 'endocrine';
 if(/Respiratory/.test(groups))return 'respiratory';
 if(/Digestive/.test(groups))return 'digestive';
 if(/Peritoneal/.test(groups)||name==='Meso-appendix')return 'connective';
 if(/nasal cavity|Nasopharynx/.test(name))return 'respiratory';
 if(/Oropharynx|Laryngopharynx|Soft palate/.test(name))return 'digestive';
 if(/kidney|ureter|bladder|urethra/i.test(name))return 'urinary';
 if(/thyroid|adrenal|pituitary|pineal/i.test(name))return 'endocrine';
 if(/testis|prostate|seminal|penis|epididymis|deferens/i.test(name))return 'reproductive';
 throw new Error(`Unmapped anatomy: ${name} (${record.tissue})`);
}
const parts=[],chunks=[],sourceGroups=new Map(),omitted=[];
let segments=[],bytes=0,sourceTriangles=0;
const used=new Set(),unidentifiedGeometry=[];
const normalize=name=>name.replace(/\s+/g,' ').replace(/\.\./g,'.').replace(/M(\d)[ -]segment/g,'M$1-segment').replace(/\(M2\)/g,'(M2-segment)').replace(/'$/,'');
const normalizedNames=new Map();
for(const name of Object.keys(catalogue)){
 const key=normalize(name);
 normalizedNames.set(key,normalizedNames.has(key)?null:name);
}
function append(values){
 const padding=(4-bytes%4)%4;if(padding){segments.push(Buffer.alloc(padding));bytes+=padding;}
 const offset=bytes;segments.push(Buffer.from(values.buffer,values.byteOffset,values.byteLength));bytes+=values.byteLength;return offset;
}
async function flush(){
 if(!bytes)return;
 const name=`${prefix}-${chunks.length}.bin`,buffer=Buffer.concat(segments),gzip=gzipSync(buffer,{level:9});
 await fs.writeFile(new URL(name,out),buffer);await fs.writeFile(new URL(name+'.gz',out),gzip);
 chunks.push({url:'/models/'+name,bytes,gzip:'/models/'+name+'.gz',gzipBytes:gzip.length});segments=[];bytes=0;
}
const sources=[];
for(const region of index.regions){
 const file=path.join(cache,region.region+'.glb');let glb;
 try{glb=await fs.readFile(file);}catch(error){if(error.code!=='ENOENT')throw error;
  const response=await fetch(`https://anatomy-atlas.brianp.chatgpt.site/geometry/${region.region}.glb?v=${region.geometry_sha256}`);
  if(!response.ok)throw new Error(`${region.region}: HTTP ${response.status}`);
  glb=Buffer.from(await response.arrayBuffer());await fs.writeFile(file,glb);
 }
 assert.equal(createHash('sha256').update(glb).digest('hex'),region.geometry_sha256,`Source revision changed: ${region.region}`);
 assert.equal(glb.readUInt32LE(0),0x46546c67);assert.equal(glb.readUInt32LE(4),2);
 const jsonLength=glb.readUInt32LE(12),doc=JSON.parse(glb.subarray(20,20+jsonLength).toString()),binary=glb.subarray(28+jsonLength);
 sources.push({region,doc,binary});
}
// The catalog has stale region assignments for some attachment markers.
// Prefer its stated chunk when actually present, otherwise the first real source.
const ownerRegions=new Map();
for(const {region,doc} of sources)for(const node of doc.nodes){
 if(!node.name)continue;
 const owner=resolve(node.name,region.region);
 if(!owner)continue;
 if(!ownerRegions.has(owner)||catalogue[owner].chunk===region.region)ownerRegions.set(owner,region.region);
}
for(const {region,doc,binary} of sources){
 const decoded=new Map();
 function view(i){
  if(decoded.has(i))return decoded.get(i);
  const v=doc.bufferViews[i],extension=v.extensions?.EXT_meshopt_compression;let buffer;
  if(extension){
   assert.equal(extension.buffer,0);buffer=new Uint8Array(extension.count*extension.byteStride);
   MeshoptDecoder.decodeGltfBuffer(buffer,extension.count,extension.byteStride,binary.subarray(extension.byteOffset??0,(extension.byteOffset??0)+extension.byteLength),extension.mode,extension.filter);
  }else{assert.equal(v.buffer,0);buffer=binary.subarray(v.byteOffset??0,(v.byteOffset??0)+v.byteLength);}
  decoded.set(i,buffer);return buffer;
 }
 const widths={SCALAR:1,VEC2:2,VEC3:3,VEC4:4};
 const formats={5120:[1,'getInt8',127],5121:[1,'getUint8',255],5122:[2,'getInt16',32767],5123:[2,'getUint16',65535],5125:[4,'getUint32',4294967295],5126:[4,'getFloat32',1]};
 function accessor(i){
  const a=doc.accessors[i];assert.ok(!a.sparse,'Sparse accessors need explicit support');
  const v=doc.bufferViews[a.bufferView],b=view(a.bufferView),d=new DataView(b.buffer,b.byteOffset,b.byteLength);
  const [size,method,divisor]=formats[a.componentType],width=widths[a.type],stride=v.byteStride??width*size;
  const values=new Float64Array(a.count*width);
  for(let row=0;row<a.count;row++)for(let col=0;col<width;col++){
   const value=d[method]((a.byteOffset??0)+row*stride+col*size,true);
   values[row*width+col]=a.normalized?Math.max(-1,value/divisor):value;
  }
  return values;
 }
 const geometry=new Map();
 function walk(i,parent,owner){
  const node=doc.nodes[i],local=node.matrix?new Matrix4().fromArray(node.matrix):new Matrix4().compose(new Vector3().fromArray(node.translation??[0,0,0]),new Quaternion().fromArray(node.rotation??[0,0,0,1]),new Vector3().fromArray(node.scale??[1,1,1]));
  const world=parent.clone().multiply(local);if(node.name)owner=resolve(node.name,region.region)??null;
  if(node.mesh!==undefined){
   if(!owner){unidentifiedGeometry.push({region:region.region,node:node.name??i});return;}
   const record=catalogue[owner];
   if(ownerRegions.get(owner)!==region.region){unidentifiedGeometry.push({region:region.region,node:owner,reason:'Another chunk supplies this structure; duplicate representation excluded'});return;}
   assert.ok(['CC BY-SA 4.0','CC BY-NC 4.0','CC BY-NC-SA 4.0'].includes(record.licence));
   const normalMatrix=new Matrix3().getNormalMatrix(world),p=new Vector3(),n=new Vector3();
   for(const primitive of doc.meshes[node.mesh].primitives){
    const materialName=doc.materials?.[primitive.material]?.name;
    const surfaceOwner=record.split?.[materialName]??owner;
    assert.ok(catalogue[surfaceOwner],`Unknown material-defined anatomy: ${surfaceOwner}`);
    let result=geometry.get(surfaceOwner);if(!result){result={positions:[],normals:[],indices:[]};geometry.set(surfaceOwner,result);}

    assert.equal(primitive.mode??4,4);const pos=accessor(primitive.attributes.POSITION),norm=accessor(primitive.attributes.NORMAL),ind=accessor(primitive.indices),base=result.positions.length/3;
    for(let v=0;v<pos.length;v+=3){
     p.fromArray(pos,v).applyMatrix4(world);const correctedReflection=surfaceOwner==='Sternocleidomastoid muscle.o1l'&&catalogue[surfaceOwner].reflected_in_x;result.positions.push(correctedReflection?-p.x:p.x,p.y+.031,p.z);
     n.fromArray(norm,v).applyMatrix3(normalMatrix).normalize();result.normals.push(correctedReflection?-n.x:n.x,n.y,n.z);
    }
    for(let v=0;v<ind.length;v+=3){const a=base+ind[v],b=base+ind[v+1],c=base+ind[v+2];result.indices.push(a,...((world.determinant()<0)!==(surfaceOwner==='Sternocleidomastoid muscle.o1l'&&!!catalogue[surfaceOwner].reflected_in_x)?[c,b]:[b,c]));}
   }
  }
  for(const child of node.children??[])walk(child,world,owner);
 }
 for(const root of doc.scenes[doc.scene??0].nodes)walk(root,new Matrix4(),null);
 for(const [name,g] of geometry){
  assert.ok(!used.has(name),`Duplicate structure ${name}`);used.add(name);
  const record=catalogue[name],system=classify(name,record);
  const pos=new Float32Array(g.positions),indices=new Uint32Array(g.indices);sourceTriangles+=indices.length/3;
  const target=Math.min(indices.length,Math.max(96,Math.floor(indices.length*.3/3)*3));
  const simplified=full?indices:MeshoptSimplifier.simplify(indices,pos,3,target,.002)[0];
  const [remap,count]=MeshoptSimplifier.compactMesh(simplified),positions=new Float32Array(count*3),normals=new Int16Array(count*3);
  const bounds=[[Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity]];
  for(let v=0;v<remap.length;v++)if(remap[v]!==0xffffffff)for(let axis=0;axis<3;axis++){
   const value=pos[v*3+axis];positions[remap[v]*3+axis]=value;normals[remap[v]*3+axis]=Math.round(Math.max(-1,Math.min(1,g.normals[v*3+axis]))*32767);
   bounds[0][axis]=Math.min(bounds[0][axis],value);bounds[1][axis]=Math.max(bounds[1][axis],value);
  }
  if(bytes>4_000_000)await flush();
  const id='DETAIL:'+name;
  parts.push({id,conceptId:id,name:displayName(record),system,
   chunk:chunks.length,positions:append(positions),normals:append(normals),indices:append(simplified),vertexCount:count,indexCount:simplified.length,bounds,
   regions:record.view_memberships??[],depth:record.depth??.5,groups:record.groups,tissue:record.tissue,sourceId:name,source:record.source,license:record.licence,schematic:!!record.schematic,mirrored:record.mirrored??null,
   provenance:{label:record.source==='o3m'?'Open 3D Model':record.source==='derived'?'Derived body surface':'Z-Anatomy',
    url:'/ATTRIBUTION.md',detail:`${record.schematic?'Schematic reference; source placement is not anatomical. ':''}${record.licence}${record.noncommercial?' · noncommercial component':''}. ${record.mirrored?'Mirrored from the opposite limb. ':''}${full?'Full-resolution geometry.':'Optimized geometry.'}`}});
  for(const group of record.groups){if(!sourceGroups.has(group))sourceGroups.set(group,[]);sourceGroups.get(group).push(id);}
 }
 console.log(`${region.region}: ${geometry.size} named meshes converted`);
}
await flush();
for(const name of Object.keys(catalogue))if(!used.has(name))omitted.push(name);
const concepts=parts.map(p=>({id:p.id,name:p.name,elements:[p.id]}));
for(const [name,elements] of sourceGroups)concepts.push({id:'DETAIL:group:'+name,name,elements});
// Canonical whole-structure handles and factual source-name aliases.
for(const [name,members] of Object.entries(representations)){
 const elements=members.flatMap(n=>parts.filter(p=>p.sourceId===n).map(p=>p.id));
 if(elements.length){const existing=concepts.find(c=>c.id==='DETAIL:'+name);if(existing)existing.elements=elements;else concepts.push({id:'DETAIL:'+name,name:displayName({label:name}),elements});}
}
for(const [name,target] of Object.entries(synonyms)){
 const c=concepts.find(c=>c.id==='DETAIL:'+target);if(c&&!concepts.some(c=>c.id==='DETAIL:'+name))concepts.push({id:'DETAIL:'+name,name:displayName({label:name}),elements:c.elements});
}
// Major-organ shortcuts remain useful even when the source models chambers separately.
for(const [name,system] of [['Heart','cardiac'],['Brain','nervous']]){
 if(concepts.some(c=>c.name.toLowerCase()===name.toLowerCase()))continue;
 const elements=parts.filter(p=>name==='Heart'?p.system===system:p.tissue==='cns'&&!p.schematic&&!/spinal|cauda/i.test(p.name)).map(p=>p.id);
 if(elements.length)concepts.push({id:'DETAIL:organ:'+name,name,elements});
}
const atlas={version:'Z-Anatomy + Open 3D Model · Brian Pridgen adaptation',quality:full?'full':'optimized',sex:'male',source:'Z-Anatomy + Open 3D Model',scope:'Detailed male reference with upper-limb integration and mirrored anatomy; mixed component licenses.',parts,concepts,chunks,
 triangles:parts.reduce((n,p)=>n+p.indexCount/3,0),sourceTriangles,optimized:{maximumRelativeError:full?0:.002,preservedMeshes:parts.length},
 provenance:{sourceUrl:'https://anatomy-atlas.brianp.chatgpt.site/',notice:'/ATTRIBUTION.md',sourceHashes:index.regions.map(r=>({region:r.region,sha256:r.geometry_sha256})),catalogueEntriesWithoutGeometry:omitted,unidentifiedGeometry}};
finalizeMaleCatalogue(atlas,catalogue,representations,synonyms);
await fs.writeFile(new URL(`atlas-${prefix}.json`,out),JSON.stringify(atlas));
console.log(JSON.stringify({parts:parts.length,concepts:concepts.length,triangles:atlas.triangles,omittedCatalogueEntries:omitted.length,downloadMB:chunks.reduce((n,c)=>n+c.gzipBytes,0)/1e6},null,2));
