/** Build a development-only atlas from the Anatomy Atlas site's CC-licensed GLBs.
 * Raw downloads, textures, and converted geometry remain in .local-models/.
 * Run: npm run import:reference-glb
 */
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import fs from 'node:fs/promises';
import path from 'node:path';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';
import {gzipSync} from 'node:zlib';
import {Matrix3,Matrix4,Quaternion,Vector3} from 'three';
import {MeshoptDecoder} from 'meshoptimizer';
import {placeSpinalCord,spinalCordTissues} from './reference-spinal-cord.mjs';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const local=path.join(root,'.local-models');
const sourceDir=path.join(local,'reference','source');
const outputDir=path.join(local,'reference');
const textureDir=path.join(outputDir,'textures');
const origin='https://anatomy-atlas.brianp.chatgpt.site';
const files=[
 ['body-surface','2b5c229a65a3581b213eea39c46937045b52d0cfc79fbf62815e2a1ce6d75f7f'],
 ['core','9766b9c38417f239b5bf70a4ac7118731a871d83f7bb8e84b2c99f4f0b52917c'],
 ['thorax','3d6f0c69d8c5bca7d47985acebd1a151793d08e7d4d4f8cd2bad0215b0ef572f'],
 // This endpoint served a different valid GLB than its query-string version when downloaded.
 ['abdomen-pelvis','c086adfb4694bff4d6c2112c62ff26e62f7826f12ca239ff268f48d7f11f638b','82e88907a779120a4edc8641944a9f57ce28cd14ffa98d380a7848efc1d1b245'],
 ['back','56b1d070bbdaad2c37d48bd21cb7df1d67920a8f0b8f5490a5abe26c9636bd92'],
 ['neck','62edc42bb946dba23ae5e8f0b7421d5e85541865228cea00e59653893e8c8dbc'],
 ['upper-limb-right','bcd77622c467aa7557b1312f760fb2b6c5da12ff110d21761050c2f73b7d38f1'],
 ['upper-limb-left','bc14c0ba914f310d82594f3e96894404b193453f9581d0e2e77a849a09a80ba1'],
 ['lower-limb-right','bbc7561f128b11f9c4703375309b9051b6c1c3247fcef82025834def31d6fbc7'],
 ['lower-limb-left','50ce8a1ef4cc324eaf36bbc43154341b27f3bb2a150fa55fefaabaf93a5c2574'],
 ['head','0b4b93cea5542470222e98d675b9d2497b9113267264a953002045899d7028cb'],
];
const shapeWidth={SCALAR:1,VEC2:2,VEC3:3,VEC4:4};
const component={5120:{bytes:1,read:'getInt8',signed:true},5121:{bytes:1,read:'getUint8'},5122:{bytes:2,read:'getInt16',signed:true},5123:{bytes:2,read:'getUint16'},5125:{bytes:4,read:'getUint32'},5126:{bytes:4,read:'getFloat32'}};
const normalizeName=name=>name.toLowerCase().replace(/[^a-z0-9]/g,'');
const sourceAtlas=JSON.parse(await fs.readFile(path.join(root,'public','models','atlas-male-complete.json'),'utf8'));
const sourceByName=new Map(),sourceByNormalized=new Map();
for(const part of sourceAtlas.parts){
 if(part.sourceId){sourceByName.set(part.sourceId,part);sourceByNormalized.set(normalizeName(part.sourceId),part);}
}
const sourceFor=name=>sourceByName.get(name)??sourceByNormalized.get(normalizeName(name));
const title=name=>name.replace(/\.l$/,' (left)').replace(/\.r$/,' (right)').replace(/^./,letter=>letter.toUpperCase());
const srgb=value=>Math.round(255*(value<=.0031308?12.92*value:1.055*value**(1/2.4)-.055));
const palette={integumentary:[205,185,171],nervous:[220,189,111],sensory:[190,164,153],muscular:[163,82,78],reproductive:[174,120,118],digestive:[168,106,95],urinary:[177,125,106],cardiac:[170,77,77],arterial:[194,77,74],venous:[105,116,167],respiratory:[169,138,147],lymphatic:[136,161,128],skeletal:[217,207,185],connective:[185,181,166],endocrine:[185,130,142],fascia:[191,174,158],attachments:[150,92,83]};
// The reference viewer exposes these source IDs as suppressed aliases. Their
// replacement heads/materials are present elsewhere in the downloaded GLBs.
const suppressed=new Set([
 'Deltoid muscle.r','Deltoid muscle.l','Pectoralis major.r','Pectoralis major.l','Trapezius muscle.r','Trapezius muscle.l','Adductor pollicis.r','Adductor pollicis.l',
 'Coracoclavicular ligament.r','Coracoclavicular ligament.l','Liver','Subacromial-subdeltoid bursa.r','Subacromial-subdeltoid bursa.l',
 'Spinal part of deltoid muscle.r','Spinal part of deltoid muscle.l','Abdominal head of pectoralis major muscle.r','Abdominal head of pectoralis major muscle.l',
 'Flexor digiti minimi brevis of hand.r','Flexor digiti minimi brevis of hand.l','Radial annular ligament.r','Radial annular ligament.l',
 'Oblique cord or radio-ulnar syndesmosis.r','Oblique cord or radio-ulnar syndesmosis.l','Ulnar collateral ligament of wrist.r','Ulnar collateral ligament of wrist.l',
 'Radial collateral ligament of wrist.r','Radial collateral ligament of wrist.l','Triangular fibro cartilage disc.r','Triangular fibro cartilage disc.l',
 'Collateral ligaments of metacarpal joints.r','Collateral ligaments of metacarpal joints.l','Collateral ligaments of metacarpophalangeal joints.r','Collateral ligaments of metacarpophalangeal joints.l',
 'Extensor digitorum - Extensor indicis tendon sheath.r','Extensor digitorum - Extensor indicis tendon sheath.l','Conoid ligament.r','Conoid ligament.l',
 'Trapezoid ligament.r','Trapezoid ligament.l','Radiocapitate ligament.r','Radiocapitate ligament.l',
]);
const licenseFor=name=>/^Kidney\.[lr]$/i.test(name)?'CC BY-NC 4.0':/^(Cochlea|Vestibule)\.[lr]$/i.test(name)?'CC BY-NC-SA 4.0':'CC BY-SA 4.0';
function systemFor(name,material,matched){
 if(/^overlays?(?:\.\d+)?$/i.test(material))return 'attachments';
 if(matched)return matched.system;
 const text=`${name} ${material}`.toLowerCase();
 if(/body surface|\bskin\b|epiderm|dermis/.test(text))return 'integumentary';
 if(/arter|aorta/.test(text))return 'arterial';
 if(/vein|venous|vena cava/.test(text))return 'venous';
 if(/nerve|gangli|plexus|spinal cord|\broot\b/.test(text))return 'nervous';
 if(/cartilage|ligament|tendon|capsule|sheath/.test(text))return 'connective';
 if(/bone|skull|vertebra|rib|humerus|femur/.test(text))return 'skeletal';
 if(/muscle|pectoral|flexor|extensor|adductor|abductor/.test(text))return 'muscular';
 if(/heart|cardiac/.test(text))return 'cardiac';
 if(/lung|bronch|trachea|larynx/.test(text))return 'respiratory';
 if(/kidney|ureter|bladder|urethra/.test(text))return 'urinary';
 if(/lymph|spleen/.test(text))return 'lymphatic';
 if(/eye|ear|cochlea|vestibul/.test(text))return 'sensory';
 if(/ovary|uterus|testis|penis|prostate/.test(text))return 'reproductive';
 return 'digestive';
}
function regionsFor(stem,bounds,name='',matched){
 const y=(bounds[0][1]+bounds[1][1])/2;
 if(stem==='body-surface')return [];
 const groups=(matched?.groups??[]).map(group=>group.toLowerCase());
 const namedRegions=[
  ['upper-left',/^left (?:upper limb|arm|hand|forearm|shoulder)(?: region)?$/],
  ['upper-right',/^right (?:upper limb|arm|hand|forearm|shoulder)(?: region)?$/],
  ['lower-left',/^left (?:lower limb|leg|foot|thigh)(?: region)?$/],
  ['lower-right',/^right (?:lower limb|leg|foot|thigh)(?: region)?$/],
  ['head-neck',/^(?:head|neck|brain|face|skull|cranium|head and neck)$/],
  ['torso',/^(?:trunk|torso|thorax|chest|abdomen|pelvis|regions of trunk)$/],
 ].filter(([,pattern])=>groups.some(group=>pattern.test(group))).map(([region])=>region);
 if(namedRegions.length)return namedRegions;
 if(stem==='head'||stem==='neck')return ['head-neck'];
 if(stem==='thorax'||stem==='abdomen-pelvis'||stem==='back')return ['torso'];
 if(stem==='core'){
  // Core also contains a few long limb neurovascular objects. Their bounds
  // are lateral to the thoracic wall and should not appear detached beside a
  // torso-only cross-section.
  if(bounds[1][0]>.22&&bounds[0][0]>.065)return ['upper-left'];
  if(bounds[0][0]<-.22&&bounds[1][0]<-.065)return ['upper-right'];
  return y>1.4?['head-neck']:['torso'];
 }
 const side=stem.endsWith('right')?'right':'left';
 const limb=stem.startsWith('upper')?'upper':'lower';
 const regions=[`${limb}-${side}`];
 const nearMidline=side==='left'?bounds[0][0]<.20:bounds[1][0]>-.20;
 const upperGirdle=/clavicle|scapula|pectoral|trapezius|serratus|rhomboid|latissimus|thorac|intercost|sternocleidomastoid|levator scapulae|supraspinatus|infraspinatus|subscapularis|teres major|teres minor|subclavian|brachial plexus/i.test(name);
 if(nearMidline&&(limb==='upper'?upperGirdle&&bounds[1][1]>1.12:bounds[1][1]>.72))regions.push('torso');
 return regions;
}
function surfaceRegion(x,y){
 if(y>1.48&&Math.abs(x)<.19)return 'head-neck';
 if(y<.77)return x<0?'lower-right':'lower-left';
 if(Math.abs(x)<.205)return 'torso';
 return x<0?'upper-right':'upper-left';
}
function surfaceRegions(positions,indices){
 const count=indices.length/3,parent=new Int32Array(count),centers=new Float32Array(count*2),band=new Uint8Array(count),owner=new Map();
 for(let i=0;i<count;i++)parent[i]=i;
 const root=i=>{while(parent[i]!==i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
 const join=(a,b)=>{a=root(a);b=root(b);if(a!==b)parent[b]=a;};
 for(let i=0;i<count;i++){
  const a=indices[i*3],b=indices[i*3+1],c=indices[i*3+2];
  const x=(positions[a*3]+positions[b*3]+positions[c*3])/3,y=(positions[a*3+1]+positions[b*3+1]+positions[c*3+1])/3;
  centers[i*2]=x;centers[i*2+1]=y;
  if(y<.77||y>1.15)continue;
  band[i]=1;
  for(const vertex of [a,b,c]){const previous=owner.get(vertex);if(previous!==undefined)join(i,previous);else owner.set(vertex,i);}
 }
 const components=new Map();
 for(let i=0;i<count;i++)if(band[i]){const id=root(i),entry=components.get(id)??{sum:0,count:0};entry.sum+=centers[i*2];entry.count++;components.set(id,entry);}
 assert.equal([...components.values()].filter(({count})=>count>1000).length,3,'The torso and two arms must separate below the axilla');
 return Array.from({length:count},(_,i)=>{
  const x=centers[i*2],y=centers[i*2+1];
  if(!band[i])return surfaceRegion(x,y);
  const group=components.get(root(i)),meanX=group.sum/group.count;
  return meanX>.17?'upper-left':meanX<-.17?'upper-right':'torso';
 });
}
async function download(url,destination){
 try{await fs.access(destination);return;}catch{}
 const temporary=destination+'.download';
 const child=spawn('curl',['--fail','--location','--silent','--show-error','--output',temporary,url],{stdio:'inherit'});
 const code=await new Promise((resolve,reject)=>{child.once('error',reject);child.once('close',resolve);});
 assert.equal(code,0,`Download failed: ${url}`);
 await fs.rename(temporary,destination);
}
await fs.mkdir(sourceDir,{recursive:true});await fs.mkdir(textureDir,{recursive:true});
await Promise.all(files.map(([stem,version])=>download(`${origin}/geometry/${stem}.glb?v=${version}`,path.join(sourceDir,`${stem}.glb`))));
await download(`${origin}/NOTICE.md`,path.join(sourceDir,'NOTICE.md'));
await MeshoptDecoder.ready;

const materials={},parts=[],concepts=new Map(),chunks=[],sourceFiles=[];
const cordExtrusionsSeen=new Set();
let segments=[],byteLength=0,triangles=0,sourceMatches=0,imagesSaved=0;
function append(values){
 const pad=(4-byteLength%4)%4;if(pad){segments.push(Buffer.alloc(pad));byteLength+=pad;}
 const offset=byteLength;segments.push(Buffer.from(values.buffer,values.byteOffset,values.byteLength));byteLength+=values.byteLength;return offset;
}
async function flush(){
 if(!byteLength)return;
 const filename=`reference/chunk-${chunks.length}.bin`,raw=Buffer.concat(segments),compressed=gzipSync(raw,{level:6,mtime:0});
 await fs.writeFile(path.join(local,filename),raw);await fs.writeFile(path.join(local,filename+'.gz'),compressed);
 chunks.push({url:`/local-models/${filename}`,bytes:raw.length,gzip:`/local-models/${filename}.gz`,gzipBytes:compressed.length});
 segments=[];byteLength=0;
}
const point=new Vector3(),normal=new Vector3();
for(const [stem,version,expectedSha=version] of files){
 const glb=await fs.readFile(path.join(sourceDir,`${stem}.glb`));
 assert.equal(glb.readUInt32LE(0),0x46546c67,`${stem}: not GLB`);
 assert.equal(glb.readUInt32LE(4),2,`${stem}: unsupported GLB version`);
 assert.equal(glb.readUInt32LE(8),glb.length,`${stem}: truncated GLB`);
 const jsonLength=glb.readUInt32LE(12),doc=JSON.parse(glb.subarray(20,20+jsonLength).toString());
 assert.equal(glb.readUInt32LE(16),0x4e4f534a,`${stem}: no JSON chunk`);
 const binStart=28+jsonLength;
 assert.equal(glb.readUInt32LE(24+jsonLength),0x004e4942,`${stem}: no binary chunk`);
 const binary=glb.subarray(binStart,binStart+glb.readUInt32LE(20+jsonLength));
 const digest=createHash('sha256').update(glb).digest('hex');
 assert.equal(digest,expectedSha,`${stem}: source file hash differs from the observed GLB`);
 sourceFiles.push({file:`${stem}.glb`,url:`${origin}/geometry/${stem}.glb?v=${version}`,sha256:digest,bytes:glb.length});
 const decodedViews=new Map();
 function viewAt(index){
  if(decodedViews.has(index))return decodedViews.get(index);
  const view=doc.bufferViews[index],compression=view.extensions?.EXT_meshopt_compression;
  let bytes;
  if(compression){
   bytes=new Uint8Array(view.byteLength);
   const encoded=binary.subarray(compression.byteOffset,compression.byteOffset+compression.byteLength);
   MeshoptDecoder.decodeGltfBuffer(bytes,compression.count,compression.byteStride,encoded,compression.mode,compression.filter);
  }else bytes=binary.subarray(view.byteOffset??0,(view.byteOffset??0)+view.byteLength);
  decodedViews.set(index,bytes);return bytes;
 }
 function accessor(index,asIndices=false){
  const a=doc.accessors[index];assert.ok(a&&!a.sparse&&a.bufferView!==undefined,`${stem}: unsupported accessor ${index}`);
  const width=shapeWidth[a.type],format=component[a.componentType],v=doc.bufferViews[a.bufferView];
  assert.ok(width&&format,`${stem}: unsupported accessor type`);
  const bytes=viewAt(a.bufferView),data=new DataView(bytes.buffer,bytes.byteOffset,bytes.byteLength);
  const stride=v.byteStride??format.bytes*width,offset=a.byteOffset??0;
  const out=asIndices?new Uint32Array(a.count*width):new Float32Array(a.count*width);
  for(let row=0;row<a.count;row++)for(let column=0;column<width;column++){
   let value=data[format.read](offset+row*stride+column*format.bytes,true);
   if(!asIndices&&a.normalized&&a.componentType!==5126){
    const maximum=2**(format.bytes*8-(format.signed?1:0))-1;
    value=format.signed?Math.max(-1,value/maximum):value/maximum;
   }
   out[row*width+column]=value;
  }
  return out;
 }
 const imagePaths=new Map();
 async function imagePath(index){
  if(imagePaths.has(index))return imagePaths.get(index);
  const image=doc.images?.[index];assert.ok(image?.bufferView!==undefined,`${stem}: external image unsupported`);
  const ext=image.mimeType==='image/png'?'png':image.mimeType==='image/jpeg'?'jpg':null;
  assert.ok(ext,`${stem}: unsupported image ${image.mimeType}`);
  const filename=`${stem}-${index}.${ext}`,output=path.join(textureDir,filename);
  await fs.writeFile(output,viewAt(image.bufferView));imagesSaved++;
  const url=`/local-models/reference/textures/${filename}`;imagePaths.set(index,url);return url;
 }
 async function materialFor(index,system,matched){
  const key=`REF:${stem}:${index??'default'}:${matched?.material??system}`;
  if(materials[key])return key;
  const m=doc.materials?.[index]??{},pbr=m.pbrMetallicRoughness??{},factor=pbr.baseColorFactor??[1,1,1,1];
  const sourceColor=sourceAtlas.materials?.[matched?.material??'']?.color;
  const fallbackColor=/tendon/i.test(m.name??'')?[209,193,171]:sourceColor??palette[system]??[180,170,160];
  // An explicit glTF baseColorFactor is part of the publisher's material,
  // including on untextured organs. The old importer silently replaced it
  // with the unrelated Z-Anatomy palette, washing out lungs and viscera.
  const record={color:pbr.baseColorFactor?factor.slice(0,3).map(srgb):fallbackColor,roughness:pbr.roughnessFactor??1,metalness:pbr.metallicFactor??0};
  if(m.alphaMode==='BLEND')record.opacity=factor[3];
  if(pbr.baseColorTexture){const t=doc.textures[pbr.baseColorTexture.index];record.map=await imagePath(t.source);}
  if(m.normalTexture){const t=doc.textures[m.normalTexture.index];record.normalMap=await imagePath(t.source);record.normalScale=m.normalTexture.scale??1;}
  materials[key]=record;return key;
 }
 function localMatrix(node){return node.matrix?new Matrix4().fromArray(node.matrix):new Matrix4().compose(new Vector3().fromArray(node.translation??[0,0,0]),new Quaternion().fromArray(node.rotation??[0,0,0,1]),new Vector3().fromArray(node.scale??[1,1,1]));}
 async function visit(index,parent,label){
  const node=doc.nodes[index],world=parent.clone().multiply(localMatrix(node)),name=node.name??label;
  if(node.mesh!==undefined){
   const raw=name??`${stem} mesh ${index}`,matched=sourceFor(raw),display=title(raw),mesh=doc.meshes[node.mesh],normalMatrix=new Matrix3().getNormalMatrix(world),mirrored=world.determinant()<0;
   if(matched)sourceMatches++;
   for(let pi=0;pi<mesh.primitives.length;pi++){
    const primitive=mesh.primitives[pi];assert.equal(primitive.mode??4,4,`${stem}:${raw}: non-triangle primitive`);
    const positionAccessor=doc.accessors[primitive.attributes.POSITION];
    const cordExtrusion=stem==='core'&&positionAccessor.min?.[1]===47397&&positionAccessor.max?.[1]===65535;
    assert.equal(spinalCordTissues.has(raw),cordExtrusion,`${stem}:${raw}: spinal cord extrusion classification changed`);
    if(cordExtrusion)cordExtrusionsSeen.add(raw);
    const sourcePositions=accessor(primitive.attributes.POSITION),sourceNormals=accessor(primitive.attributes.NORMAL),sourceIndices=accessor(primitive.indices,true);
    assert.equal(sourcePositions.length,sourceNormals.length,`${stem}:${raw}: normal count`);
    const vertexCount=sourcePositions.length/3,positions=new Float32Array(sourcePositions.length),normals=new Int16Array(sourceNormals.length),bounds=[[Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity]];
    for(let v=0;v<vertexCount;v++){
     point.fromArray(sourcePositions,v*3).applyMatrix4(world);
     normal.fromArray(sourceNormals,v*3).applyMatrix3(normalMatrix).normalize();
     if(spinalCordTissues.has(raw))placeSpinalCord(point,normal);
     for(let axis=0;axis<3;axis++){
      const value=point.getComponent(axis);positions[v*3+axis]=value;
      normals[v*3+axis]=Math.round(Math.max(-1,Math.min(1,normal.getComponent(axis)))*32767);
      bounds[0][axis]=Math.min(bounds[0][axis],value);bounds[1][axis]=Math.max(bounds[1][axis],value);
     }
    }
    const indices=new Uint32Array(sourceIndices);
    if(mirrored)for(let i=0;i<indices.length;i+=3){const temp=indices[i+1];indices[i+1]=indices[i+2];indices[i+2]=temp;}
    const classified=systemFor(raw,doc.materials?.[primitive.material]?.name??'',matched);
    // All non-derived skin meshes in these GLBs are named topographic region
    // overlays. The reference viewer hides that tissue class by default.
    const system=stem!=='body-surface'&&classified==='integumentary'?'regions':classified;
    const materialId=await materialFor(primitive.material,system,matched),material=doc.materials?.[primitive.material]??{},textureInfo=material.pbrMetallicRoughness?.baseColorTexture??material.normalTexture;
    const transform=textureInfo?.extensions?.KHR_texture_transform;
    let uvs;
    if(primitive.attributes.TEXCOORD_0!==undefined){
     const sourceUvs=accessor(primitive.attributes.TEXCOORD_0);uvs=new Float32Array(sourceUvs.length);
     const offset=transform?.offset??[0,0],scale=transform?.scale??[1,1],rotation=transform?.rotation??0,c=Math.cos(rotation),s=Math.sin(rotation);
     for(let v=0;v<vertexCount;v++){
      const u=sourceUvs[v*2]*scale[0],w=sourceUvs[v*2+1]*scale[1];
      uvs[v*2]=offset[0]+c*u-s*w;uvs[v*2+1]=offset[1]+s*u+c*w;
     }
    }
    let colors;
    if(primitive.attributes.COLOR_0!==undefined){
     const sourceColors=accessor(primitive.attributes.COLOR_0),width=sourceColors.length/vertexCount;
     colors=new Uint8Array(vertexCount*3);
     for(let v=0;v<vertexCount;v++)for(let axis=0;axis<3;axis++)colors[v*3+axis]=Math.round(Math.max(0,Math.min(1,sourceColors[v*width+axis]))*255);
     materials[materialId].vertexColors=true;
    }
    let variants=[{region:null,indices,bounds}];
    // The reference viewer isolates regions of its derived outer surface.
    // Keep the original triangles and normals, but put them in regional pieces
    // so a torso section does not accidentally include hands and legs.
    if(stem==='body-surface'){
     const regional=new Map();
     const labels=surfaceRegions(positions,indices);
     for(let k=0;k<indices.length;k+=3){
      const a=indices[k],b=indices[k+1],c=indices[k+2];
      const region=labels[k/3];
      let entry=regional.get(region);if(!entry){entry={region,indices:[],bounds:[[Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity]]};regional.set(region,entry);}
      entry.indices.push(a,b,c);
      for(const v of [a,b,c])for(let axis=0;axis<3;axis++){
       const value=positions[v*3+axis];entry.bounds[0][axis]=Math.min(entry.bounds[0][axis],value);entry.bounds[1][axis]=Math.max(entry.bounds[1][axis],value);
      }
     }
     variants=[...regional.values()].map(entry=>({...entry,indices:new Uint32Array(entry.indices)}));
    }
    const conceptId=`REF:${raw}`;
    if(!concepts.has(conceptId))concepts.set(conceptId,{id:conceptId,name:display,elements:[]});
    for(const variant of variants){
     if(byteLength>5_000_000)await flush();
     const id=`REF:${stem}:${index}:${pi}${variant.region?`:${variant.region}`:''}`;
     const part={id,conceptId,name:variant.region?`${display} · ${variant.region}`:mesh.primitives.length===1?display:`${display} · ${material.name??`surface ${pi+1}`}`,system,chunk:chunks.length,positions:append(positions),normals:append(normals),indices:append(variant.indices),vertexCount,indexCount:variant.indices.length,bounds:variant.bounds,material:materialId,sourceId:raw,suppressed:suppressed.has(raw),license:licenseFor(raw),provenance:{label:'Anatomy Atlas adapted anatomy',url:`${origin}/NOTICE.md`,detail:`${stem}.glb; ${licenseFor(raw)}; source object ${raw}${spinalCordTissues.has(raw)?'; fitted to the spinal dura centerline':''}.`}};
     if(uvs)part.uvs=append(uvs);if(colors)part.colors=append(colors);
     if(matched){if(matched.groups)part.groups=matched.groups;if(matched.depth!==undefined)part.depth=matched.depth;}
     part.regions=variant.region?[variant.region]:regionsFor(stem,variant.bounds,raw,matched);
     parts.push(part);concepts.get(conceptId).elements.push(id);triangles+=variant.indices.length/3;
    }
   }
  }
  for(const child of node.children??[])await visit(child,world,name);
 }
 for(const node of doc.scenes[doc.scene??0].nodes)await visit(node,new Matrix4(),undefined);
 decodedViews.clear();
 console.log(`${stem}: ${doc.nodes.length} nodes, ${doc.meshes.length} meshes, ${parts.length} total parts`);
}
await flush();
assert.deepEqual([...cordExtrusionsSeen].sort(),[...spinalCordTissues].sort(),'A spinal cord extrusion was missed by the canal correction');
assert.equal(new Set(parts.map(part=>part.id)).size,parts.length,'Duplicate part IDs');
const importedIds=new Set(parts.map(part=>part.sourceId));
for(const name of suppressed)assert.ok(importedIds.has(name),`Suppressed source alias missing: ${name}`);
assert.equal(parts.filter(part=>part.system==='regions').length,696,'Topographic skin overlays changed');
assert.equal(parts.filter(part=>part.system==='integumentary').length,6,'Derived body surface partition changed');
const atlas={version:'Anatomy Atlas adapted GLB local import',sex:'male',source:'Anatomy Atlas / Z-Anatomy / Open 3D Model',scope:'Local educational reference with separately licensed components',quality:'full',parts,concepts:[...concepts.values()],chunks,materials,triangles,provenance:{sourceUrl:origin,noticeUrl:`${origin}/NOTICE.md`,sourceFiles,sourceMatches,imagesSaved,localOnly:true,licenseSummary:'CC BY-SA 4.0; kidney and inner-ear exceptions are NonCommercial (see NOTICE.md).'}};
await fs.writeFile(path.join(local,'reference.json'),JSON.stringify(atlas));
await fs.writeFile(path.join(outputDir,'import-manifest.json'),JSON.stringify({sourceFiles,sourceMatches,parts:parts.length,concepts:concepts.size,triangles,chunks:chunks.length,imagesSaved},null,2));
console.log(JSON.stringify({parts:parts.length,concepts:concepts.size,triangles,chunks:chunks.length,imagesSaved,sourceMatches},null,2));
