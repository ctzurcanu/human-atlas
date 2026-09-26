/** Convert the local Zygote Body mirror into the viewer's atlas format.
 * Source and output geometry stay outside version control and deployment assets.
 * Usage: npm run import:local-models -- [MIRROR_CONTENT_DIRECTORY]
 */
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {gzipSync} from 'node:zlib';
import JSON5 from 'json5';

const source=path.resolve(process.argv[2]??'/Users/christiantzurcanu/Documents/dev/fetch/content');
const out=new URL('../.local-models/',import.meta.url);
const stride=8;
const systems={2:'integumentary',3:'muscular',4:'connective',5:'skeletal',7:'arterial',8:'lymphatic',9:'nervous'};

function label(raw){
 const side=raw.startsWith('l_')?' (left)':raw.startsWith('r_')?' (right)':'';
 const words=raw.replace(/^[lr]_/, '').replaceAll('_',' ').replace(/\s+/g,' ').trim();
 return words.charAt(0).toUpperCase()+words.slice(1)+side;
}

function categories(metadata){
 const dag=new Map(metadata.dag),names=new Map([...metadata.nodes,...metadata.leafs]);
 const result=new Map();
 function visit(id,layer,ancestors){
  const name=names.get(id)??'';
  const route=[...ancestors,name];
  if(metadata.leafs.some(([leaf])=>leaf===id))result.set(name,{layer,route});
  for(const child of dag.get(id)??[])visit(child,layer,route);
 }
 for(const layer of metadata.layers)visit(layer,layer,[]);
 return result;
}

function systemFor(raw,category){
 const {layer,route}=category??{};
 if(layer===7){
  const context=route?.join(' ')??raw;
  if(/heart|cardiac/i.test(context))return 'cardiac';
  return /vein|venous/i.test(context)?'venous':'arterial';
 }
 if(layer===6){
  const context=route?.join(' ')??raw;
  if(/heart|cardiac/i.test(context))return 'cardiac';
  if(/brain|spinal|nerve|gangli/i.test(context))return 'nervous';
  if(/eye|ear|olfact|retina/i.test(context))return 'sensory';
  if(/thyroid|adrenal|pituitary|pineal|thymus/i.test(context))return 'endocrine';
  if(/kidney|ureter|bladder|urethra|urinary/i.test(context))return 'urinary';
  if(/uter|ovar|vagin|testi|prostat|penis|epididym|vas_deferens|genital|reproductive|fallopian/i.test(context))return 'reproductive';
  if(/lung|bronch|trachea|larynx|respirat|nasal/i.test(context))return 'respiratory';
  if(/spleen|lymph/i.test(context))return 'lymphatic';
  return 'digestive';
 }
 if(layer===2&&/eye|iris|pupil|ear|cochlea/i.test(raw))return 'sensory';
 return systems[layer]??'connective';
}

function decode(text,mesh,params){
 const count=mesh.attribRange[1],attributes=new Float32Array(count*stride);
 let offset=mesh.attribRange[0];
 for(let component=0;component<stride;component++){
  let accumulator=0;
  for(let vertex=0;vertex<count;vertex++){
   const value=text.charCodeAt(offset++);
   assert.ok(Number.isFinite(value),'Truncated compressed attributes');
   accumulator+=(value>>>1)^-(value&1);
   attributes[vertex*stride+component]=params.decodeScales[component]*(accumulator+params.decodeOffsets[component]);
  }
 }
 const length=mesh.indexRange[1]*3,indices=new Uint16Array(length);
 let highest=0;
 for(let i=0;i<length;i++){
  const delta=text.charCodeAt(offset++);
  assert.ok(Number.isFinite(delta),'Truncated compressed indices');
  indices[i]=highest-delta;
  if(delta===0)highest++;
  assert.ok(indices[i]<count,`Index ${indices[i]} exceeds ${count} vertices`);
 }
 assert.equal(mesh.lengths.reduce((a,b)=>a+b,0),length,'Named geometry does not cover index range');
 return {attributes,indices};
}

async function importSex(sex){
 const stem=`adult_${sex}2`;
 const manifestDir=path.join(source,'www.zygotebody.com/models',stem);
 const manifestName=(await fs.readdir(manifestDir)).find(name=>name.endsWith('.js'));
 assert.ok(manifestName,`Missing ${sex} model manifest`);
 const manifestText=await fs.readFile(path.join(manifestDir,manifestName),'utf8');
 const manifest=JSON5.parse(manifestText.slice(manifestText.indexOf('{'),manifestText.lastIndexOf('}')+1));
 const modelDir=path.join(source,'cdn.zygotebody.com/models',stem),files=await fs.readdir(modelDir);
 const metadataName=files.find(name=>name.startsWith('entity_metadata')&&name.endsWith('.json'));
 assert.ok(metadataName,`Missing ${sex} entity metadata`);
 const metadata=JSON.parse(await fs.readFile(path.join(modelDir,metadataName),'utf8'));
 const category=categories(metadata),leafNames=new Set(metadata.leafs.map(([,name])=>name));
 const geometries=[];const globalMin=[Infinity,Infinity,Infinity],globalMax=[-Infinity,-Infinity,-Infinity];
 for(const [url,meshes] of Object.entries(manifest.urls)){
  const basename=path.basename(url.split('?')[0],'.utf8');
  const file=files.find(name=>name.startsWith(`${basename}__q_`)&&name.endsWith('.utf8'));
  assert.ok(file,`Missing local mesh ${basename}`);
  const compressed=await fs.readFile(path.join(modelDir,file),'utf8');
  for(const mesh of meshes){
   const decoded=decode(compressed,mesh,manifest.decodeParams);
   geometries.push({mesh,...decoded});
   for(let i=0;i<decoded.attributes.length;i+=stride)for(let axis=0;axis<3;axis++){
    const value=decoded.attributes[i+axis];
    globalMin[axis]=Math.min(globalMin[axis],value);globalMax[axis]=Math.max(globalMax[axis],value);
   }
  }
 }
 const textureDir=path.join(source,'cdn.zygotebody.com/models/Adult_Maps');
 const textureFiles=await fs.readdir(textureDir),materials={};
 let textureIndex=0;
 for(const name of new Set(geometries.map(({mesh})=>mesh.material))){
  const definition=manifest.materials[name];
  assert.ok(definition,`Missing ${sex} material ${name}`);
  const record={color:definition.Kd??[255,255,255]};
  if(definition.map_Kd){
   const url=new URL(definition.map_Kd,'https://cdn.zygotebody.com');
   const basename=path.basename(url.pathname),ext=path.extname(basename),stem=basename.slice(0,-ext.length);
   const token=url.searchParams.get('token');
   const file=textureFiles.find(candidate=>candidate.startsWith(`${stem}__q_token_${token}_`)&&candidate.endsWith(ext));
   assert.ok(file,`Missing local texture for ${sex} material ${name}`);
   const output=`${sex}-texture-${textureIndex++}${ext.toLowerCase()}`;
   await fs.copyFile(path.join(textureDir,file),new URL(output,out));
   record.map=`/local-models/${output}`;
  }
  materials[name]=record;
 }
 const scale=1.75/(globalMax[1]-globalMin[1]);
 const midpointX=(globalMax[0]+globalMin[0])/2,midpointZ=(globalMax[2]+globalMin[2])/2;
 const parts=[],chunks=[],conceptMap=new Map();let segments=[],bytes=0,triangles=0;
 function append(values){
  const padding=(4-bytes%4)%4;if(padding){segments.push(Buffer.alloc(padding));bytes+=padding;}
  const offset=bytes;segments.push(Buffer.from(values.buffer,values.byteOffset,values.byteLength));bytes+=values.byteLength;return offset;
 }
 async function flush(){
  if(!bytes)return;
  const name=`${sex}-${chunks.length}.bin`,buffer=Buffer.concat(segments),gzip=gzipSync(buffer,{level:6});
  await fs.writeFile(new URL(name,out),buffer);await fs.writeFile(new URL(name+'.gz',out),gzip);
  chunks.push({url:`/local-models/${name}`,bytes,gzip:`/local-models/${name}.gz`,gzipBytes:gzip.length});
  segments=[];bytes=0;
 }
 for(const {mesh,attributes,indices} of geometries){
  let cursor=0;
  for(let segment=0;segment<mesh.names.length;segment++){
   const raw=mesh.names[segment],length=mesh.lengths[segment];
   assert.ok(leafNames.has(raw),`Unrecognized ${sex} entity: ${raw}`);
   assert.ok(length%3===0,`Non-triangular ${raw}`);
   const remap=new Map(),sourceIndices=indices.subarray(cursor,cursor+length);
   const compactIndices=new Uint32Array(length);
   for(let i=0;i<length;i++){
    const old=sourceIndices[i];
    if(!remap.has(old))remap.set(old,remap.size);
    compactIndices[i]=remap.get(old);
   }
   const positions=new Float32Array(remap.size*3),normals=new Int16Array(remap.size*3),uvs=new Float32Array(remap.size*2);
   const bounds=[[Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity]];
   for(const [old,next] of remap){
    uvs[next*2]=attributes[old*stride+3];uvs[next*2+1]=attributes[old*stride+4];
    for(let axis=0;axis<3;axis++){
     const value=axis===0?(attributes[old*stride]-midpointX)*scale:axis===1?(attributes[old*stride+1]-globalMin[1])*scale:(attributes[old*stride+2]-midpointZ)*scale;
     positions[next*3+axis]=value;
     normals[next*3+axis]=Math.round(Math.max(-1,Math.min(1,attributes[old*stride+5+axis]))*32767);
     bounds[0][axis]=Math.min(bounds[0][axis],value);bounds[1][axis]=Math.max(bounds[1][axis],value);
    }
   }
   if(bytes>4_000_000)await flush();
   const id=`LOCAL:${sex}:${raw}:${parts.length}`,conceptId=`LOCAL:${sex}:${raw}`;
   parts.push({id,conceptId,name:label(raw),system:systemFor(raw,category.get(raw)),material:mesh.material,chunk:chunks.length,
    positions:append(positions),normals:append(normals),uvs:append(uvs),indices:append(compactIndices),
    vertexCount:remap.size,indexCount:length,bounds,sourceId:raw});
   if(!conceptMap.has(conceptId))conceptMap.set(conceptId,{id:conceptId,name:label(raw),elements:[]});
   conceptMap.get(conceptId).elements.push(id);
   triangles+=length/3;cursor+=length;
  }
 }
 await flush();
 const concepts=[...conceptMap.values()];
 // Whole-organ handles make segmented source meshes easy to find and focus.
 for(const [name,pattern] of [
  ['Stomach',/^stomach\b/i],['Heart',/^heart\b/i],['Brain',/^brain\b/i],
  ['Lung',/^(left|right) lung\b/i],['Kidney',/^(left|right) kidney\b/i],
 ]){
  const elements=parts.filter(part=>pattern.test(part.name)).map(part=>part.id);
  if(elements.length&&!concepts.some(concept=>concept.name===name))concepts.push({id:`LOCAL:${sex}:organ:${name}`,name,elements});
 }
 const atlas={version:'local-zygote-2',sex,source:'Local Zygote Body mirror',scope:'Local use only',parts,concepts,chunks,materials,triangles};
 await fs.writeFile(new URL(`${sex}.json`,out),JSON.stringify(atlas));
 console.log(`${sex}: ${parts.length} parts, ${concepts.length} concepts, ${triangles.toLocaleString()} triangles, ${chunks.length} chunks`);
}

await fs.mkdir(out,{recursive:true});
for(const sex of ['male','female'])await importSex(sex);
