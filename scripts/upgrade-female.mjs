/** Import only female-source geometry from a pinned, CC BY 4.0 adaptation.
 * Run: node scripts/upgrade-female.mjs
 * Original sources and adaptations are recorded in public/ATTRIBUTION.md.
 * Repack the allowed meshes: excluded male geometry never enters our output.
 */
import fs from 'node:fs/promises';
import assert from 'node:assert/strict';
import {gzipSync,gunzipSync} from 'node:zlib';

const revision='5bb5713aab18d7fe9380c3339eb09f173491ea06';
const base=`https://raw.githubusercontent.com/slorksmo/Human-Atlas/${revision}/public/models/`;
const out=new URL('../public/models/',import.meta.url);
async function download(name){
 const response=await fetch(base+name);
 if(!response.ok)throw new Error(`${name}: HTTP ${response.status}`);
 return Buffer.from(await response.arrayBuffer());
}
const input=JSON.parse((await download('atlas-female.json')).toString());
assert.equal(input.sex,'female');
const femaleParts=input.parts.filter(p=>p.system!=='borrowed');
assert.equal(femaleParts.length,1040);
assert.ok(femaleParts.every(p=>/^HRAF\d+$|^VHF\d+$/.test(p.id)));
assert.equal(femaleParts.filter(p=>p.system==='donor-muscle').length,76);

// Use the study's full rectus femoris on each side, avoiding duplicate surfaces.
// Keep the HRA concepts searchable by resolving them to the corresponding muscle.
const replacements=new Map([['HRAF0394','VHF0009'],['HRAF0396','VHF0047']]);
for(const [from,to] of replacements){
 assert.equal(femaleParts.find(p=>p.id===from)?.name,femaleParts.find(p=>p.id===to)?.name);
}
const retained=femaleParts.filter(p=>!replacements.has(p.id));
const ids=new Set(retained.map(p=>p.id));
const concepts=input.concepts.map(c=>({...c,elements:[...new Set(c.elements.map(id=>replacements.get(id)??id).filter(id=>ids.has(id)))]})).filter(c=>c.elements.length);
const buffers=new Map();
for(const index of new Set(retained.map(p=>p.chunk))){
 const chunk=input.chunks[index];
 const buffer=gunzipSync(await download(chunk.gzip.split('/').pop()));
 assert.equal(buffer.length,chunk.bytes);
 buffers.set(index,buffer);
 console.log(`Loaded female geometry chunk ${index}`);
}
const chunks=[],parts=[];let segments=[],bytes=0;
async function flush(){
 if(!bytes)return;
 const buffer=Buffer.concat(segments),compressed=gzipSync(buffer,{level:9});
 const name=`female-expanded-${chunks.length}.bin`;
 await fs.writeFile(new URL(name,out),buffer);
 await fs.writeFile(new URL(name+'.gz',out),compressed);
 chunks.push({url:'/models/'+name,bytes,gzip:'/models/'+name+'.gz',gzipBytes:compressed.length});
 segments=[];bytes=0;
}
function append(buffer,start,length){
 assert.ok(start>=0&&start+length<=buffer.length);
 const padding=(4-bytes%4)%4;
 if(padding){segments.push(Buffer.alloc(padding));bytes+=padding;}
 const offset=bytes;segments.push(buffer.subarray(start,start+length));bytes+=length;return offset;
}
for(const part of retained){
 if(bytes>4_000_000)await flush();
 const buffer=buffers.get(part.chunk),study=part.system==='donor-muscle';
 const carried=part.id.startsWith('HRAF')&&Number(part.id.slice(4))>=956;
 const provenance=study?{
  label:'Visible Human Female · lower-limb study',url:'https://digitalcommons.du.edu/visiblehuman/1/',
  detail:'Andreassen et al. (2023). Female-source muscle fitted to the HRA reference bones; approximate placement, with source fit residuals of 17–36 mm.'
 }:{label:`Human Reference Atlas ${carried?'v1.5':'v1.10'}`,url:`https://purl.humanatlas.io/ref-organ/united-female/${carried?'v1.5':'v1.10'}`,
  detail:carried?'Female pelvic geometry retained from the earlier HRA release.':'Female reference organ assembly.'};
 parts.push({...part,system:study?'muscular':part.system==='brain'?'nervous':part.system,provenance,chunk:chunks.length,
  positions:append(buffer,part.positions,part.vertexCount*12),
  normals:append(buffer,part.normals,part.vertexCount*6),
  indices:append(buffer,part.indices,part.indexCount*4)});
}
await flush();
const atlas={version:'HRA female v1.10 + Visible Human Female lower-limb muscles',sex:'female',
 source:'Human Reference Atlas / Visible Human Female',
 scope:'Female-source reference with expanded lower-limb muscles. Upper-body skeleton and muscle coverage remains incomplete.',
 parts,concepts,chunks,triangles:parts.reduce((sum,p)=>sum+p.indexCount/3,0),
 provenance:{geometryRevision:revision,femaleOnly:true,excludedMaleMeshes:180,replacedDuplicateMuscles:2,lowerLimbFit:input.donorMuscle}};
await fs.writeFile(new URL('atlas-female.json',out),JSON.stringify(atlas));
console.log(JSON.stringify({parts:parts.length,concepts:concepts.length,triangles:atlas.triangles,systems:Object.fromEntries([...new Set(parts.map(p=>p.system))].map(id=>[id,parts.filter(p=>p.system===id).length]))},null,2));
