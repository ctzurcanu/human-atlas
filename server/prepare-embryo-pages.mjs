import {originalEmbryoFiles} from './embryo-publication-source.mjs';
import {readFile,writeFile,mkdir,unlink} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {resolve,basename} from 'node:path';
import {gzipSync,gunzipSync} from 'node:zlib';
const root=resolve(import.meta.dirname,'..'),source=resolve(root,'.local-models/embryo-cs23'),out=resolve(root,'outputs/embryo-pages-package');
const original=await readFile(resolve(source,'atlas.json')),atlas=JSON.parse(original),files=[];
const prefix='models/embryo-3month/';await mkdir(resolve(out,prefix),{recursive:true});
async function put(path,bytes){await writeFile(resolve(out,path),bytes);files.push({path,bytes:bytes.length,sha256:createHash('sha256').update(bytes).digest('hex')});}
for(const chunk of atlas.chunks.slice(0,8)){
 if(!chunk.url.startsWith('/local-models/embryo-cs23/'))throw Error('Unexpected nonlocal chunk');
 const name=basename(chunk.url);if(!/^mesh-\d+-[a-f0-9]+\.bin$/.test(name))throw Error('Invalid chunk');
 const raw=await readFile(resolve(source,name)),compressed=await readFile(resolve(source,name+'.gz'));
 if(!gunzipSync(compressed).equals(raw)||raw.length!==chunk.bytes)throw Error('Chunk gzip/size mismatch');
 await put(prefix+name+'.gz',compressed);
 chunk.url='/'+prefix+name+'.gz';chunk.gzip=chunk.url;
}
const sourcePreservation=await originalEmbryoFiles(root,atlas,put,prefix);
for(const material of Object.values(atlas.materials))for(const key of ['map','normalMap'])if(material[key])throw Error('External material asset needs explicit packaging: '+material[key]);
atlas.displayName='Embryo 3month';atlas.version='Embryo 3month · CS23 specimen 9226 · partial joined assembly';
atlas.publication={status:'source-preserving-noncommercial-candidate',sourceLicense:'CC BY-NC-ND 4.0',sourceLicenseUrl:'https://www.3dembryoatlas.com/blank',stage:'CS23, 56–60 days after fertilization; display name is not source stage evidence'};
const catalogue=Buffer.from(JSON.stringify(atlas)+'\n');await put(prefix+'atlas.json',catalogue);await put(prefix+'atlas.json.gz',gzipSync(catalogue,{mtime:0}));
const attribution=`# Embryo 3month — local publication candidate\n\n3D Atlas of Human Embryology, Carnegie stage 23, specimen 9226, Department of Medical Biology, Academic Medical Center, University of Amsterdam (2016). Source: https://www.3dembryoatlas.com/blank — CC BY-NC-ND 4.0.\n\nHuman Reference Atlas placental context: Kristen Browne and Heidi Schlehlein / HuBMAP, CC BY 4.0; see the site's ATTRIBUTION.md.\n\nCS23 source: every original face and source material color is retained, with technical browser packing and a common rigid placement. HRA context is a separate neighboring model, with its cord and placental geometry adjusted independently under CC BY 4.0. Source stage is 56–60 days after fertilization; the requested display name does not establish a three-month specimen. Mixed-stage surfaces, open rims, inherited defects and wall/lumen representations remain provisional.\n\nThe CS23 source retains CC BY-NC-ND 4.0. This is a source-preserving noncommercial publication candidate, with separate source attribution and license terms; no permission for commercial use or altered CS23 geometry is claimed.\n`;
await put(prefix+'ATTRIBUTION.md',Buffer.from(attribution));
const manifest={model:'Embryo 3month',catalogueSha256:createHash('sha256').update(catalogue).digest('hex'),sourceCatalogueSha256:createHash('sha256').update(original).digest('hex'),parts:atlas.parts.length,triangles:atlas.triangles,sourcePreservation,releasePermission:{status:'source-preserving-noncommercial-candidate',source:'https://www.3dembryoatlas.com/blank',reason:'All original CS23 faces/material colors preserved. CC BY-NC-ND terms, including noncommercial use, remain in force; separately licensed HRA context stays independent.'},files};
if(process.argv.includes('--stage-public')){
 // Remove only previously generated raw duplicates whose recorded hash still matches.
 let old;try{old=JSON.parse(await readFile(resolve(root,'public',prefix,'release-manifest.json'),'utf8'));}catch(error){if(error.code!=='ENOENT')throw error;}
 for(const item of old?.files??[])if(item.path.startsWith(prefix)&&item.path.endsWith('.bin')){const file=resolve(root,'public',item.path),bytes=await readFile(file);if(createHash('sha256').update(bytes).digest('hex')!==item.sha256)throw Error('Generated buffer was changed; refusing cleanup');await unlink(file);}

 for(const item of files){const bytes=await readFile(resolve(out,item.path));if(createHash('sha256').update(bytes).digest('hex')!==item.sha256)throw Error('Staging hash mismatch');await mkdir(resolve(root,'public',prefix),{recursive:true});await writeFile(resolve(root,'public',item.path),bytes);}
 await writeFile(resolve(root,'public',prefix,'release-manifest.json'),JSON.stringify(manifest,null,2)+'\n');
}
await writeFile(resolve(out,'package.json'),JSON.stringify(manifest,null,2)+'\n');console.log(JSON.stringify({directory:out,parts:atlas.parts.length,triangles:atlas.triangles,files:files.length,releasePermission:manifest.releasePermission},null,2));
