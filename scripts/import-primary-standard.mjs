/** Rebuild the standard male atlas from official BodyParts3D 4.0 files.
 * Usage: node scripts/import-primary-standard.mjs [--cache DIR] [--out DIR]
 * Optionally supply --archive ZIP and --tables DIR to use verified local files.
 */
import fs from 'node:fs/promises';
import fssync from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import {createHash} from 'node:crypto';
import {spawn} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const option=(flag,fallback)=>{const i=process.argv.indexOf(flag);return i<0?fallback:path.resolve(process.argv[i+1]);};
const cache=option('--cache',path.join(os.tmpdir(),'human-atlas-bodyparts3d-primary'));
const out=option('--out',path.join(root,'public','models'));
const archive=option('--archive',path.join(cache,'isa_BP3D_4.0_obj_99.zip'));
const tableDir=option('--tables',path.join(cache,'tables'));
const stage=path.join(cache,'output');
const base='https://dbarchive.biosciencedbc.jp/data/bodyparts3d/LATEST/';
const files={
 'isa_BP3D_4.0_obj_99.zip':'40665852c49f218326590e204db91064a1ecfc3c6f8cbd7bbbcaac62c7cd409e',
 'isa_parts_list_e.txt':'ab7796deedd49205e77f3609a1cb8c53e2bbee14ecb5c9a6ca05227469780513',
 'isa_element_parts.txt':'a3de74423f943b0d724ae8f59b3a817f87c423a544f8db98113b1980817cbeaf',
 'partof_parts_list_e.txt':'9224080557053e6f1322f1e13ab27f0ecde0db19bb3b505f0631afad230eeebd',
 'partof_element_parts.txt':'3f5f6df1028eb122b30de77c711597b6bb8e5541658e5985859fd228adbf88ea',
};
async function run(command,args){const child=spawn(command,args,{cwd:root,stdio:'inherit'});const code=await new Promise((resolve,reject)=>{child.once('error',reject);child.once('close',resolve)});if(code!==0)throw new Error(`${command} exited with ${code}`);}
async function sha256(filename){const hash=createHash('sha256');for await(const block of fssync.createReadStream(filename))hash.update(block);return hash.digest('hex');}
await fs.mkdir(tableDir,{recursive:true});await fs.mkdir(stage,{recursive:true});await fs.mkdir(out,{recursive:true});
for(const [name,expected] of Object.entries(files)){
 const filename=name.endsWith('.zip')?archive:path.join(tableDir,name);
 if(!fssync.existsSync(filename)){
  const temporary=filename+'.download';
  await run('curl',['--fail','--location','--silent','--show-error','--output',temporary,base+name]);
  if(await sha256(temporary)!==expected)throw new Error(`Publisher source changed: ${name}`);
  await fs.rename(temporary,filename);
 }
 if(await sha256(filename)!==expected)throw new Error(`Publisher source hash differs: ${name}`);
}
const objects=path.join(cache,'isa_BP3D_4.0_obj_99');
if(!fssync.existsSync(path.join(objects,'FJ1252.obj')))await run('unzip',['-q',archive,'-d',cache]);
await run(process.execPath,[path.join(root,'scripts','build-bodyparts3d-metadata.mjs'),tableDir,stage]);
await run('python3',[path.join(root,'scripts','convert-anatomy.py'),objects,path.join(stage,'concept-map.json'),path.join(stage,'system-map.json'),stage]);
await run(process.execPath,[path.join(root,'scripts','optimize-anatomy.mjs'),'atlas.json',stage]);
await run(process.execPath,[path.join(root,'scripts','compress-models.mjs'),stage]);
const manifestPath=path.join(stage,'atlas.json'),atlas=JSON.parse(await fs.readFile(manifestPath));
atlas.provenance={sourceUrl:'https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html',sourceHashes:files,systemMapping:'Human Atlas curated, scripts/data/bodyparts3d-elements.json'};
await fs.writeFile(manifestPath,JSON.stringify(atlas));
await run(process.execPath,[path.join(root,'scripts','validate-atlas.mjs'),'atlas.json',stage]);
for(const chunk of atlas.chunks)for(const url of [chunk.url,chunk.gzip]){const name=path.basename(url);await fs.copyFile(path.join(stage,name),path.join(out,name));}
await fs.copyFile(manifestPath,path.join(out,'atlas.json'));
if((await fs.stat(path.join(out,'atlas-z-anatomy.json')).catch(()=>null))?.isFile())
 await run(process.execPath,[path.join(root,'scripts','build-complete-male.mjs'),out]);
console.log(`Installed ${atlas.parts.length} independently rebuilt BodyParts3D meshes in ${out}`);
