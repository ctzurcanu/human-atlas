/** Rebuild the direct-source male atlas. No intermediary exports or catalogs.
 * Usage: node scripts/import-primary-male.mjs [--cache DIR] [--out DIR]
 * Set BLENDER to the Blender executable if it is not on PATH.
 */
import {spawn} from 'node:child_process';
import fs from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const option=(name,fallback)=>{const i=process.argv.indexOf(name);return i>=0?path.resolve(process.argv[i+1]):fallback;};
const cache=option('--cache',path.join(os.tmpdir(),'human-atlas-primary-models'));
const out=option('--out',path.join(root,'public','models'));
const stage=path.join(cache,'output-male');
const blender=process.env.BLENDER??(process.platform==='darwin'?'/Applications/Blender.app/Contents/MacOS/Blender':'blender');
async function run(command,args){const child=spawn(command,args,{cwd:root,stdio:'inherit'});const code=await new Promise((resolve,reject)=>{child.once('error',reject);child.once('close',resolve)});if(code!==0)throw new Error(`${command} exited with ${code}`);}
await fs.mkdir(out,{recursive:true});
await fs.mkdir(stage,{recursive:true});
await run(process.execPath,[path.join(root,'scripts','fetch-primary-models.mjs'),cache]);
await run(process.execPath,[path.join(root,'scripts','fetch-open3d-upper-limb.mjs'),cache]);
await run(blender,['-b',path.join(cache,'Startup.blend'),'--factory-startup','--disable-autoexec','--python',path.join(root,'scripts','build-z-anatomy.py'),'--','--out',stage,'--source-sha256','e029688545627bd0214b269e1063143abb580aad72b2c2445d6d8a9a0d9da736']);
await run(blender,['-b','--factory-startup','--disable-autoexec','--python',path.join(root,'scripts','build-open3d-upper-limb.py'),'--','--glb',path.join(cache,'upper-limb.glb'),'--atlas',path.join(stage,'atlas-z-anatomy.json'),'--zip-sha256','5af0190a6d7bf47393447ac30021e4f3ba619721c7f3a620c39a895947078432']);
const manifest=path.join(stage,'atlas-z-anatomy.json');
await run(process.execPath,[path.join(root,'scripts','prune-duplicate-upper-limb.mjs'),manifest]);
await run(process.execPath,[path.join(root,'scripts','validate-primary-atlas.mjs'),manifest]);
const atlas=JSON.parse(await fs.readFile(manifest));
for(const chunk of atlas.chunks)for(const url of [chunk.url,chunk.gzip]){const name=path.basename(url);await fs.copyFile(path.join(stage,name),path.join(out,name));}
for(const texture of atlas.provenance.additionalSources[0].preservedTextures){const name=path.basename(texture.url);await fs.copyFile(path.join(stage,name),path.join(out,name));}
await fs.copyFile(manifest,path.join(out,'atlas-z-anatomy.json'));
if((await fs.stat(path.join(out,'atlas.json')).catch(()=>null))?.isFile())
 await run(process.execPath,[path.join(root,'scripts','build-complete-male.mjs'),out]);
console.log(`Installed ${atlas.parts.length} publisher-source male surfaces in ${out}`);
