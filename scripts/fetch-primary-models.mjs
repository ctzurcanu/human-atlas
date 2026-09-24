/** Fetch the original publisher's Z-Anatomy Blender file and verify both hashes.
 * Usage: node scripts/fetch-primary-models.mjs [CACHE_DIRECTORY]
 * The archive is held outside public/ and no files from another atlas are used.
 */
import {createHash} from 'node:crypto';
import fs from 'node:fs';
import fsp from 'node:fs/promises';
import os from 'node:os';
import path from 'node:path';
import {pipeline} from 'node:stream/promises';
import {spawn} from 'node:child_process';

const source={
 url:'https://raw.githubusercontent.com/Z-Anatomy/Models-of-human-anatomy/refs/heads/master/Z-Anatomy.zip',
 archiveSha256:'e029688545627bd0214b269e1063143abb580aad72b2c2445d6d8a9a0d9da736',
 blendSha256:'9f08a17ea0115fed80b2a73ecdf0a1bc2ab2f6956f37c593ce23d513ea35afcd',
 member:'Z-Anatomy/Startup.blend',
 license:'CC BY-SA 4.0 with separately licensed components; see the original publisher notice',
};
const directory=path.resolve(process.argv[2]??path.join(os.tmpdir(),'human-atlas-primary-models'));
await fsp.mkdir(directory,{recursive:true});
const archive=path.join(directory,'Z-Anatomy.zip'),blend=path.join(directory,'Startup.blend');
async function sha256(filename){const hash=createHash('sha256');for await(const block of fs.createReadStream(filename))hash.update(block);return hash.digest('hex');}
async function valid(filename,expected){try{return await sha256(filename)===expected;}catch(error){if(error.code==='ENOENT')return false;throw error;}}
if(!await valid(archive,source.archiveSha256)){
 const temporary=archive+'.download';
 const curl=spawn('curl',['--fail','--location','--silent','--show-error','--output',temporary,source.url],{stdio:['ignore','ignore','inherit']});
 try{const code=await new Promise((resolve,reject)=>{curl.once('error',reject);curl.once('close',resolve)});
  if(code!==0)throw new Error(`Original archive download failed (curl exit ${code}).`);
  if(!await valid(temporary,source.archiveSha256))throw new Error('Original archive SHA-256 differs from the verified edition. Refusing to import.');
  await fsp.rename(temporary,archive);
 }finally{await fsp.rm(temporary,{force:true});}
}
if(!await valid(blend,source.blendSha256)){
 const temporary=blend+'.extract';
 const unzip=spawn('unzip',['-p',archive,source.member],{stdio:['ignore','pipe','inherit']});
 const completed=new Promise((resolve,reject)=>{unzip.once('error',reject);unzip.once('close',resolve)});
 try{await pipeline(unzip.stdout,fs.createWriteStream(temporary));
  const code=await completed;
  if(code!==0||!await valid(temporary,source.blendSha256))throw new Error('Original Blender file SHA-256 differs from the verified edition.');
  await fsp.rename(temporary,blend);
 }finally{await fsp.rm(temporary,{force:true});}
}
console.log(JSON.stringify({sourceUrl:source.url,archiveSha256:source.archiveSha256,blendSha256:source.blendSha256,license:source.license,blend},null,2));
