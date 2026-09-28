import {createHash} from 'node:crypto';
import {createReadStream} from 'node:fs';
import {readdir,stat} from 'node:fs/promises';
import {join} from 'node:path';

const supported=/\.(json|bin|gz|png|jpe?g|webp)$/i;
async function collect(directory,prefix,files){
 let entries;try{entries=await readdir(directory,{withFileTypes:true});}catch(error){if(error.code==='ENOENT')return;throw error;}
 for(const entry of entries){
  const file=join(directory,entry.name),key=prefix+'/'+entry.name;
  if(entry.isDirectory())await collect(file,key,files);
  else if(entry.isFile()&&supported.test(entry.name))files.push({file,key});
 }
}

/** Deployment hashes stay stable across viewer-only builds; dev versions follow file edits. */
export async function assetManifest(publicDirectory,localDirectory){
 const files=[];
 await collect(join(publicDirectory,'models'),'models',files);
 await collect(join(publicDirectory,'assets'),'assets',files);
 if(localDirectory)await collect(localDirectory,'local-models',files);
 const versions={};let cursor=0;
 await Promise.all(Array.from({length:6},async()=>{
  while(cursor<files.length){
   const {file,key}=files[cursor++],info=await stat(file);let version;
   if(localDirectory)version=`dev:${info.size}:${info.mtimeMs}`;
   else{const hash=createHash('sha256');for await(const chunk of createReadStream(file))hash.update(chunk);version=hash.digest('hex');}
   versions[key]={version,bytes:info.size};
  }
 }));
 return {schema:1,files:Object.fromEntries(Object.entries(versions).sort(([a],[b])=>a.localeCompare(b)))};
}
