import {createHash} from 'node:crypto';
import {readdir,readFile,rename,unlink,writeFile} from 'node:fs/promises';
import {resolve,sep} from 'node:path';
import {gzipSync,gunzipSync} from 'node:zlib';

/** Compact only the generated Pages site. Model loaders support standalone
 * gzip banks, and some catalogues share the same bank. Verify before pruning. */
export async function compactPagesGeometry(directory) {
 const version=bytes=>({version:createHash('sha256').update(bytes).digest('hex'),bytes:bytes.length});
 async function replace(file,bytes) {
  const temporary=file+'.'+process.pid+'.tmp';
  await writeFile(temporary,bytes);
  await rename(temporary,file);
 }
 const manifestFile=resolve(directory,'asset-manifest.json');
 const manifest=JSON.parse(await readFile(manifestFile,'utf8'));
 const rawFiles=new Set();
 const models=resolve(directory,'models');
 for (const entry of await readdir(models,{recursive:true,withFileTypes:true})) {
  if (!entry.isFile() || !entry.name.endsWith('.json')) continue;
  const file=resolve(entry.parentPath,entry.name), atlas=JSON.parse(await readFile(file,'utf8'));
  if (!Array.isArray(atlas.parts) || !Array.isArray(atlas.chunks)) continue;
  let changed=false;
  for (const chunk of atlas.chunks??[]) {
   if (!chunk.url?.startsWith('/models/') || !chunk.url.endsWith('.bin')) continue;
   const rawFile=resolve(directory,'.'+chunk.url);
   if (!rawFile.startsWith(models+sep)) throw Error('Invalid Pages bank path');
   const zipped=chunk.url+'.gz', gzipFile=rawFile+'.gz';
   if (!rawFiles.has(rawFile)) {
    const raw=await readFile(rawFile);
    if (raw.length!==chunk.bytes) throw Error('Incomplete Pages bank '+chunk.url);
    const compressed=chunk.gzip===zipped ? await readFile(gzipFile) : gzipSync(raw,{level:6});
    if (!gunzipSync(compressed).equals(raw)) throw Error('Pages gzip mismatch '+zipped);
    if (chunk.gzip!==zipped) await replace(gzipFile,compressed);
    manifest.files[zipped.slice(1)]=version(compressed);
    rawFiles.add(rawFile);
   }
   chunk.url=zipped; chunk.gzip=zipped; chunk.gzipBytes=manifest.files[zipped.slice(1)].bytes;
   changed=true;
  }
  if (changed) {
   const bytes=Buffer.from(JSON.stringify(atlas));
   await replace(file,bytes);
   manifest.files[file.slice(directory.length+1).split(sep).join('/')]=version(bytes);
  }
 }
 for (const file of rawFiles) {
  await unlink(file);
  delete manifest.files[file.slice(directory.length+1).split(sep).join('/')];
 }
 await replace(manifestFile,JSON.stringify(manifest));
}
