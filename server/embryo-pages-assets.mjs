import {readFile,readdir,mkdir,link,copyFile} from 'node:fs/promises';
import {resolve,sep} from 'node:path';
import {createHash} from 'node:crypto';

/** Explicit staging only; ordinary Pages builds never consume local review assets. */
export function embryoPagesAssets(mode){
 const directory=process.env.EMBRYO_PAGES_PACKAGE;
 return {name:'embryo-pages-package',apply:'build',async writeBundle(options){
  if(mode!=='embryo-local-preview'||!directory)return;
  async function linkedCopy(source,target){
   await mkdir(target,{recursive:true});
   for(const item of await readdir(source,{withFileTypes:true})){const from=resolve(source,item.name),to=resolve(target,item.name);
    if(item.isDirectory())await linkedCopy(from,to);else if(item.isFile()){try{await link(from,to);}catch(error){if(error.code==='EXDEV')await copyFile(from,to);else if(error.code!=='EEXIST')throw error;}}
   }
  }
  await linkedCopy(resolve(import.meta.dirname,'../public'),resolve(options.dir));
 },async generateBundle(_options,bundle){
  if(!directory)return;
  const root=resolve(directory),manifest=JSON.parse(await readFile(resolve(root,'package.json'),'utf8'));
  if(mode!=='embryo-local-preview'&&manifest.releasePermission?.status!=='verified'&&!(manifest.releasePermission?.status==='source-preserving-noncommercial-candidate'&&manifest.sourcePreservation?.allOriginalFacesRetained))throw Error('Embryo publication blocked: adapted CS23 source requires redistribution permission. Use embryo-local-preview for local review only.');
  const files={};
  for(const item of manifest.files){
   const file=resolve(root,item.path);
   if(!file.startsWith(root+sep)||!/^models\/embryo-3month\/[a-zA-Z0-9._-]+$/.test(item.path))throw Error('Invalid embryo package path');
   const bytes=await readFile(file),hash=createHash('sha256').update(bytes).digest('hex');
   if(hash!==item.sha256)throw Error('Stale embryo package: '+item.path);
   this.emitFile({type:'asset',fileName:item.path,source:bytes});files[item.path]={version:hash,bytes:bytes.length};
  }
  // Include staged assets in the same content-version cache as public assets.
  const asset=bundle['asset-manifest.json'];if(!asset||asset.type!=='asset')throw Error('Missing base asset manifest');
  const assets=JSON.parse(String(asset.source));Object.assign(assets.files,files);asset.source=JSON.stringify(assets);
  this.emitFile({type:'asset',fileName:'embryo-release-status.json',source:JSON.stringify({model:manifest.model,catalogueSha256:manifest.catalogueSha256,releasePermission:manifest.releasePermission,localPreviewOnly:mode==='embryo-local-preview'},null,2)});
 }};
}
