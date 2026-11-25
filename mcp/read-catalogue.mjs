import {readFileSync,statSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {build} from 'esbuild';
import {MODELS} from './atlas-core.mjs';

// Bundle the viewer's TypeScript catalogue rules for Node. The Worker receives the resulting data only.
const bundle=await build({entryPoints:[fileURLToPath(new URL('./viewer-catalogue.ts',import.meta.url))],bundle:true,platform:'node',format:'esm',write:false});
const {selectableCatalogue}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
function catalogueUrl(model,filename){
 const directory=model.startsWith('local-')?'../.local-models/':'../public/models/';
 return new URL(directory+filename,import.meta.url);
}
export function catalogueRevision(model,filename){
 const {ino,size,mtimeNs}=statSync(catalogueUrl(model,filename),{bigint:true});
 return [ino,size,mtimeNs].join(':');
}
export function readCatalogue(model,filename){
 const atlas=JSON.parse(readFileSync(catalogueUrl(model,filename),'utf8'));
 return model==='male-ta98'?selectableCatalogue(atlas):atlas;
}

export function publicModelCounts(){
 const counts={},byFile=new Map();
 for(const [model,file] of Object.entries(MODELS)){
  if(model.startsWith('local-'))continue;
  if(!byFile.has(file)){
   const atlas=readCatalogue(model,file);
   byFile.set(file,{total:atlas.totalParts??atlas.parts.length,selectable:atlas.parts.filter(p=>!p.suppressed).length});
  }
  counts[model]=byFile.get(file);
 }
 return counts;
}
