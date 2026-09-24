import {fileURLToPath} from 'node:url';
import {createReadStream} from 'node:fs';
import {stat} from 'node:fs/promises';
import {resolve,sep} from 'node:path';
import {defineConfig,type Plugin} from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/postcss';
const path=(relative:string)=>fileURLToPath(new URL(relative,import.meta.url));
function localModels():Plugin{
 const directory=path('./.local-models');
 return {name:'local-models-dev-only',apply:'serve',configureServer(server){
  server.middlewares.use('/local-models',async(req,res,next)=>{
   let name:string;
   try{name=decodeURIComponent((req.url??'').split('?')[0]).replace(/^\/+/, '');}catch{res.statusCode=400;res.end();return;}
   const file=resolve(directory,name);
   if(!name||!file.startsWith(directory+sep)||!/\.(json|bin|bin\.gz|jpe?g|png|webp)$/.test(name)){next();return;}
   try{
    const info=await stat(file);if(!info.isFile()){next();return;}
    res.setHeader('Content-Type',name.endsWith('.json')?'application/json':/\.jpe?g$/.test(name)?'image/jpeg':name.endsWith('.png')?'image/png':name.endsWith('.webp')?'image/webp':'application/octet-stream');
    res.setHeader('Content-Length',info.size);
    createReadStream(file).pipe(res);
   }catch{next();}
  });
 }};
}
export default defineConfig({base:process.env.VITE_BASE_PATH||'/',root:path('./web'),publicDir:path('./public'),plugins:[react(),localModels()],resolve:{alias:{'@':path('./')}},css:{postcss:{plugins:[tailwindcss()]}},server:{watch:{usePolling:true}},build:{outDir:path('./dist'),emptyOutDir:true}});
