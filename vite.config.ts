import {fileURLToPath} from 'node:url';
import {createReadStream} from 'node:fs';
import {stat} from 'node:fs/promises';
import {resolve,sep} from 'node:path';
import {defineConfig,type Plugin} from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/postcss';
import {attachConnectRelay,requestIp} from './server/connect-relay.mjs';
import {connectionIdentity} from './server/connect-address.mjs';
import {connectIdentityCors} from './server/connect-origin.mjs';
import {assetManifest} from './server/asset-manifest.mjs';
const path=(relative:string)=>fileURLToPath(new URL(relative,import.meta.url));
function persistentAssets():Plugin{
 return {name:'atlas-asset-versions',configureServer(server){
  server.middlewares.use('/asset-manifest.json',async(_req,res)=>{
   try{res.setHeader('Content-Type','application/json');res.setHeader('Cache-Control','no-cache');res.end(JSON.stringify(await assetManifest(path('./public'),path('./.local-models'))));}
   catch{res.statusCode=503;res.end();}
  });
 },async generateBundle(){this.emitFile({type:'asset',fileName:'asset-manifest.json',source:JSON.stringify(await assetManifest(path('./public')))});}};
}
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
function connectRelay():Plugin{
 const identity=async(req:import('node:http').IncomingMessage,res:import('node:http').ServerResponse)=>{if(!connectIdentityCors(req,res)){res.statusCode=403;res.end();return;}res.setHeader('Content-Type','application/json');res.setHeader('Cache-Control','no-store');res.end(JSON.stringify({ip:requestIp(req),...await connectionIdentity(req)}));};
 return {name:'atlas-connect-relay',configureServer(server){server.middlewares.use('/atlas-connect/identity',identity);if(server.httpServer)attachConnectRelay(server.httpServer);},configurePreviewServer(server){server.middlewares.use('/atlas-connect/identity',identity);attachConnectRelay(server.httpServer);}};
}
export default defineConfig({base:process.env.VITE_BASE_PATH||'/',root:path('./web'),publicDir:path('./public'),plugins:[react(),localModels(),connectRelay(),persistentAssets()],resolve:{alias:{'@':path('./')}},css:{postcss:{plugins:[tailwindcss()]}},server:{allowedHosts:true,watch:{usePolling:true}},build:{outDir:path('./dist'),emptyOutDir:true}});
