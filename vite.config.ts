import {embryoPagesAssets} from './server/embryo-pages-assets.mjs';
import {createHash} from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {createReadStream} from 'node:fs';
import {stat,writeFile,readFile} from 'node:fs/promises';
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
 },async generateBundle(){this.emitFile({type:'asset',fileName:'asset-manifest.json',source:JSON.stringify(await assetManifest(path('./public')))});this.emitFile({type:'asset',fileName:'.nojekyll',source:''});}};
}
function localModels():Plugin{
 const directory=path('./.local-models');
 return {name:'local-models-dev-only',apply:'serve',configureServer(server){
  server.middlewares.use('/wall-texture-evidence',async(req,res)=>{
   const origin=String(req.headers.origin??'');
   if(req.method!=='POST'||!/^http:\/\/(localhost|127\.0\.0\.1):\d+$/.test(origin)){res.statusCode=403;res.end();return;}
   try{
    const chunks:Buffer[]=[];let size=0;
    for await(const chunk of req){const bytes=Buffer.isBuffer(chunk)?chunk:Buffer.from(chunk);size+=bytes.length;if(size>512000)throw Error('Too large');chunks.push(bytes);}
    const data=JSON.parse(Buffer.concat(chunks).toString('utf8'));
    const planBytes=await readFile(resolve(directory,'ta98-review/wall-texture-gpu-plan.json'));
    if(data.status!=='measured-local-wall-three-gpu'||data.planSha256!==createHash('sha256').update(planBytes).digest('hex')||data.records?.length!==JSON.parse(planBytes.toString('utf8')).expectedCases)throw Error('Invalid evidence');
    await writeFile(path('./reports/ta98-female/wall-texture-gpu.json'),JSON.stringify(data,null,2)+'\n');res.statusCode=201;res.end('saved');
   }catch{res.statusCode=400;res.end('Invalid review evidence');}
  });
  server.middlewares.use('/local-models',async(req,res,next)=>{
   let name:string;
   try{name=decodeURIComponent((req.url??'').split('?')[0]).replace(/^\/+/, '');}catch{res.statusCode=400;res.end();return;}
   const file=resolve(directory,name);
   if(!name||!file.startsWith(directory+sep)||!/\.(json|gltf|bin|bin\.gz|jpe?g|png|webp)$/.test(name)){next();return;}
   try{
    let served=file;
    if(name==='ta98-runtime.json'||name==='ta98-female-runtime.json'){
     const runtime=JSON.parse(await readFile(file,'utf8'));
     const hash=createHash('sha256').update(await readFile(resolve(directory,name==='ta98-female-runtime.json'?'ta98-female.json':'ta98-male.json'))).digest('hex');
     if(runtime.delivery?.sourceAtlasSha256!==hash){res.statusCode=409;res.end('Runtime is stale. Rebuild it from the current TA98 master.');return;}
    }
    if(name==='ta98-runtime.json'||name==='ta98-female-runtime.json')try{await stat(file+'.gz');served=file+'.gz';res.setHeader('Content-Encoding','gzip');}catch{}
    const info=await stat(served);if(!info.isFile()){next();return;}
    res.setHeader('Content-Type',name.endsWith('.json')?'application/json':/\.jpe?g$/.test(name)?'image/jpeg':name.endsWith('.png')?'image/png':name.endsWith('.webp')?'image/webp':'application/octet-stream');
    if(name.startsWith('dermatomes/'))res.setHeader('Cache-Control','no-store');
    res.setHeader('Content-Length',info.size);
    createReadStream(served).pipe(res);
   }catch{next();}
  });
 }};
}
function connectRelay():Plugin{
 const identity=async(req:import('node:http').IncomingMessage,res:import('node:http').ServerResponse)=>{if(!connectIdentityCors(req,res)){res.statusCode=403;res.end();return;}res.setHeader('Content-Type','application/json');res.setHeader('Cache-Control','no-store');res.end(JSON.stringify({ip:requestIp(req),...await connectionIdentity(req)}));};
 return {name:'atlas-connect-relay',configureServer(server){server.middlewares.use('/atlas-connect/identity',identity);if(server.httpServer)attachConnectRelay(server.httpServer);},configurePreviewServer(server){server.middlewares.use('/atlas-connect/identity',identity);attachConnectRelay(server.httpServer);}};
}
function recordingUpload():Plugin{return {name:'atlas-recording-upload',apply:'serve',configureServer(server){server.middlewares.use('/recording-upload',async(req,res,next)=>{if(req.method!=='POST'){next();return;}const chunks:Buffer[]=[];for await(const chunk of req)chunks.push(Buffer.isBuffer(chunk)?chunk:Buffer.from(chunk));await writeFile(resolve('outputs/atlas-demo-dark.webm'),Buffer.concat(chunks));res.statusCode=201;res.end('saved');});}};}
export default defineConfig(({mode})=>({base:process.env.VITE_BASE_PATH||'/',root:path('./web'),publicDir:path('./public'),plugins:[react(),localModels(),connectRelay(),persistentAssets(),embryoPagesAssets(mode),recordingUpload()],resolve:{alias:{'@':path('./')}},css:{postcss:{plugins:[tailwindcss()]}},server:{allowedHosts:true,watch:{usePolling:true}},build:{copyPublicDir:mode!=='embryo-local-preview',outDir:path('./dist'),emptyOutDir:true}}));
