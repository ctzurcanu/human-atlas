import {maleTa98Directory,stagedMaleTa98} from './server/male-ta98-pages';
import {menuItemCounts} from './server/model-menu';
import {embryoPagesAssets} from './server/embryo-pages-assets.mjs';
import {createHash} from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {createReadStream} from 'node:fs';
import {stat,writeFile,readFile} from 'node:fs/promises';
import {resolve,sep} from 'node:path';
import {promisify} from 'node:util';
import {gzip} from 'node:zlib';
import {defineConfig,type Plugin,type UserConfig} from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/postcss';
import {attachConnectRelay,requestIp} from './server/connect-relay.mjs';
import {connectionIdentity} from './server/connect-address.mjs';
import {connectIdentityCors} from './server/connect-origin.mjs';
import {assetManifest} from './server/asset-manifest.mjs';
import {prepareComparisonPage} from './server/comparison-page.mjs';
import {catalogueAttribution,localAttributionPath,projectCatalogueAttribution} from './server/catalogue-attribution.mjs';
const path=(relative:string)=>fileURLToPath(new URL(relative,import.meta.url));
const compress=promisify(gzip);
function persistentAssets(maleDirectory?:string):Plugin{
 return {name:'atlas-asset-versions',configureServer(server){
  server.middlewares.use('/asset-manifest.json',async(_req,res)=>{
   try{res.setHeader('Content-Type','application/json');res.setHeader('Cache-Control','no-cache');res.end(JSON.stringify(await assetManifest(path('./public'),path('./.local-models'),maleDirectory)));}
   catch{res.statusCode=503;res.end();}
  });
 },async generateBundle(){this.emitFile({type:'asset',fileName:'asset-manifest.json',source:JSON.stringify(await assetManifest(path('./public')))});this.emitFile({type:'asset',fileName:'.nojekyll',source:''});}};
}
function localModels(sources:Awaited<ReturnType<typeof catalogueAttribution>>):Plugin{
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
   let file=resolve(directory,name);
   if(!name||!file.startsWith(directory+sep)||!/\.(json|gltf|bin|bin\.gz|jpe?g|png|webp)$/.test(name)){next();return;}
   try{
    try{await stat(file);}catch{
     file=resolve(directory,localAttributionPath(name,sources));
     if(!file.startsWith(directory+sep)){next();return;}
    }
    let served=file;
    let catalogue:unknown;
    const acceptsGzip=String(req.headers['accept-encoding']??'').split(',').some(value=>value.trim().split(';')[0]==='gzip'&&!/;\s*q=0(?:\.0*)?\s*$/i.test(value));
    if(name==='ta98-runtime.json'||name==='ta98-female-runtime.json'){
     const runtime=JSON.parse(await readFile(file,'utf8'));
     const hash=createHash('sha256').update(await readFile(resolve(directory,name==='ta98-female-runtime.json'?'ta98-female.json':'ta98-male.json'))).digest('hex');
     if(runtime.delivery?.sourceAtlasSha256!==hash){res.statusCode=409;res.end('Runtime is stale. Rebuild it from the current TA98 master.');return;}
     catalogue=runtime;
    }
    if(name.endsWith('.json')){
     catalogue??=JSON.parse(await readFile(file,'utf8'));
     const projected=projectCatalogueAttribution(catalogue,sources);
     if(projected!==catalogue){
      res.setHeader('Content-Type','application/json');
      res.setHeader('Cache-Control','no-cache');
      res.setHeader('Vary','Accept-Encoding');
      const body=JSON.stringify(projected);
      if(acceptsGzip){
       const compressed=await compress(body);
       res.setHeader('Content-Encoding','gzip');
       res.end(compressed);
      }else res.end(body);
      return;
     }
    }
    // Keep large local geometry transfers small; Fetch decodes Content-Encoding
    // before the model loader verifies the original uncompressed byte count.
    if(name.endsWith('.bin')&&acceptsGzip)try{await stat(file+'.gz');served=file+'.gz';res.setHeader('Content-Encoding','gzip');res.setHeader('Vary','Accept-Encoding');}catch{}
    const info=await stat(served);if(!info.isFile()){next();return;}
    res.setHeader('Content-Type',name.endsWith('.json')?'application/json':/\.jpe?g$/.test(name)?'image/jpeg':name.endsWith('.png')?'image/png':name.endsWith('.webp')?'image/webp':'application/octet-stream');
    if(name.startsWith('dermatomes/'))res.setHeader('Cache-Control','no-store');
    res.setHeader('Content-Length',info.size);
    const stream=createReadStream(served);
    // A demand change or iframe navigation can cancel a large download.
    // pipe() alone leaves its paused source open after the response closes.
    const cancel=()=>stream.destroy();
    req.once('aborted',cancel);res.once('close',cancel);
    stream.once('close',()=>{req.off('aborted',cancel);res.off('close',cancel);});
    stream.once('error',error=>{
     if(res.headersSent){res.destroy(error);return;}
     res.statusCode=500;res.removeHeader('Content-Length');res.removeHeader('Content-Encoding');res.end('Could not read the local anatomy asset.');
    });
    stream.pipe(res);
   }catch{next();}
  });
 }};
}
function connectRelay():Plugin{
 const identity=async(req:import('node:http').IncomingMessage,res:import('node:http').ServerResponse)=>{if(!connectIdentityCors(req,res)){res.statusCode=403;res.end();return;}res.setHeader('Content-Type','application/json');res.setHeader('Cache-Control','no-store');res.end(JSON.stringify({ip:requestIp(req),...await connectionIdentity(req)}));};
 return {name:'atlas-connect-relay',configureServer(server){server.middlewares.use('/atlas-connect/identity',identity);if(server.httpServer)attachConnectRelay(server.httpServer);},configurePreviewServer(server){server.middlewares.use('/atlas-connect/identity',identity);attachConnectRelay(server.httpServer);}};
}
function recordingUpload():Plugin{return {name:'atlas-recording-upload',apply:'serve',configureServer(server){server.middlewares.use('/recording-upload',async(req,res,next)=>{if(req.method!=='POST'){next();return;}const chunks:Buffer[]=[];for await(const chunk of req)chunks.push(Buffer.isBuffer(chunk)?chunk:Buffer.from(chunk));await writeFile(resolve('outputs/atlas-demo-dark.webm'),Buffer.concat(chunks));res.statusCode=201;res.end('saved');});}};}
export default defineConfig(async({mode,command}):Promise<UserConfig>=>{
 const root=resolve(path('./')),dev=command==='serve';
 const maleDirectory=maleTa98Directory(root,dev,mode);
 const counts=await menuItemCounts(root,dev,maleDirectory,process.env.VITE_EMBRYO_CATALOGUE);
 await prepareComparisonPage(path('./'));const sources=await catalogueAttribution(path('./ATTRIBUTION.md'));return {define:{'import.meta.env.VITE_MODEL_ITEM_COUNTS':JSON.stringify(counts)},base:process.env.VITE_BASE_PATH||'/',root:path('./web'),publicDir:path('./public'),plugins:[react(),localModels(sources),connectRelay(),persistentAssets(maleDirectory),stagedMaleTa98(root,maleDirectory),embryoPagesAssets(mode) as Plugin,recordingUpload()],resolve:{alias:{'@':path('./')}},css:{postcss:{plugins:[tailwindcss()]}},server:{allowedHosts:true,watch:{usePolling:true}},build:{copyPublicDir:mode!=='embryo-local-preview',outDir:path('./dist'),emptyOutDir:true,rollupOptions:{input:{atlas:path('./web/index.html'),compare:path('./web/compare.html')}}}};});
