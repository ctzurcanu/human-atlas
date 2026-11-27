import {createReadStream,existsSync} from 'node:fs';
import {readdir,readFile,stat} from 'node:fs/promises';
import {join,resolve,sep} from 'node:path';
import type {Plugin} from 'vite';
import {createHash} from 'node:crypto';
import {compactPagesGeometry} from './pages-geometry.mjs';

const prefix='models/male-ta98';
export function maleTa98Directory(root:string, dev:boolean, mode:string) {
 const installed=join(root,'public',prefix);
 const staged=join(root,'outputs/male-ta98-pages',prefix);
 if (existsSync(join(installed,'atlas.json'))) return installed;
 if ((dev || mode==='male-ta98-preview') && existsSync(join(staged,'atlas.json'))) return staged;
 if (mode==='male-ta98-preview') throw Error('Prepare Male TA98 with npm run prepare:male-ta98 first.');
}

export function stagedMaleTa98(root:string, directory:string|undefined):Plugin {
 const staged=directory && directory!==join(root,'public',prefix) ? directory : undefined;
 return {name:'male-ta98-package',enforce:'post',
  configureServer(server) {
   if (!staged) return;
   server.middlewares.use('/'+prefix,async(req,res,next)=>{
    if (req.method!=='GET' && req.method!=='HEAD') return next();
    const name=(req.url??'').split('?')[0];
    if (!/^\/(?:[a-zA-Z0-9._-]+\/)*[a-zA-Z0-9_-][a-zA-Z0-9._-]*$/.test(name)) return next();
    const file=resolve(staged,'.'+name);
    if (!file.startsWith(staged+sep)) return next();
    try {
     const info=await stat(file);
     if (!info.isFile()) return next();
     const mime:Record<string,string>={json:'application/json',png:'image/png',jpg:'image/jpeg',jpeg:'image/jpeg',webp:'image/webp'};
     res.setHeader('Content-Type',mime[name.split('.').pop()!]??'application/octet-stream');
     res.setHeader('Content-Length',info.size);
     res.setHeader('Cache-Control',name.endsWith('.json')?'no-cache':'public, max-age=31536000, immutable');
     if (req.method==='HEAD') return res.end();
     const stream=createReadStream(file);
     res.once('close',()=>stream.destroy());stream.on('error',error=>res.destroy(error));stream.pipe(res);
    } catch (error) { if ((error as NodeJS.ErrnoException).code==='ENOENT') return next(); next(error as Error); }
   });
  },
  async generateBundle(_,bundle) {
   if (!staged) return;
   const manifest=bundle['asset-manifest.json'];
   if (manifest?.type!=='asset') throw Error('Missing asset version manifest');
   const versions=JSON.parse(String(manifest.source));
   for (const entry of await readdir(staged,{recursive:true,withFileTypes:true})) {
    if (!entry.isFile()) continue;
    const file=join(entry.parentPath,entry.name),relative=file.slice(staged.length+1).split(sep).join('/');
    if (!/\.(json|bin|gz|png|jpe?g|webp)$/.test(relative)) continue;
    const fileName=prefix+'/'+relative,bytes=await readFile(file);
    this.emitFile({type:'asset',fileName,source:bytes});
    versions.files[fileName]={version:createHash('sha256').update(bytes).digest('hex'),bytes:bytes.length};
   }
   manifest.source=JSON.stringify(versions);
  },
  async writeBundle(options) {
   if (directory && options.dir) {
    const output=resolve(options.dir);
    if (output===resolve(root,'public') || output.startsWith(resolve(root,'public')+sep)) throw Error('Expected a generated site directory');
    await compactPagesGeometry(output);
   }
  },
 };
}
