import fs from 'node:fs';
import path from 'node:path';
import {builtinModules,createRequire} from 'node:module';
import {build as esbuild} from 'esbuild';

const root=path.resolve(import.meta.dirname,'..');
const builtins=new Set(builtinModules.map(name=>name.replace(/^node:/,'')));
const extensions=['.tsx','.ts','.jsx','.js','.mjs','.mts','.json'];
const config=JSON.parse(fs.readFileSync(path.join(root,'tsconfig.json'),'utf8'));

// Native esbuild directory resolution can stall on this large local checkout.
// Resolve and read test inputs through Node; esbuild still compiles and bundles
// the actual application files, with the same project compiler settings.
export function build(options){
 const files={name:'atlas-node-file-resolution',setup(builder){
  builder.onResolve({filter:/.*/},args=>{
   if(builtins.has(args.path.replace(/^node:/,'')))return {path:args.path,external:true};
   const requested=args.path.startsWith('@/')?path.join(root,args.path.slice(2)):args.path;
   if(requested.startsWith('.')||path.isAbsolute(requested)){
    const base=path.resolve(args.resolveDir||root,requested);
    const choices=[base,...extensions.map(ext=>base+ext),...extensions.map(ext=>path.join(base,'index'+ext))];
    for(const file of choices){try{if(fs.statSync(file).isFile())return {path:file};}catch(error){if(!['ENOENT','ENOTDIR'].includes(error.code))throw error;}}
    throw new Error('Unresolved local test import: '+base);
   }
   return {path:createRequire(path.join(args.resolveDir||root,'package.json')).resolve(requested)};
  });
  builder.onLoad({filter:/.*/},args=>{
   const ext=path.extname(args.path).slice(1);
   const loader=ext==='mts'?'ts':['ts','tsx','jsx','json','css'].includes(ext)?ext:'js';
   return {contents:fs.readFileSync(args.path,'utf8'),loader,resolveDir:path.dirname(args.path)};
  });
 }};
 return esbuild({...options,tsconfigRaw:options.tsconfigRaw??config,plugins:[...(options.plugins??[]),files]});
}
