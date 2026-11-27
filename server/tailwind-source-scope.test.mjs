import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,mkdirSync} from 'node:fs';
import {resolve,relative,isAbsolute,dirname} from 'node:path';
import {performance} from 'node:perf_hooks';
import test from 'node:test';
import postcss from 'postcss';
import tailwindcss from '@tailwindcss/postcss';

const root=process.cwd(),cssPath=resolve(root,'app/globals.css');
const sourceDirectories=['app','components','web'].map(dir=>resolve(root,dir));
const inside=(file,folder)=>{const rel=relative(folder,file);return rel===''||(!rel.startsWith('..')&&!isAbsolute(rel));};

test('actual UI stylesheet generates utilities without a workspace-wide source dependency',async()=>{
 const started=performance.now(),source=readFileSync(cssPath,'utf8');
 assert.match(source,/@import\s+['"]tailwindcss['"]\s+source\(none\)/);
 const result=await postcss([tailwindcss()]).process(source,{from:cssPath});
 assert.match(result.css,/\.inline-flex\s*\{/,'Shared button/component utilities must remain generated');
 assert.match(result.css,/\.flex\s*\{/,'Application flex utilities must remain generated');
 const directories=result.messages.filter(message=>message.type==='dir-dependency');
 assert.ok(directories.length>0,'Inspect actual Tailwind source dependencies, rather than only the configured strings');
 const scopes=directories.map(message=>{
  const pattern=message.glob??'**',index=[...pattern].findIndex(char=>'*?{['.includes(char));
  const prefix=index<0?pattern:pattern.slice(0,index),scope=resolve(message.dir,prefix);
  assert.ok(sourceDirectories.some(dir=>inside(scope,dir)),`Source dependency escapes UI directories: ${scope} (${pattern})`);
  return {dir:relative(root,resolve(message.dir)),glob:pattern,scope:relative(root,scope)};
 });
 const forbidden=['.local-models','reports','scripts','tests','tmp','outputs'].map(dir=>resolve(root,dir));
 for(const message of result.messages.filter(message=>message.type==='dependency')){
  assert.ok(!forbidden.some(dir=>inside(resolve(message.file),dir)),`Generated geometry/review data became a CSS dependency: ${message.file}`);
 }
 for(const dir of sourceDirectories)assert.ok(scopes.some(scope=>inside(resolve(root,scope.scope),dir)),`Missing UI source directory ${dir}`);
 const evidence={schema:'tailwind-ui-source-scope/v1',passed:true,elapsed_seconds:(performance.now()-started)/1000,
  css_bytes:Buffer.byteLength(result.css),source_dependency_scopes:scopes,
  dependency_messages:result.messages.filter(message=>message.type==='dependency').length,
  rss_bytes:process.memoryUsage().rss,process_peak_rss_bytes:process.resourceUsage().maxRSS*1024,
  utilities_preserved:['flex','inline-flex'],workspace_wide_globs:0};
 const destination=process.env.ATLAS_TAILWIND_SCOPE_EVIDENCE;
 if(destination){mkdirSync(dirname(resolve(destination)),{recursive:true});writeFileSync(destination,JSON.stringify(evidence,null,2)+'\n');}
 console.log(JSON.stringify(evidence));
});
