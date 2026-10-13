import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
const compiled=await build({entryPoints:['app/theme.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {followSystemTheme}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
test('System appearance switches existing renderers and capture theme without rebuilding the model',()=>{
 let listener,metaColor,classDark;
 const root={dataset:{theme:'dark'},classList:{toggle:(_name,value)=>{classDark=value;}}};
 const document={documentElement:root,querySelector:()=>({setAttribute:(_name,value)=>{metaColor=value;}})};
 const media={matches:false,addEventListener:(_name,fn)=>{listener=fn;},removeEventListener:(_name,fn)=>{assert.equal(fn,listener);listener=null;}};
 const stop=followSystemTheme(document,media);
 assert.equal(root.dataset.theme,'light');assert.equal(classDark,false);assert.equal(metaColor,'#f2f3f3');
 media.matches=true;listener();assert.equal(root.dataset.theme,'dark');assert.equal(classDark,true);assert.equal(metaColor,'#141b23');
 media.matches=false;listener();assert.equal(root.dataset.theme,'light');assert.equal(classDark,false);
 stop();assert.equal(listener,null);
});
