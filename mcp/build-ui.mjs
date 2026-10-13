import {readFile,writeFile,mkdir} from 'node:fs/promises';
import {fileURLToPath,pathToFileURL} from 'node:url';
import {build} from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/postcss';
import {build as buildScript} from 'esbuild';
import {createHash} from 'node:crypto';
export const UI_FILE = new URL('../cloudflare/generated/mcp-view.html', import.meta.url);
export const IFRAME_FILE = new URL('../cloudflare/generated/mcp-iframe-view.html', import.meta.url);
export const APP_FILE = new URL('../cloudflare/generated/mcp-app.html', import.meta.url);
export const ASSETS_FILE = new URL('../cloudflare/generated/mcp-app-assets.mjs', import.meta.url);
let building;
let buildingIframe;
export function buildIframeUi() {
 return buildingIframe ??= (async () => {
  const bundle = await buildScript({entryPoints: [fileURLToPath(new URL('./iframe-main.mjs', import.meta.url))],bundle: true,platform: 'browser',format: 'iife',minify: true,write: false});
  const script = bundle.outputFiles[0].text.replace(/<\/script/gi, '<\\/script');
  const html = (await readFile(new URL('./iframe-view.html', import.meta.url),'utf8')).replace('<!-- ATLAS_SCRIPT -->', () => `<script>${script}</script>`);
  await mkdir(new URL('./', IFRAME_FILE), {recursive: true});
  await writeFile(IFRAME_FILE, html);
  return html;
 })();
}
export function buildMcpUi() {
 return building ??= (async () => {
  const path = value => fileURLToPath(new URL(value, import.meta.url));
  const bundle = await build({
   configFile: false, logLevel: 'silent', root: path('../'), publicDir: false,
   base: 'https://ctzurcanu.github.io/human-atlas/',
   define: {'process.env.NODE_ENV': JSON.stringify('production')},
   plugins: [react()], resolve: {alias: {'@': path('../')}},
   css: {postcss: {plugins: [tailwindcss()]}},
   build: {write: false, emptyOutDir: false, outDir: path('../cloudflare/generated/ui-build'), target: 'es2022', cssCodeSplit: false,
    lib: {entry: path('./main.tsx'), name: 'AtlasMcpApp', formats: ['iife'], fileName: 'atlas-app'},
   },
  });
  const output = (Array.isArray(bundle) ? bundle : [bundle]).flatMap(item => item.output);
  const script = output.find(item => item.type === 'chunk' && item.isEntry)?.code;
  const styles = output.filter(item => item.type === 'asset' && item.fileName.endsWith('.css')).map(item => item.source).join('\n');
  if (!script || !styles) throw new Error('MCP app build is missing its renderer or styles.');
  if (output.some(item => item.type === 'chunk' && !item.isEntry)) throw new Error('MCP app must be a single self-contained bundle.');
  const html = (await readFile(new URL('./view.html', import.meta.url), 'utf8'))
   .replace('<!-- ATLAS_STYLE -->', () => `<style>${styles.replace(/<\/style/gi, '<\\/style')}</style>`)
   .replace('<!-- ATLAS_SCRIPT -->', () => `<script>${script.replace(/<\/script/gi, '<\\/script')}</script>`);
  await mkdir(new URL('./', UI_FILE), {recursive: true});
  await writeFile(UI_FILE, html);
  // Keep the MCP resource small; serve the unchanged renderer as normal web assets.
  const origin='https://human-atlas-connect.ctzurcanu.workers.dev';
  const hash=value=>createHash('sha256').update(value).digest('hex').slice(0,16);
  const scriptPath=`/mcp-app/viewer-${hash(script)}.js`,stylePath=`/mcp-app/viewer-${hash(styles)}.css`;
  const app=(await readFile(new URL('./view.html',import.meta.url),'utf8'))
   .replace('<!-- ATLAS_STYLE -->',()=>`<link rel="stylesheet" href="${origin}${stylePath}">`)
   .replace('<!-- ATLAS_SCRIPT -->',()=>`<script src="${origin}${scriptPath}" defer></script>`);
  await writeFile(APP_FILE,app);
  await writeFile(ASSETS_FILE,`export default ${JSON.stringify({[scriptPath]:{type:'application/javascript; charset=utf-8',body:script},[stylePath]:{type:'text/css; charset=utf-8',body:styles}})};\n`);
  return html;
 })();
}
if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
 const html = await buildMcpUi();
 console.log(`MCP atlas app: ${(Buffer.byteLength(html)/1024/1024).toFixed(2)} MiB, no nested frames`);
}
