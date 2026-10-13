// Local browser fixture: enforce Claude's no-nested-frame policy on the real app.
import {createServer} from 'node:http';
import {readFile} from 'node:fs/promises';
import {APP_FILE,IFRAME_FILE} from './build-ui.mjs';
import {anatomyView} from './atlas-data.mjs';

const host = `<!doctype html><html><head><meta charset="utf-8"><title>Atlas MCP sandbox check</title>
<style>body{font:15px system-ui;background:#eee;margin:12px}iframe{width:100%;height:610px;border:0}pre{white-space:pre-wrap}</style></head><body>
<label>Model <select id="model"><option>cell</option><option>male-detail</option><option>female</option><option>embryo</option></select></label>
<button id="show">Show model</button><button id="theme">Change host theme</button>
<iframe id="app" title="MCP Apps sandbox" sandbox="allow-scripts allow-same-origin allow-downloads"></iframe><pre id="checks"></pre>
<script>
const frame=document.getElementById('app'),checks=document.getElementById('checks');let theme='light',shown=false;
document.getElementById('model').value=new URLSearchParams(location.search).get('model')??'cell';
frame.src='/app'+location.search;
const send=data=>frame.contentWindow.postMessage({jsonrpc:'2.0',...data},'*');
async function show(){const params=new URLSearchParams(location.search);params.set('model',document.getElementById('model').value);const result=await(await fetch('/result?'+params)).json();send({method:'ui/notifications/tool-result',params:result});shown=true;}
addEventListener('message',event=>{if(event.source!==frame.contentWindow)return;const message=event.data;
if(message.method==='ui/initialize')send({id:message.id,result:{protocolVersion:'2026-01-26',hostInfo:{name:'Claude policy fixture',version:'1'},hostCapabilities:{},hostContext:{theme}}});
if(message.method==='ui/notifications/initialized')show();});
document.getElementById('show').onclick=show;
document.getElementById('theme').onclick=()=>{theme=theme==='light'?'dark':'light';send({method:'ui/notifications/host-context-changed',params:{theme}});};
setInterval(()=>{try{const d=frame.contentDocument;checks.textContent=JSON.stringify({shown,nestedFrames:d.querySelectorAll('iframe').length,canvas:d.querySelectorAll('canvas').length,theme:d.documentElement.dataset.theme,status:d.getElementById('status')?.hidden?'hidden':d.getElementById('status')?.textContent,alerts:[...d.querySelectorAll('[role="alert"]')].map(x=>x.textContent),loading:[...d.querySelectorAll('.loading')].map(x=>x.textContent),violations:frame.contentWindow.atlasViolations??[],errors:frame.contentWindow.atlasErrors??[]},null,2);}catch(e){checks.textContent=String(e);}},1000);
</script></body></html>`;
const csp="default-src 'none'; script-src 'self' 'unsafe-inline' https://human-atlas-connect.ctzurcanu.workers.dev; style-src 'self' 'unsafe-inline' https://human-atlas-connect.ctzurcanu.workers.dev; img-src 'self' data: blob: https://ctzurcanu.github.io; font-src 'self' data: https://ctzurcanu.github.io; connect-src https://ctzurcanu.github.io https://human-atlas-connect.ctzurcanu.workers.dev wss://human-atlas-connect.ctzurcanu.workers.dev; frame-src 'none'; base-uri 'none'; object-src 'none'";
const instrumentation=`<script>window.atlasViolations=[];window.atlasErrors=[];addEventListener('securitypolicyviolation',e=>atlasViolations.push({directive:e.violatedDirective,uri:e.blockedURI}));addEventListener('error',e=>atlasErrors.push(e.message));addEventListener('unhandledrejection',e=>atlasErrors.push(String(e.reason)));</script>`;
createServer(async(req,res)=>{
 try {
  const url=new URL(req.url,'http://localhost:3018');
  if(url.pathname==='/app') {const legacy=url.searchParams.has('legacy');res.setHeader('Content-Security-Policy',legacy?csp.replace("frame-src 'none'","frame-src https://ctzurcanu.github.io"):csp);res.setHeader('Content-Type','text/html');res.end((await readFile(legacy?IFRAME_FILE:APP_FILE,'utf8')).replace('<body>','<body>'+instrumentation));}
  else if(url.pathname==='/result') {const axis=url.searchParams.get('cut');res.setHeader('Content-Type','application/json');res.end(JSON.stringify({structuredContent:anatomyView({model:url.searchParams.get('model')??'cell',labels:false,...(['axial','sagittal','coronal'].includes(axis)?{section:{axis,position:.5,flip:false},view:'side',skinOpacity:.25}:{}),controls:['model','search','sections','systems','details','open']})}));}
  else if(url.pathname==='/') {res.setHeader('Content-Type','text/html');res.end(host);}
  else {res.statusCode=404;res.end();}
 } catch(error){res.statusCode=500;res.end(String(error));}
}).listen(3018,'127.0.0.1',()=>console.log('MCP sandbox fixture: http://localhost:3018'));
