import assert from 'node:assert/strict';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {StreamableHTTPClientTransport} from '@modelcontextprotocol/sdk/client/streamableHttp.js';
import {unpackRenderer} from '../mcp/packed-ui.mjs';
const endpoint=new URL(process.argv[2]??'https://human-atlas-connect.ctzurcanu.workers.dev/mcp');
const client=new Client({name:'human-atlas-live-check',version:'1.0.0'});
await client.connect(new StreamableHTTPClientTransport(endpoint));
try{
 const {tools}=await client.listTools();assert.equal(tools.length,3);
 for(const query of ['Stomach','heart','left kidney','TA98:','muscle']){
  const start=performance.now();
  const search=await client.callTool({name:'search_anatomy',arguments:{query}});
  assert.ok(!search.isError);assert.ok(search.structuredContent.matches.length>0);
  console.log(`Search ${query}: ${search.structuredContent.matches.length} matches, ${Math.round(performance.now()-start)} ms round trip`);
 }
 for(const model of ['male-detail','male-full','male','female','embryo','cell']){
  const options=await client.callTool({name:'get_anatomy_options',arguments:{model}});assert.ok(!options.isError);
  assert.ok(options.structuredContent.models.every(id=>!id.startsWith('local-')));
  const view=await client.callTool({name:'show_anatomy',arguments:{model}});assert.ok(!view.isError);
  assert.equal(new URL(view.structuredContent.url).origin,'https://ctzurcanu.github.io');
 }
 const stomach=await client.callTool({name:'show_anatomy',arguments:{structure:'Stomach',isolate:true}});
 assert.equal(stomach.structuredContent.structures[0].id,'ZA:Stomach');
 const bad=await client.callTool({name:'show_anatomy',arguments:{structure:'nonexistent-example-anatomy'}});assert.equal(bad.isError,true);
 const {resources}=await client.listResources();
 const ui=(await client.readResource({uri:resources[0].uri})).contents[0];
 assert.equal(ui.mimeType,'text/html;profile=mcp-app');assert.ok(ui.text.length<3*1024*1024);
 assert.ok(!/<script\b[^>]*\bsrc=/i.test(ui.text));assert.ok(!/<link\b/i.test(ui.text));
 const payload=await unpackRenderer(ui.text.match(/<script id="atlas-renderer" type="application\/gzip">([^<]+)<\/script>/)[1]);
 assert.ok(payload.script.includes('ui/initialize'));assert.ok(payload.styles.length>0);
 const original=(await client.readResource({uri:'ui://human-atlas/anatomy-view'})).contents[0];
 assert.equal(original.text,ui.text);assert.deepEqual(original._meta.ui.csp.frameDomains,[]);
 const v3=(await client.readResource({uri:'ui://human-atlas/anatomy-view-v3'})).contents[0];assert.equal(v3.text,ui.text);
 const previous=(await client.readResource({uri:'ui://human-atlas/anatomy-view-v2'})).contents[0];assert.equal(previous.text,ui.text);
 console.log(`Public MCP verified: ${endpoint}`);
}finally{await client.close();}
