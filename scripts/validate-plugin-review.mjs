import assert from 'node:assert/strict';
import {readFile,writeFile} from 'node:fs/promises';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {StreamableHTTPClientTransport} from '@modelcontextprotocol/sdk/client/streamableHttp.js';
const client=new Client({name:'human-atlas-review-check',version:'1.0.0'});
const cases=JSON.parse(await readFile('plugin-review/test-cases.json','utf8'));
const evidence={scope:'Direct MCP tool checks; conversational tool selection and all eight prompts still need manual ChatGPT evaluation.',endpoint:'https://human-atlas-connect.ctzurcanu.workers.dev/mcp',checks:[]};
await client.connect(new StreamableHTTPClientTransport(new URL(evidence.endpoint)));
try{
 const {tools}=await client.listTools();
 assert.deepEqual(tools.map(t=>t.name).sort(),['get_anatomy_options','search_anatomy','show_anatomy']);
 for(const tool of tools)for(const field of ['readOnlyHint','destructiveHint','openWorldHint'])assert.equal(typeof tool.annotations[field],'boolean');
 async function call(index,name,args,check){
  const result=await client.callTool({name,arguments:args});check(result);
  evidence.checks.push({case:cases.positive[index].description,tool:name,arguments:args,status:'passed',output:result.structuredContent});
 }
 await call(0,'show_anatomy',{structure:'Stomach',isolate:true},r=>{assert.ok(!r.isError);const u=new URL(r.structuredContent.url);assert.equal(u.searchParams.get('isolate'),'1');assert.equal(u.searchParams.get('focus'),'1');assert.equal(r.structuredContent.structures[0].id,'ZA:Stomach');});
 await call(1,'get_anatomy_options',{},r=>{assert.ok(!r.isError);assert.deepEqual(r.structuredContent.models,['male-detail','male-full','male','female','embryo','cell']);});
 await call(2,'show_anatomy',{model:'female',systems:['skeletal'],view:'front'},r=>{assert.ok(!r.isError);const u=new URL(r.structuredContent.url);assert.equal(u.searchParams.get('model'),'female');assert.equal(u.searchParams.get('layers'),'skeletal');assert.equal(u.searchParams.get('view'),'front');});
 await call(3,'search_anatomy',{query:'left kidney',model:'male-detail'},r=>{assert.ok(!r.isError);assert.ok(r.structuredContent.matches.some(m=>/left/i.test(m.name)));});
 await call(4,'show_anatomy',{structures:['Stomach','Pancreas'],isolate:true},r=>{assert.ok(!r.isError);assert.equal(r.structuredContent.structures.length,2);assert.equal(new URL(r.structuredContent.url).searchParams.get('isolate'),'1');});
 for(const [index,args]of [[0,{structure:'nonexistent-example-anatomy',isolate:true}],[1,{model:'local-male'}]]){
  const result=await client.callTool({name:'show_anatomy',arguments:args});assert.equal(result.isError,true);
  evidence.checks.push({case:cases.negative[index].description,tool:'show_anatomy',arguments:args,status:'passed',output:result.content});
 }
 evidence.checks.push({case:cases.negative[2].description,status:'capability check passed; conversational scope response pending',evidence:'The only exposed tools search the public anatomy catalogue, list options and generate views. No patient-record or diagnosis tool is exposed.'});
 await writeFile('plugin-review/tool-checks.json',JSON.stringify(evidence,null,2)+'\n');
 console.log('Five positive tool checks and two negative input checks passed; no patient-record tool is exposed. Evidence: plugin-review/tool-checks.json');
}finally{await client.close();}
