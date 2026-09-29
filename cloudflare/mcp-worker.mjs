import {WebStandardStreamableHTTPServerTransport} from '@modelcontextprotocol/sdk/server/webStandardStreamableHttp.js';
import {configureAtlasData} from '../mcp/atlas-core.mjs';
import {createServer} from '../mcp/server-core.mjs';
import data from './generated/mcp-data.json';
import uiHtml from '../mcp/view.html';

configureAtlasData((_model,file)=>{const atlas=data.catalogues[file];if(!atlas)throw new Error('Model is not publicly available.');return atlas;},data.terminology);
export async function handleMcp(request,env){
 const origin=request.headers.get('Origin');
 const allowed=new Set([new URL(request.url).origin,'https://chatgpt.com','https://chat.openai.com',...(env.ALLOWED_VIEWER_ORIGINS??'').split(',')]);
 if(origin&&!allowed.has(origin))return new Response('Forbidden',{status:403});
 const headers=origin?{'Access-Control-Allow-Origin':origin,'Vary':'Origin'}:{};
 if(request.method==='OPTIONS')return new Response(null,{status:204,headers:{...headers,'Access-Control-Allow-Methods':'POST, GET, DELETE, OPTIONS','Access-Control-Allow-Headers':'Content-Type, Accept, MCP-Protocol-Version, MCP-Session-Id'}});
 if(Number(request.headers.get('Content-Length')??0)>65536)return new Response('Request too large',{status:413});
 const server=createServer({uiHtml,publicOnly:true});
 const transport=new WebStandardStreamableHTTPServerTransport({enableJsonResponse:true});
 await server.connect(transport);
 try{
  let parsedBody;
  if(request.method==='POST'){
   const text=await request.text();
   if(text.length>65536)return new Response('Request too large',{status:413});
   try{parsedBody=JSON.parse(text);}catch{return Response.json({jsonrpc:'2.0',id:null,error:{code:-32700,message:'Parse error'}},{status:400,headers});}
  }
  const result=await transport.handleRequest(request,{parsedBody});
  const response=new Response(result.body,result);
  for(const [key,value]of Object.entries(headers))response.headers.set(key,value);
  return response;
 }finally{await server.close();}
}
