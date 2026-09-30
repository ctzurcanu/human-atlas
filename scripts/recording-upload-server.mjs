import {createServer} from 'node:http';
import {writeFile,mkdir} from 'node:fs/promises';
import {resolve} from 'node:path';

const output=resolve('outputs/atlas-demo-dark.webm');
await mkdir(resolve('outputs'),{recursive:true});
createServer(async (request,response)=>{
  response.setHeader('Access-Control-Allow-Origin','http://localhost:3016');
  response.setHeader('Access-Control-Allow-Methods','POST, OPTIONS');
  response.setHeader('Access-Control-Allow-Headers','Content-Type');
  if(request.method==='OPTIONS'){response.writeHead(204);response.end();return;}
  if(request.method!=='POST'||request.url!=='/record'){response.writeHead(404);response.end('Not found');return;}
  const chunks=[];for await(const chunk of request)chunks.push(chunk);
  await writeFile(output,Buffer.concat(chunks));
  response.writeHead(201,{'Content-Type':'text/plain'});response.end(output);
}).listen(3120,'127.0.0.1',()=>console.log(`Recording upload server: ${output}`));
