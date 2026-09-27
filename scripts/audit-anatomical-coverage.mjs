import {readFileSync} from 'node:fs';
import {build} from 'esbuild';

const file=process.argv[2]??'public/models/atlas-male-complete.json';
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);
const parts=JSON.parse(readFileSync(file,'utf8')).parts;
const sampleSize=Number(process.argv[3]??250);
const sampled=sampleSize>0&&sampleSize<parts.length?parts.filter((_,index)=>index%Math.ceil(parts.length/sampleSize)===0):parts;
const bySystem={};
for(const part of sampled){
 const entry=bySystem[part.system]??={total:0,linked:0,unresolved:[]};entry.total++;
 if(anatomicalRelations(part,parts).length)entry.linked++;
 else if(sampled===parts||entry.unresolved.length<10)entry.unresolved.push({id:part.id,name:part.name});
 bySystem[part.system]=entry;
}
const total=sampled.length,linked=Object.values(bySystem).reduce((count,entry)=>count+entry.linked,0);
console.log(JSON.stringify({file,atlasPieces:parts.length,auditedPieces:total,linked,unresolved:total-linked,bySystem},null,2));
