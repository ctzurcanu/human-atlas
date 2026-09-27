import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';

const models=process.argv.slice(2);
const catalogs=models.length?models:[
 'public/models/atlas-male-complete.json',
 'public/models/atlas.json',
 'public/models/atlas-hra-female.json',
 'public/models/atlas-embryo.json',
 'public/models/atlas-cell.json',
 '.local-models/reference.json',
 '.local-models/male.json',
 '.local-models/female.json',
];
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);

for(const file of catalogs){
 const parts=JSON.parse(readFileSync(file,'utf8')).parts.filter(part=>!part.suppressed);
 const relations=new Map(parts.map(part=>[part.id,anatomicalRelations(part,parts)]));
 const missing=[];
 let count=0;
 for(const source of parts)for(const relation of relations.get(source.id)){
  count++;
  if(relations.get(relation.target.id)?.some(back=>back.target.id===source.id))continue;
  missing.push(`${source.name} [${source.id}] --${relation.kind}--> ${relation.target.name} [${relation.target.id}]`);
 }
 assert.equal(missing.length,0,`${file} has ${missing.length} one-way relationship links:\n${missing.slice(0,20).join('\n')}`);
 console.log(`${file}: ${parts.length} pieces, ${count} reciprocal links`);
}
