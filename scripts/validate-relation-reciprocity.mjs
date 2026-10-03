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
// An unrelated return link (e.g. opposite side) cannot prove that a route,
// supply, component or attachment relationship has its proper inverse.
const inverseKinds={
 before:['after'],after:['before'],innervation:['innervates'],innervates:['innervation'],
 arterial:['supplies'],supplies:['arterial'],venous:['drains'],drains:['venous'],
 articulates:['articulates'],connects:['connectedBy'],connectedBy:['connects'],joint:['joint'],
 continuous:['continuous'],covers:['coveredBy'],coveredBy:['covers'],partOf:['contains'],contains:['partOf'],
 originFor:['originSites'],originSites:['originFor'],insertionFor:['insertionSites'],insertionSites:['insertionFor'],
 counterpart:['counterpart'],adjacent:['adjacent'],lymphaticDrainage:['lymphaticTributaries'],lymphaticTributaries:['lymphaticDrainage'],origin:['muscles'],insertion:['muscles'],
 cartilages:['bones'],bones:['cartilages','tendons','fascia'],
 tendons:['muscles','bones','fascia'],muscles:['tendons','origin','insertion'],fascia:['tendons','bones'],
};

for(const file of catalogs){
 const parts=JSON.parse(readFileSync(file,'utf8')).parts.filter(part=>!part.suppressed);
 const relations=new Map(parts.map(part=>[part.id,anatomicalRelations(part,parts)]));
 const missing=[];
 let count=0;
 for(const source of parts)for(const relation of relations.get(source.id)){
  count++;
  assert.ok(inverseKinds[relation.kind],`Unspecified inverse kind: ${relation.kind}`);
  // The muscle inspector calls its attached tendon "tendons", while the
  // tendon inspector calls the same physical attachment "connects".
  const expected=relation.kind==='connects'&&relation.target.system==='muscular'?['connectedBy','tendons']:inverseKinds[relation.kind];
  if(relations.get(relation.target.id)?.some(back=>back.target.id===source.id&&expected.includes(back.kind)))continue;
  missing.push(`${source.name} [${source.id}] --${relation.kind}--> ${relation.target.name} [${relation.target.id}]`);
 }
 assert.equal(missing.length,0,`${file} has ${missing.length} one-way relationship links:\n${missing.slice(0,20).join('\n')}`);
 console.log(`${file}: ${parts.length} pieces, ${count} reciprocal links`);
}
