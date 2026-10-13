import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
import {existsSync,readFileSync} from 'node:fs';
const compiled=await build({stdin:{contents:"export * from './app/bowel-wall-peeling';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {bowelWallPeelSteps,wallPeelDepth}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const layers=['external-connective-coat','longitudinal-muscle','circular-muscle','submucosa','muscularis-mucosae','lamina-propria','epithelium'];
function fixture(){
 const parts=['native','connector'].flatMap(region=>layers.map(layer=>({id:`OESOPHWALL:male:${region}:${layer}`,name:layer,pairedTissue:true,bowelWallSource:region,tissueWallAssembly:'oesophageal:male',tissueWallLayer:layer,chunk:0,positions:0,normals:300,geometryOrigin:[0,1,0]})));
 return {parts,concepts:[],bowelWallPrototype:{partIds:parts.map(p=>p.id)}};
}
const choice=ids=>({id:'hierarchy:guest:ta98:A05.4.01.001',elements:ids});
test('whole oesophagus peels seven coats across both connected regions',()=>{
 const atlas=fixture(),steps=bowelWallPeelSteps(atlas,choice(atlas.parts.map(p=>p.id)));
 assert.equal(steps.length,7);assert.ok(steps.every(s=>s.ids.length===2));assert.equal(steps[0].name,'External connective coat');
 assert.equal(wallPeelDepth(steps,[...steps[0].ids,...steps[1].ids]),2);
 assert.equal(wallPeelDepth(steps,[steps[0].ids[0]]),0);
});
test('regional, muscular and mucosal selections retain their own peel scope',()=>{
 const atlas=fixture();
 for(const [predicate,n] of [[p=>p.bowelWallSource==='connector',7],[p=>['longitudinal-muscle','circular-muscle'].includes(p.tissueWallLayer),2],[p=>['muscularis-mucosae','lamina-propria','epithelium'].includes(p.tissueWallLayer),3]]){
  const ids=atlas.parts.filter(predicate).map(p=>p.id),steps=bowelWallPeelSteps(atlas,choice(ids));assert.equal(steps.length,n);assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(ids));
 }
});
test('partial coats, mixed banks and multi-selection do not masquerade as a complete wall',()=>{
 const atlas=fixture(),ids=atlas.parts.map(p=>p.id);
 assert.equal(bowelWallPeelSteps(atlas,choice(ids.slice(1))),null);
 atlas.parts[1].positions=12;assert.equal(bowelWallPeelSteps(atlas,choice(ids)),null);
 atlas.parts[1].positions=0;assert.equal(bowelWallPeelSteps(atlas,{...choice(ids),id:'selection-set'}),null);
});
const packagePath='.local-models/ta98-review/oesophageal-layer-bank/male/runtime-package/package.json';
test('actual packed cells expose seven regional coats and scoped aggregate controls',{skip:!existsSync(packagePath)},()=>{
 const packed=JSON.parse(readFileSync(packagePath,'utf8'));
 const atlas={parts:packed.parts,concepts:[],bowelWallPrototype:{partIds:packed.parts.map(p=>p.id)}};
 assert.equal(atlas.parts.length,14);
 for(const source of [null,'ZA:Oesophagus','GASTRIC:male:cardiac-connection']){
  const region=atlas.parts.filter(p=>!source||p.bowelWallSource===source);
  const steps=bowelWallPeelSteps(atlas,choice(region.map(p=>p.id)));
  assert.equal(steps.length,7);
  assert.deepEqual(steps.map(step=>step.ids.length),Array(7).fill(source?1:2));
  for(const [keys,count] of [[layers.slice(1,3),2],[layers.slice(4),3]]){
   const members=region.filter(p=>keys.includes(p.tissueWallLayer)).map(p=>p.id);
   const aggregate=bowelWallPeelSteps(atlas,choice(members));
   assert.equal(aggregate.length,count);
   assert.deepEqual(new Set(aggregate.flatMap(step=>step.ids)),new Set(members));
  }
 }
});
