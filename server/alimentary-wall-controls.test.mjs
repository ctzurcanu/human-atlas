import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export {resolveGuestHierarchy,parseGuestHierarchy} from './app/guest-hierarchy';export {bowelWallPeelSteps,wallPeelDepth} from './app/bowel-wall-peeling';export {localTa98LayerCatalogue} from './app/local-ta98-layer-catalogue';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const hierarchy=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
for(const sex of ['female','male'])test(`${sex}: actual intestinal hierarchy groups peel every stored tissue and preserve source coverage`,()=>{
 const atlas=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 const nodes=new Map();function walk(n){nodes.set(n.id,n);n.children.forEach(walk)}api.resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 for(const [organ,bank,count] of [['small-intestine','BOWELWALL',36],['large-intestine','CAECOWALL',60]]){
  const node=nodes.get(`ALIMENTARY:${organ}:wall-assembly`);assert.ok(node);
  assert.equal(node.parts.length,count);assert.ok(node.parts.every(p=>p.id.includes(`:${bank}:`)&&p.sectionSolid&&p.pairedTissue&&!p.suppressed));
  const choice={id:`hierarchy:guest:ta98:${node.id}`,name:node.name,elements:node.parts.map(p=>p.id)};
  const steps=api.bowelWallPeelSteps(atlas,choice);assert.equal(steps?.length,8);
  assert.equal(new Set(steps.flatMap(s=>s.ids)).size,count);
  for(let stage=0;stage<=8;stage++){const hidden=steps.slice(0,stage).flatMap(s=>s.ids);assert.equal(api.wallPeelDepth(steps,hidden),stage)}
  assert.equal(api.bowelWallPeelSteps(atlas,{...choice,elements:choice.elements.slice(1)}),null,'partial groups cannot masquerade as complete walls');
 }
 assert.ok(!nodes.get('ALIMENTARY:large-intestine:wall-assembly').parts.some(p=>/rectum|valve/i.test(p.name)),'missing rectal/valve layers are not fabricated');
 for(const code of ['A05.6.01.002','A05.6.01.004','A05.6.01.005','A05.6.01.006','A05.6.01.008','A05.6.01.009','A05.6.01.010','A05.7.01.002','A05.7.01.004','A05.7.01.005','A05.7.01.006','A05.7.01.007','A05.7.03.012','A05.7.03.017']){
  assert.ok(nodes.get(code)?.directParts.length,`${code}: generated tissues must also be accessible on the canonical TA98 coat node`);
  assert.ok(nodes.get(code).directParts.every(p=>p.pairedTissue&&p.sectionSolid),'a source exterior must not be inferred as a finite coat');
 }
});
test('staged female rectum has one complete assembly and eight distinct peel stages in TA98',()=>{
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('bowel=rectal'),true),'/local-models/ta98-review/rectal-wall-main-atlas/female/atlas.json');
 assert.equal(api.localTa98LayerCatalogue('local-female-ta98',new URLSearchParams('bowel=rectal'),false),undefined);
 assert.equal(api.localTa98LayerCatalogue('local-ta98',new URLSearchParams('bowel=rectal'),true),'/local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json');
 const atlas=JSON.parse(fs.readFileSync('.local-models/ta98-review/rectal-wall-main-atlas/female/atlas.json'));
 const nodes=new Map();function walk(n){nodes.set(n.id,n);n.children.forEach(walk)}api.resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 const wall=nodes.get('CAECOWALL:female:rectum:wall-assembly');assert.equal(wall?.parts.length,8);
 assert.ok(atlas.parts.find(p=>p.id==='LOCAL:female:large_intestine_rectum:80')?.suppressed);
 const steps=api.bowelWallPeelSteps(atlas,{id:'hierarchy:guest:ta98:'+wall.id,name:wall.name,elements:wall.parts.map(p=>p.id)});
 assert.equal(steps?.length,8);assert.equal(new Set(steps.flatMap(s=>s.ids)).size,8);
 assert.match(steps[0].name,/distribution unresolved/);assert.match(steps[1].name,/allocation unresolved/);
 for(const [key,count] of [['mucosa',3],['muscularis',2]]){
  const node=nodes.get('CAECOWALL:female:rectum:'+key);assert.equal(node.parts.length,count);
  assert.equal(api.bowelWallPeelSteps(atlas,{id:'hierarchy:guest:ta98:'+node.id,name:node.name,elements:node.parts.map(p=>p.id)})?.length,count);
 }
 const all=nodes.get('ALIMENTARY:large-intestine:wall-assembly');assert.equal(all.parts.length,68);
 assert.ok(all.children.some(c=>c.id===wall.id),'rectum controls belong to the large-intestine wall group');
 assert.equal(nodes.get('A05.7.04.010').directParts.filter(p=>p.caecocolicRegion==='rectum').length,1);
 assert.equal(nodes.get('A05.7.04.015').directParts.filter(p=>p.caecocolicRegion==='rectum').length,1);
});
