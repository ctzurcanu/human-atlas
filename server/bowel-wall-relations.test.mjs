import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {build} from 'esbuild';

const built=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const functional=new Set(['arterial','venous','innervation','supplies','drains','innervates']);
for(const sex of ['female','male'])test(`${sex}: finite bowel tissues inherit native organ routes with intermediate structures`,async()=>{
 const atlas=JSON.parse(await readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`,'utf8'));
 const walls=atlas.parts.filter(p=>p.bowelWallSource&&/:BOWELWALL:(jejunum|ileum):/.test(p.id));
 assert.equal(walls.length,24);
 for(const wall of walls){
  const links=anatomicalRelations(wall,atlas.parts).filter(r=>functional.has(r.kind));
  // The male ileum is a separate LOCAL source; ZA collectors cannot establish registration.
  if(sex==='male'&&wall.id.includes(':ileum:'))assert.ok(!links.some(r=>r.target.id.startsWith('ZA:')));
  else{
   assert.ok(links.some(r=>r.kind==='arterial'&&r.via.includes(wall.id.includes(':jejunum:')?'Jejunal arteries':'Ileal arteries')));
   assert.ok(links.some(r=>r.kind==='venous'));
   assert.ok(links.some(r=>r.kind==='innervation'&&r.via.includes('Enteric plexuses')));
   for(const link of links){
    assert.match(link.note,/Organ-level association/);assert.match(link.note,/unverified/);
    assert.ok(anatomicalRelations(link.target,atlas.parts).some(r=>r.target.id===wall.id&&r.note===link.note));
   }
  }
 }
 const wall=walls.find(p=>p.id.includes(':jejunum:'));
 for(const invalid of [{...wall,bowelWallSource:undefined},{...wall,bowelWallSource:'other'},{...wall,sourceId:'wrong'},{...wall,id:'HRA:donor-wall'},{...wall,conceptId:'SIW:ileum:serosa'},{...wall,id:wall.id.split(':BOWELWALL:')[0]+':BOWELWALL:jejunum:circular-folds:serosa',conceptId:'SIW:jejunum:circular-folds'}]){
  const parts=atlas.parts.map(p=>p.id===wall.id?invalid:p);
  assert.equal(anatomicalRelations(invalid,parts).filter(r=>functional.has(r.kind)).length,0);
 }
});
test('female duodenal tissues inherit pancreaticoduodenal supply and the biliary route without inventing a papillary partition',async()=>{
 const atlas=JSON.parse(await readFile('.local-models/ta98-review/bowel-wall-main-prototype/female/atlas.json'));
 const walls=atlas.parts.filter(p=>p.bowelWallSource&&p.id.includes(':BOWELWALL:duodenum:'));assert.equal(walls.length,12);
 const routes=[['LOCAL:female:anterior_superior_pancreaticoduodenal_artery:1678','arterial','supplies','Duodenal branches'],['LOCAL:female:gastroduodenal_artery:1723','arterial','supplies','Superior pancreaticoduodenal arteries'],['LOCAL:female:superior_mesenteric_artery:1670','arterial','supplies','Inferior pancreaticoduodenal artery'],['LOCAL:female:anterior_inferior_pancreaticoduodenal_vein:1758','venous','drains','Duodenal venules'],['LOCAL:female:superior_mesenteric_vein:1755','venous','drains','Pancreaticoduodenal veins'],['LOCAL:female:l_vagus_nerve:894','innervation','innervates','Enteric plexuses'],['LOCAL:female:bile_duct:1031','before','after','Major duodenal papilla']];
 for(const wall of walls)for(const [id,kind,reverse,via] of routes){const target=atlas.parts.find(p=>p.id===id),links=anatomicalRelations(wall,atlas.parts),link=links.find(r=>r.kind===kind&&r.target.id===id);assert.ok(link,`${wall.id} → ${id}`);assert.ok(link.via.includes(via));assert.match(link.note,/Organ-level/);assert.match(link.note,/unverified/);assert.ok(anatomicalRelations(target,atlas.parts).some(r=>r.kind===reverse&&r.target.id===wall.id&&r.note===link.note));}
 const wall=walls.find(p=>p.id.endsWith(':epithelium'));
 for(const invalid of [{...wall,bowelWallSource:undefined},{...wall,sourceId:'wrong'},{...wall,id:'HRA:unregistered-duodenal-layer'},{...wall,conceptId:'SIW:jejunum:epithelium'}]){const parts=atlas.parts.map(p=>p.id===wall.id?invalid:p);const links=anatomicalRelations(invalid,parts);assert.ok(!links.some(r=>functional.has(r.kind)||r.target.id==='LOCAL:female:bile_duct:1031'));}
});
test('male duodenal preview has no functional inheritance from an absent ZA source to a LOCAL carrier',async()=>{
 const atlas=JSON.parse(await readFile('.local-models/ta98-review/bowel-wall-main-prototype/male/atlas.json'));
 const walls=atlas.parts.filter(p=>p.bowelWallSource&&p.id.includes(':BOWELWALL:duodenum:'));assert.equal(walls.length,12);assert.ok(!atlas.parts.some(p=>p.id==='ZA:Duodenum'));
 for(const wall of walls){assert.equal(wall.bowelWallSource,'LOCAL:male:stomach_duodenum:1265');assert.equal(wall.sourceId,'stomach_duodenum');assert.match(wall.coverageLimitation,/separate LOCAL/);assert.ok(!anatomicalRelations(wall,atlas.parts).some(r=>functional.has(r.kind)));}
});
