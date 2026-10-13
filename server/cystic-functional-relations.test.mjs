import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/cystic-functional-relations';export * from './app/anatomical-relations';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
for(const sex of ['male','female'])test(`${sex} cystic duct exposes scoped upstream arterial and neural references`,()=>{
 const a=JSON.parse(fs.readFileSync(`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`));
 const duct=a.parts.find(p=>p.id===(sex==='male'?'FJ3080':'LOCAL:female:cystic_duct:authored'));
 const graph=api.cysticFunctionalRelations(a.parts),links=graph.get(duct.id);
 assert.ok(links.some(r=>r.kind==='arterial'&&r.via.includes('Arterial branch to cystic duct')));
 assert.ok(links.some(r=>r.kind==='innervation'&&r.via.includes('Hepatic plexus')));
 assert.equal(links.filter(r=>r.kind==='innervation').length,4);
 assert.ok(links.filter(r=>r.target.name.startsWith('Sympathetic nerves')).every(r=>r.via.includes('Coeliac plexus')&&r.note.endsWith('duct-wall fibers are not modeled.')));
 assert.ok(links.every(r=>!r.target.id.startsWith('HRA:')));
 const publicLinks=api.anatomicalRelations(duct,a.parts);
 for(const r of links){
  assert.ok(publicLinks.some(x=>x.kind===r.kind&&x.target.id===r.target.id));
  assert.ok(api.anatomicalRelations(r.target,a.parts).some(x=>x.target.id===duct.id));
  if(r.via)assert.equal(r.viaModeled,false);
  for(const change of [{sourceOffset:[.02,0,0]},{suppressed:true},{provenance:{mirrored:true}}]){
   const changed=a.parts.map(p=>p===r.target?{...p,...change}:p);
   assert.ok(!(api.cysticFunctionalRelations(changed).get(duct.id)??[]).some(x=>x.target.id===r.target.id));
  }
 }
 if(sex==='male')assert.ok(links.some(r=>r.kind==='adjacent'&&r.target.name==='Cystic node'));
 for(const r of links)assert.ok(graph.get(r.target.id).some(x=>x.target.id===duct.id));
 for(const change of [{suppressed:true},{sourceId:'foreign'},{conceptId:'wrong'},{sourceOffset:[0,0,0]},{pairedTissue:false},{sectionSolid:false},{tissueDrawGroups:[]}]){
  const altered=a.parts.map(p=>p===duct?{...p,...change}:p);assert.equal(api.cysticFunctionalRelations(altered).has(duct.id),false);
 }
});
