import test from 'node:test';
import assert from 'node:assert/strict';
import {existsSync,readFileSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {build} from 'esbuild';
const reportPath='reports/installed-rectal-donor-models.json';
const compiled=await build({stdin:{contents:"export * from './app/viewer-state'; export * from './app/bowel-wall-peeling'; export * from './app/guest-hierarchy'; export * from './app/guest-choice'; export * from './app/caecocolic-wall-relations'; export * from './app/requested-selection-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
test('both ordinary local TA98 models replace one rectum with eight peelable hollow coats and reciprocal sigmoid rims',{skip:!existsSync(reportPath)},()=>{
 const report=JSON.parse(readFileSync(reportPath));
 for(const row of report.catalogues){
  const atlas=JSON.parse(readFileSync(row.path)),sex=atlas.rectalWallIntegration.sex;
  const ids=atlas.rectalWallIntegration.partIds,parts=ids.map(id=>atlas.parts.find(p=>p.id===id));assert.equal(parts.length,8);assert.ok(parts.every(p=>p&&!p.suppressed&&p.pairedTissue&&p.sectionSolid));
  const source=atlas.rectalWallIntegration.sourceRectum;assert.equal(atlas.parts.find(p=>p.id===source).suppressed,true);assert.deepEqual(new Set(api.resolveSelection(atlas,[source])),new Set(ids));
  const choice=api.requestedSelectionChoice(atlas,[source],ids);assert.ok(choice);assert.equal(choice.id,atlas.parts.find(p=>p.id===source).conceptId);assert.deepEqual(new Set(choice.elements),new Set(ids));assert.equal(api.bowelWallPeelSteps(atlas,choice).length,8);
  assert.equal(api.requestedSelectionChoice(atlas,[source],ids.slice(1)),null,'a partial wall cannot retain the whole-organ identity');
  assert.equal(api.requestedSelectionChoice(atlas,[source,'unrelated'],ids),null,'an explicit mixed request cannot become one whole organ');
  assert.ok(!api.resolveSelection(atlas,atlas.defaultHidden??[]).some(id=>ids.includes(id)),'retired hidden references must not hide their installed rectal replacement');
  assert.equal(api.bowelWallPeelSteps(atlas,{id:source,elements:ids}).length,8);
  const parent='CAECOWALL:'+sex+':rectum';for(const [key,count] of [['muscularis',2],['mucosa',3]]){const c=atlas.concepts.find(c=>c.id===parent+':'+key);assert.equal(api.bowelWallPeelSteps(atlas,{id:c.id,elements:c.elements}).length,count);}
  const tree=api.resolveGuestHierarchy(atlas,api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json')))),nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};tree.forEach(walk);const wall=nodes.get(parent+':wall-assembly');assert.ok(wall);assert.deepEqual(new Set(wall.parts.map(p=>p.id)),new Set(ids));assert.equal(wall.children.length,5);assert.equal(nodes.get(parent+':muscularis').children.length,2);assert.equal(nodes.get(parent+':mucosa').children.length,3);
  const bank=atlas.chunks[parts[0].chunk];const raw=readFileSync(bank.url.replace('/local-models/','.local-models/'));assert.equal(raw.length,bank.bytes);assert.deepEqual(gunzipSync(readFileSync(bank.gzip.replace('/local-models/','.local-models/'))),raw);
  const graph=api.caecocolicWallRelations(atlas.parts);
  for(const p of parts){assert.equal(p.positions,parts[0].positions);assert.equal(p.normals,parts[0].normals);assert.deepEqual(atlas.materials[p.material].tissueAppearance,p.tissueAppearance);assert.ok((graph.get(p.id)??[]).some(r=>r.target.caecocolicRegion==='sigmoid'));
   const index=new Uint32Array(raw.buffer.slice(raw.byteOffset+p.indices,raw.byteOffset+p.indices+p.indexCount*4));assert.ok(index.every(i=>i<p.vertexCount));
   const exterior=p.tissueDrawGroups.filter(g=>!g.neighbor);assert.ok(exterior.length); // physiological outlet has tissue annuli, not a lumen disk
   for(const g of p.tissueDrawGroups.filter(g=>g.neighbor?.includes(':sigmoid:'))){const other=atlas.parts.find(q=>q.id===g.neighbor),h=other.tissueDrawGroups.find(h=>h.neighbor===p.id);assert.equal(h.count,g.count);const a=index.slice(g.start,g.start+g.count),b=new Uint32Array(raw.buffer.slice(raw.byteOffset+other.indices+h.start*4,raw.byteOffset+other.indices+(h.start+h.count)*4));for(let i=0;i<a.length;i+=3)assert.deepEqual([a[i],a[i+2],a[i+1]],Array.from(b.slice(i,i+3)));}
  }
  if(sex==='female')assert.equal(atlas.parts.find(p=>p.id==='HRA:VH_F_rectum').suppressed,true);
 }
});
