import test from 'node:test';
import assert from 'node:assert/strict';
import {existsSync, readFileSync} from 'node:fs';
import {gunzipSync} from 'node:zlib';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const reportPath='reports/ta98-female/installed-gallbladder-portal-course.json';
const built=await build({stdin:{contents:"export * from './app/viewer-state'; export * from './app/anatomical-relations'; export * from './app/requested-selection-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const read=path=>JSON.parse(readFileSync(path));

test('all four local female catalogues deliver one native portal course without extra draws, textures or anatomy',{skip:!existsSync(reportPath)},()=>{
 const report=read(reportPath);assert.equal(report.catalogues.length,4);
 assert.equal(new Set(report.catalogues.map(r=>r.binaryUrl)).size,1);
 for(const row of report.catalogues){
  const atlas=read(row.path),before=JSON.parse(gunzipSync(readFileSync(row.snapshot)));
  const part=atlas.parts.find(p=>p.id===report.partId),source=before.parts.find(p=>p.id===report.partId);
  assert.equal(atlas.parts.filter(p=>p.id===report.partId&&!p.suppressed).length,1);
  assert.deepEqual(api.resolveSelection(atlas,[part.conceptId]),[part.id]);
  // Later organ deliveries may add their own records. Bind this repair to
  // the literal portal record and bank, rather than freezing the whole atlas.
  const withoutBufferOffsets=p=>Object.fromEntries(Object.entries(p).filter(([k])=>!['chunk','indices','normals','positions','uvs'].includes(k)));
  assert.deepEqual(withoutBufferOffsets(part),withoutBufferOffsets(source));
  assert.deepEqual(atlas.materials[part.material],before.materials[source.material]);
  assert.deepEqual(part.bounds,source.bounds);
  assert.equal(part.vertexCount,source.vertexCount);assert.equal(part.indexCount,source.indexCount);
  assert.equal(atlas.portalGallbladderCourse.wholeSourceRepairComplete,false);
  const chunk=atlas.chunks[part.chunk],raw=readFileSync(chunk.url.replace('/local-models/','.local-models/'));
  assert.equal(createHash('sha256').update(raw).digest('hex'),row.binarySha256);
  assert.deepEqual(gunzipSync(readFileSync(chunk.gzip.replace('/local-models/','.local-models/'))),raw);
  assert.deepEqual(gunzipSync(readFileSync(row.path+'.gz')),readFileSync(row.path));
 }
});

test('portal correction keeps reciprocal gallbladder venous relations in the native frame',{skip:!existsSync(reportPath)},()=>{
 const a=read('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json');
 const source='LOCAL:female:gallbladder:1030',ids=api.resolveSelection(a,[source]),choice=api.requestedSelectionChoice(a,[source],ids);
 const gall=ids.map(id=>a.parts.find(p=>p.id===id)),portal=a.parts.find(p=>p.id==='LOCAL:female:portal_vein:1795');
 const links=api.anatomicalConceptFunctionalRelations(choice,gall,a.parts).filter(r=>r.target.id===portal.id);
 assert.ok(links.length);
 assert.ok(api.anatomicalRelations(portal,a.parts).some(r=>r.targetConcept?.id===choice.id&&r.targetConcept.elements.length===ids.length));
 const liver=a.parts.find(p=>p.id==='LOCAL:female:liver:1032');
 assert.ok(api.anatomicalRelations(portal,a.parts).some(r=>r.kind==='after'&&r.target.id===liver.id));
 assert.ok(api.anatomicalRelations(liver,a.parts).some(r=>r.kind==='before'&&r.target.id===portal.id));
 assert.ok(!api.anatomicalRelations(portal,a.parts).some(r=>r.kind==='after'&&r.target.id==='LOCAL:female:hepatic_veins:1794'));
 assert.deepEqual(api.resolveSelection(a,[portal.id]),[portal.id]);
});
