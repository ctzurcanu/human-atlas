import test from 'node:test';
import assert from 'node:assert/strict';
import {existsSync,readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
import {build} from 'esbuild';
const currentReport='reports/installed-oesophageal-models-as-is.json';
const reportPath=existsSync(currentReport)?currentReport:'reports/ta98-male/approved-oesophageal-layer-installation.json';
const hash=buffer=>createHash('sha256').update(buffer).digest('hex');
const compiled=await build({stdin:{contents:"export * from './app/viewer-state'; export * from './app/bowel-wall-peeling';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {resolveSelection,bowelWallPeelSteps}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
test('installed local oesophageal catalogues deliver valid binary routes, procedural materials and replacement scopes',{skip:!existsSync(reportPath)},()=>{
 const report=JSON.parse(readFileSync(reportPath,'utf8'));
 for(const row of report.catalogues){
  const raw=readFileSync(row.path);const atlas=JSON.parse(raw);const sex=row.sex??'male';
  if(report.models)assert.equal(atlas.oesophagealWallGeometry.candidateSha256,report.models.find(m=>m.sex===sex).candidateSha256);
  else assert.equal(hash(raw),row.newSha256);
  const parts=atlas.parts.filter(p=>!p.suppressed&&p.tissueWallAssembly==='oesophageal:'+sex);assert.equal(parts.length,sex==='male'?14:7);
  const ids=parts.map(p=>p.id),whole=sex==='male'?'ZA:Oesophagus':'LOCAL:female:esophagus:1043';assert.deepEqual(new Set(resolveSelection(atlas,[whole])),new Set(ids));
  if(sex==='male'){const distal=resolveSelection(atlas,['GASTRIC:male:cardiac-connection']);assert.equal(distal.length,7);assert.ok(distal.every(id=>id.includes(':connector:')));}
  assert.equal(bowelWallPeelSteps(atlas,{id:whole,elements:ids}).length,7);
  for(const p of parts){
   const material=atlas.materials[p.material];assert.ok(material?.tissueAppearance);assert.deepEqual(material.tissueAppearance,p.tissueAppearance);
   assert.equal(material.map,undefined);assert.equal(material.normalMap,undefined);
   const chunk=atlas.chunks[p.chunk];assert.ok(chunk.url.startsWith('/local-models/'));
   const path=chunk.url.replace('/local-models/','.local-models/');const bytes=readFileSync(path);assert.equal(bytes.length,chunk.bytes);assert.deepEqual(gunzipSync(readFileSync(path+'.gz')),bytes);
   for(const [offset,count] of [[p.positions,p.vertexCount*12],[p.normals,p.vertexCount*6],[p.indices,p.indexCount*4],[p.tissueChart,p.vertexCount*6],[p.normalSigns,p.vertexCount]])assert.ok(offset>=0&&offset+count<=bytes.length);
  }
 }
 for(const [path,expected] of Object.entries(report.femaleCataloguesPreserved??{}))assert.equal(hash(readFileSync(path)),expected);
});
