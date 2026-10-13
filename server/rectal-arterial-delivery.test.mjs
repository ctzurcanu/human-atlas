import test from 'node:test';
import assert from 'node:assert/strict';
import {existsSync, readFileSync} from 'node:fs';
import {build} from 'esbuild';

const reportPath='reports/installed-male-rectal-arterial-course.json';
const compiled=await build({stdin:{contents:"export * from './app/viewer-state'; export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
test('superior rectal artery resolves to one fitted source in the ordinary male TA98 tree',{skip:!existsSync(reportPath)},()=>{
 const report=JSON.parse(readFileSync(reportPath));
 assert.equal(report.catalogues.length,2);
 for(const row of report.catalogues){
  const atlas=JSON.parse(readFileSync(row.path));
  for(const id of ['FMA14832','FJ3656','ZA:Superior anorectal artery'])
   assert.deepEqual(api.resolveSelection(atlas,[id]),['ZA:Superior anorectal artery']);
  const tree=api.resolveGuestHierarchy(atlas,api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json'))));
  const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};tree.forEach(walk);
  const node=nodes.get('A12.2.12.073');assert.ok(node);
  assert.deepEqual(node.directParts.map(p=>p.id),['ZA:Superior anorectal artery']);
  assert.equal(atlas.parts.find(p=>p.id==='FJ3656').suppressed,true);
  assert.ok(!api.resolveSelection(atlas,atlas.defaultHidden??[]).includes('ZA:Superior anorectal artery'),'a hidden legacy alternative must not hide its active replacement');
 }
});
const femaleReportPath='reports/installed-female-rectal-arterial-course.json';
test('the best-present female source is delivered once under the superior rectal artery TA98 control',{skip:!existsSync(femaleReportPath)},()=>{
 const report=JSON.parse(readFileSync(femaleReportPath));assert.equal(report.catalogues.length,3);
 for(const row of report.catalogues){
  const atlas=JSON.parse(readFileSync(row.path));
  assert.deepEqual(api.resolveSelection(atlas,['LOCAL:female:superior_rectal_artery']),[report.partId]);
  assert.equal(atlas.parts.filter(p=>p.id===report.partId&&!p.suppressed).length,1);
  const tree=api.resolveGuestHierarchy(atlas,api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json'))));
  const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};tree.forEach(walk);
  assert.deepEqual(nodes.get('A12.2.12.073').directParts.map(p=>p.id),[report.partId]);
  assert.ok(!api.resolveSelection(atlas,atlas.defaultHidden??[]).includes(report.partId));
 }
});
