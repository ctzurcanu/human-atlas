import assert from 'node:assert/strict';import test from 'node:test';import {readFile,writeFile,mkdir} from 'node:fs/promises';import {build} from 'esbuild';
const built=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].text).toString('base64'));
const described=await build({entryPoints:['app/structure-description.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {localDescription}=await import('data:text/javascript;base64,'+Buffer.from(described.outputFiles[0].text).toString('base64'));
const report={source:'https://www.nhlbi.nih.gov/health/heart/blood-flow',geometryChanged:false,geometryAcceptance:'Not established by relationship metadata',models:[]};
for(const model of ['ta98-male','ta98-female-runtime'])test(`${model}: native valve context and reciprocal links`,async()=>{
 const atlas=JSON.parse(await readFile(`.local-models/${model}.json`)),parts=atlas.parts;
 const valves=parts.filter(p=>p.system==='cardiac'&&(/valve|coronary leaflet/.test(p.name.toLowerCase()))&&(p.id.startsWith('ZA:')||p.id.startsWith('LOCAL:')));
 assert.ok(valves.length>=4);
 const rows=[];
 for(const source of valves){const description=localDescription(source.name,source,[]);assert.ok(description.text.includes('backward blood flow'));assert.ok(description.url.startsWith('https://www.nhlbi.nih.gov/'));if(/\b(?:leaflet|cusp)\b/i.test(source.name))assert.ok(description.text.startsWith('This is a leaflet'));}
 for(const source of valves){
  const links=anatomicalRelations(source,parts).filter(r=>r.kind==='adjacent'&&/valve region separates/.test(r.note??''));assert.ok(links.length,source.name);
  for(const r of links){assert.ok(!r.target.id.startsWith('HRA:'));assert.ok(anatomicalRelations(r.target,parts).some(back=>back.kind==='adjacent'&&back.target.id===source.id),'Reciprocal native adjacency');assert.ok(r.note.includes('remain unverified'));if(r.target.name==='Ventricles')assert.ok(r.note.includes('separate left/right ventricular representation is missing'));}
  rows.push({id:source.id,name:source.name,links:links.map(r=>({id:r.target.id,name:r.target.name,note:r.note}))});
 }
 report.models.push({model,reviewedValveRecords:valves.length,relations:rows});
 await mkdir('reports/ta98-completeness',{recursive:true});await writeFile('reports/ta98-completeness/cardiac-relations.json',JSON.stringify(report,null,2)+'\n');
});
test('No cross-frame, suppressed, unrelated-system or unsided chamber substitutions',()=>{
 const p=(id,name,system='cardiac')=>({id,conceptId:id,name,system});
 const valve=p('LOCAL:female:mitral','Mitral valve');
 const good=p('LOCAL:female:atrium','Atrium (left)'),right=p('LOCAL:female:right','Atrium (right)'),brain=p('LOCAL:female:brain','Left ventricle','nervous'),foreign=p('LOCAL:male:atrium','Left atrium'),donor=p('HRA:VH_F_left_ventricle','Heart left ventricle'),hidden={...p('LOCAL:female:hidden','Left ventricle'),suppressed:true};
 const parts=[valve,good,right,brain,foreign,donor,hidden];const links=anatomicalRelations(valve,parts).filter(r=>r.kind==='adjacent');assert.deepEqual(links.map(r=>r.target.id),[good.id]);
 const semilunar=p('LOCAL:female:pulmonary','Pulmonary valve'),papillary=p('LOCAL:female:papillary','Anterior papillary muscle of right ventricle'),rv=p('LOCAL:female:rv','Right ventricle');const pp=[semilunar,papillary,rv];
 assert.ok(anatomicalRelations(papillary,pp).some(r=>r.kind==='adjacent'&&r.target.id===rv.id));assert.ok(!anatomicalRelations(papillary,pp).some(r=>r.target.id===semilunar.id),'No chordal/papillary relation to semilunar valves');
});

test('Bicuspid is a whole valve label, not a cusp selection',()=>{const part={id:'LOCAL:female:bicuspid',name:'Bicuspid atrioventricular valve (left)',system:'cardiac'};assert.ok(localDescription(part.name,part,[]).text.startsWith('The mitral valve controls flow'));});
