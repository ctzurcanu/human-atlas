// Metadata and actual frontend control availability only; never copies geometry.
import assert from 'node:assert/strict';import fs from 'node:fs';import path from 'node:path';import{createHash}from'node:crypto';import{build}from'esbuild';
const root=process.cwd(),sha=b=>createHash('sha256').update(b).digest('hex'),work=path.join(root,'.local-models/geometry-workbench/circulation-census-v1');
const inventory=JSON.parse(fs.readFileSync('reports/ta98-circulation-inventory.json'));
const missing=inventory.terms.filter(r=>r.male_applicable&&!r.resolved_direct_part_ids.length);assert.equal(missing.length,815);
const codes=new Set(missing.map(r=>r.id));
const compiled=await build({stdin:{contents:`export * from './app/guest-hierarchy';export * from './app/reviewed-catalogue-identities';`,resolveDir:root},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const definition=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
const sources=[['female-current','.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json','model:local-female-ta98:assembled-metres'],['male-original','.local-models/male.json','model:local-male:assembled-metres'],['female-original','.local-models/female.json','model:local-female:assembled-metres'],['reference-original','.local-models/reference.json','model:local-reference:assembled-metres'],['male-preinstalled-ta98','.local-models/ta98-male.json','model:local-ta98:assembled-metres'],['female-preinstalled-ta98','.local-models/ta98-female.json','model:local-female-ta98:assembled-metres']];
const out=[];
for(const [family,file,frame]of sources){
 const raw=fs.readFileSync(file),a=JSON.parse(raw),before=sha(JSON.stringify(a)),parts=new Map(a.parts.map(p=>[p.id,p])),hidden=new Set(a.defaultHidden??[]);
 const roots=api.resolveGuestHierarchy(api.withReviewedCatalogueIdentities(a),definition),seen=new Set(),hits=[];
 const visit=n=>{if(seen.has(n.id))return;seen.add(n.id);if(codes.has(n.id)&&n.directParts.length)hits.push({term:n.id,name:n.name,partIds:n.directParts.map(p=>p.id),sources:n.directParts.map(p=>({id:p.id,name:p.name,conceptId:p.conceptId,system:p.system,suppressed:!!p.suppressed,defaultHidden:hidden.has(p.id),bounds:p.bounds,chunkUrl:a.chunks[p.chunk].url,sourceId:p.sourceId,license:p.license??null,provenance:p.provenance??null,sourceFrameRegistrationEstablished:!p.id.startsWith('HRA:')}))});n.children.forEach(visit);};roots.forEach(visit);
 assert.equal(sha(JSON.stringify(a)),before);assert.equal(sha(fs.readFileSync(file)),sha(raw));
 out.push({family,path:file,sourceSha256:sha(raw),declaredFrame:frame,frameIsOnlyThisCatalogueRenderFrame:true,sourceDescription:a.source,sex:a.sex,hitCount:hits.length,hits});
}
const file=path.join(work,'missing-male-source-catalogue-crosscheck.json');assert.ok(!fs.existsSync(file));fs.writeFileSync(file,JSON.stringify({schema:'circulation-missing-source-catalogue-crosscheck/v1',maleMissingTerms:815,sourceInventorySha256:sha(fs.readFileSync('reports/ta98-circulation-inventory.json')),sources:out,noGeometryDecodedOrCopied:true,noTransferOrAnatomicalAcceptance:true},null,2)+'\n',{flag:'wx'});
console.log(JSON.stringify({output:file,sources:out.map(r=>({family:r.family,hits:r.hitCount})),noGeometryDecodedOrCopied:true}));
