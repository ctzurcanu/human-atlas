import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const b=await build({stdin:{contents:"export * from './app/biliary-junction';export * from './app/digestive-duct-relations';export * from './app/gallbladder-wall';export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(b.outputFiles[0].contents).toString('base64'));
const read=p=>JSON.parse(fs.readFileSync(p));
const a=read('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json');
const cystic=a.parts.find(p=>p.id==='LOCAL:female:cystic_duct:authored'),trunk=a.parts.find(p=>p.id==='LOCAL:female:bile_duct:1031');
test('female native cystic wall has reciprocal open trunk and four Ga tissue joins',()=>{
 assert.ok(api.biliaryJunctionInterface(cystic,trunk));
 const graph=api.digestiveDuctRelations(a.parts);
 assert.ok(graph.get(cystic.id).some(r=>r.kind==='continuous'&&r.target.id===trunk.id));
 assert.ok(graph.get(trunk.id).some(r=>r.kind==='continuous'&&r.target.id===cystic.id));
 const coats=a.gallbladderWallPreview.partIds.map(id=>a.parts.find(p=>p.id===id));
 const neckCoats=coats.filter(p=>p.tissueDrawGroups.some(g=>g.neighbor===cystic.id));assert.equal(neckCoats.length,4);
 for(const coat of coats){assert.equal(api.gallbladderWallInterface(coat,cystic),neckCoats.includes(coat));assert.equal(api.gallbladderWallInterface(coat,trunk),false);}
 for(const change of [{suppressed:true},{sourceOffset:[0,0,0]},{chunk:cystic.chunk+1},{positions:4},{sourceId:'foreign'},{conceptId:'foreign'},{pairedTissue:false},{sectionSolid:false},{tissueDrawGroups:[]}])assert.equal(api.biliaryJunctionInterface({...cystic,...change},trunk),false);
 const proof=read(`reports/ta98-female/${a.femaleGallbladderRegions?'gallbladder-regions':a.cysticFoldIntegration?'cystic-folds':'cystic-local-connection'}-packed-validation.json`);
 assert.equal(proof.literalReciprocalContacts.length,a.femaleGallbladderRegions?30:a.cysticFoldIntegration?9:8);assert.ok(proof.continuousMacroLumen);
 assert.ok(proof.macroBoundaries.every(r=>r.components===1&&r.openTerminalLoops===1&&r.euler===1&&r.selfIntersections===0));
 assert.equal(proof.catalogueSha256,read(a.femaleGallbladderRegions?'.local-models/ta98-review/female-gallbladder-regions/construction.json':a.cysticFoldIntegration?'.local-models/ta98-review/cystic-folds/female/construction.json':'reports/ta98-female/cystic-local-connection-installed.json').catalogueSha256);
});
test('TA98 .011 selects the distinct native cell and retains explicit finer-feature limits',()=>{
 const c=a.concepts.find(c=>c.id==='LOCAL:female:cystic_duct');assert.deepEqual(c.elements,a.cysticFoldIntegration?[cystic.id,a.cysticFoldIntegration.foldPartId]:[cystic.id]);assert.equal(c.ta98Term,'A05.8.02.011');
 const tree=api.resolveGuestHierarchy(a,api.parseGuestHierarchy(read('public/assets/ta98.json')));let node;
 const walk=n=>{if(n.id==='A05.8.02.011')node=n;n.children.forEach(walk);};tree.forEach(walk);
 assert.deepEqual(new Set(node.directParts.map(p=>p.id)),new Set(c.elements));assert.match(api.digestiveDuctDescription(cystic).text,a.cysticFoldIntegration?/Three raised crescentic/:/spiral folds.*under review/);
 const retired=a.parts.find(p=>p.id==='HRA:VH_F_cystic_duct');assert.ok(retired.suppressed&&a.defaultHidden.includes(retired.id));
 assert.equal(c.partOf,undefined);
 const packing=read('reports/ta98-female/cystic-terminal-neck-lowered-packing.json');assert.equal(packing.newPrimaryInnerContacts,0);assert.equal(packing.newPrimaryOuterContacts,0);
 assert.ok(read('reports/ta98-female/cystic-local-connection-installed.json').finiteCells.every(c=>c.components===1&&c.openEdges===0&&c.selfIntersections===0));
});
