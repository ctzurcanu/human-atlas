import assert from 'node:assert/strict';
import test from 'node:test';
import {existsSync,readFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const built=await build({stdin:{contents:`
 export * from './app/lingual-arterial-relations';
 export {hierarchyEntries} from './app/anatomy-hierarchy';
 export {buildSystemTree} from './app/system-tree';
 export {nativePartFrame,ta98ModeledParent} from './app/ta98-modeled-parent';
 `,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));

test('Metadata-only unit fixtures retain old layouts and reject incomplete or mixed common-bank markers',()=>{
 const roles={lingual:'A12.2.05.015',suprahyoid:'A12.2.05.016','dorsal-lingual-0':'A12.2.05.017','dorsal-lingual-1':'A12.2.05.017',sublingual:'A12.2.05.018','deep-lingual':'A12.2.05.019'};
 const origins={left:[.016594554721975487,1.623810152874436,.04152941869862918],right:[-.01985116613921688,1.6236474710160156,.043219807199565606]};
 const marker={carotidFiniteWallRevision:'left-four-cell-common-bank-v2',carotidFiniteWallBankSha256:'e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36'};
 const family=(side,common=false)=>{
  const parent=`ZA:External carotid artery.${side[0]}`,ids=Object.keys(roles).map(role=>`TA98ARTERIAL:male:${role}:${side}:finite-wall`),main=ids[0];
  const defaults={system:'arterial',pairedTissue:true,sectionSolid:true,normalSigns:600,positions:0,normals:100,uvs:200,vertexCount:1000,chunk:side==='left'?0:1,
   geometryOrigin:common?[0,1.625,0]:origins[side],sectionAssembly:`${side}-lingual-arterial-family:male`,lingualArterialWallRevision:`male-${side}-lingual-finite-network-v${side==='left'?5:2}`,
   provenance:{url:'https://github.com/Z-Anatomy/Models-of-human-anatomy'},...(common?marker:{})};
  return [{...defaults,id:parent,conceptId:parent,sourceId:`External carotid artery.${side[0]}`,indexCount:72,tissueDrawGroups:[{start:0,count:72,neighbor:main}]},
   ...Object.entries(roles).map(([role,term],i)=>({...defaults,id:ids[i],conceptId:role.startsWith('dorsal-lingual-')?`TA98ARTERIAL:male:dorsal-lingual:${side}:whole`:ids[i],sourceId:'independent-lingual-course/native-ECA-and-recipient-constrained',ta98Term:term,
    indexCount:role==='lingual'?432:72,tissueDrawGroups:(role==='lingual'?[parent,...ids.slice(1)]:[main]).map((neighbor,i)=>({start:i*72,count:72,neighbor}))}))];
 };
 // These are guard unit fixtures, not packed geometry or anatomical evidence.
 const original=[...family('left'),...family('right')],common=[...family('left',true),...family('right')];
 for(const parts of [original,common])for(const side of ['left','right'])assert.equal(api.registeredLingualFamily(parts,side)?.size,6);
 const id='TA98ARTERIAL:male:lingual:left:finite-wall';
 for(const change of [{carotidFiniteWallRevision:undefined},{carotidFiniteWallBankSha256:undefined},{carotidFiniteWallRevision:null},{carotidFiniteWallBankSha256:''},
  {carotidFiniteWallRevision:'unreviewed'},{carotidFiniteWallBankSha256:'bad'},{geometryOrigin:origins.left},{uvs:undefined},{positions:4},{tissueDrawGroups:[]}]){
  const parts=common.map(p=>p.id===id?{...p,...change}:p);
  assert.equal(api.registeredLingualFamily(parts,'left'),undefined);
  assert.equal(api.registeredLingualFamily(parts,'right')?.size,6);
 }
 assert.equal(api.registeredLingualFamily(original.map(p=>p.id===id?{...p,...marker}:p),'left'),undefined);
 assert.equal(api.registeredLingualFamily(common.map(p=>p.id===id?original.find(q=>q.id===id):p),'left'),undefined);
 assert.equal(api.registeredLingualFamily([...family('left'),...family('right',true)],'right'),undefined);
 const bifId='TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall';
 const bif={...common[1],id:bifId,conceptId:bifId,ta98Term:'A12.2.04.009',name:'Carotid bifurcation (left)',
  sourceId:'Constructed left carotid bifurcation tissue region',license:'CC BY-SA 4.0',sectionAssembly:'left-carotid-lingual-coupled-finite-wall:male'};
 assert.equal(api.isInstalledCarotidBifurcationPart(bif),true);assert.equal(api.nativePartFrame(bif),'ZA:');
 for(const change of [{id:bifId.replace(':left:',':right:')},{conceptId:'other'},{ta98Term:'A12.2.04.008'},{sourceId:'other'},
  {license:'unverified'},{provenance:{url:'https://example.com'}},{pairedTissue:false},{sectionSolid:false},{normalSigns:undefined},
  {sectionAssembly:'left-lingual-arterial-family:male'},{carotidFiniteWallBankSha256:'bad'},{geometryOrigin:origins.left},{sourceOffset:[0,0,0]}]){
  const bad={...bif,...change};assert.equal(api.isInstalledCarotidBifurcationPart(bad),false);assert.equal(api.nativePartFrame(bad),undefined);
 }
});

const base='.local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-completion-v1';
const proposedFile=[`${base}/catalogue-adapter-v1/proposed-v3.json`,`${base}/catalogue-adapter-v1/proposed-v2.json`,`${base}/catalogue-adapter-v1/proposed.json`].find(existsSync);
const oldBankFile='.local-models/geometry-workbench/circulatory-completion/male-left-lingual-finite-network-v5/owned-bank-v2/tissue.json';
const masterFile='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
if(!proposedFile||!existsSync(oldBankFile)||!existsSync(masterFile)){
 test('Audited carotid/lingual bank checks require private local anatomy',{skip:'Local geometry is intentionally absent'},()=>{});
}else{
const bankSha='e4393d63924d8bcada2655edfe520dca71a66dc9a0709a138649b37213ba7c36';
const revision='left-four-cell-common-bank-v2';
const roles=['lingual','suprahyoid','dorsal-lingual-0','dorsal-lingual-1','sublingual','deep-lingual'];
const leftIds=['ZA:External carotid artery.l',...roles.map(role=>`TA98ARTERIAL:male:${role}:left:finite-wall`)];
const leftSet=new Set(leftIds);
const sourceBytes=readFileSync(proposedFile),masterBytes=readFileSync(masterFile);
const proposed=JSON.parse(sourceBytes),old=JSON.parse(readFileSync(oldBankFile));
// The first actual ten-cell proposal predates the adapter's metadata fix.
// Model that requested fix in memory; retain its actual indices/bank buffers.
const fixture=/\/proposed-v[23]\.json$/.test(proposedFile)?proposed:{...proposed,parts:proposed.parts.map(p=>leftSet.has(p.id)?{...p,
 sectionAssembly:'left-lingual-arterial-family:male',carotidFiniteWallRevision:revision,
 carotidFiniteWallBankSha256:bankSha}:p)};
const legacy={...fixture,parts:fixture.parts.map(p=>{
 if(!leftSet.has(p.id))return p;
 const row=Object.values(old.parts).find(row=>row.id===p.id);
 assert.ok(row,`Original seven-cell record ${p.id}`);
 const restored={...p,...row,sectionAssembly:'left-lingual-arterial-family:male'};
 delete restored.carotidFiniteWallRevision;delete restored.carotidFiniteWallBankSha256;
 return restored;
})};
const graphSignature=parts=>[...api.lingualArterialRelations(parts)].flatMap(([source,rows])=>
 rows.map(row=>[source,row.kind,row.target.id,row.targetConcept?.id??'',row.note])).sort((a,b)=>JSON.stringify(a).localeCompare(JSON.stringify(b)));
const replace=(parts,id,change)=>parts.map(p=>p.id===id?{...p,...change}:p);

test('Actual proposed ten-cell bank retains both lingual systems and all reciprocal recipient links',()=>{
 const parent=fixture.parts.find(p=>p.id===leftIds[0]);
 const bankPath=fixture.chunks[parent.chunk].url.replace(/^\/local-models\//,'.local-models/');
 const bank=readFileSync(bankPath);
 assert.equal(createHash('sha256').update(bank).digest('hex'),bankSha);
 const ten=fixture.parts.filter(p=>p.chunk===parent.chunk);
 assert.equal(ten.length,10,'Use the actual four carotid plus six coupled lingual cell records');
 for(const side of ['left','right']){
  assert.equal(api.registeredLingualFamily(legacy.parts,side)?.size,6,'Original banks remain valid');
  assert.equal(api.registeredLingualFamily(fixture.parts,side)?.size,6,'Both sides remain registered');
 }
 for(const p of fixture.parts.filter(p=>leftSet.has(p.id)&&p.id!==leftIds[0]))assert.equal(api.nativePartFrame(p),'ZA:');
 assert.deepEqual(graphSignature(fixture.parts),graphSignature(legacy.parts),'Registered branch, tongue, muscle and gland graph is unchanged');
 const systems=api.buildSystemTree(api.hierarchyEntries(fixture),'systems',fixture);
 const hits=[];const visit=nodes=>{for(const n of nodes){if(n.name==='Lingual artery')hits.push(n);if(n.kind==='group')visit(n.nodes);}};
 for(const system of systems)visit(system.nodes);
 assert.equal(hits.length,1);assert.equal(hits[0].parts.length,12);assert.equal(hits[0].nodes.length,5);
 // Check the actual retained reciprocal 72-index annuli, including ECA's
 // additional bifurcation interface, without requiring equal sign buffers.
 const family=api.registeredLingualFamily(fixture.parts,'left'),main=family.get('lingual');
 const rim=(p,neighbor)=>{
  const g=p.tissueDrawGroups.find(g=>g.neighbor===neighbor);assert.equal(g.count,72);
  return Array.from({length:24},(_,i)=>Array.from({length:3},(_,k)=>bank.readUInt32LE(p.indices+(g.start+i*3+k)*4)));
 };
 for(const other of [parent,...[...family.values()].filter(p=>p!==main)]){
  const first=rim(main,other.id),second=rim(other,main.id);
  assert.deepEqual(first,second.map(([a,b,c])=>[a,c,b]));
  for(const tri of first)for(const i of tri)assert.equal(bank.readInt8(main.normalSigns+i),-bank.readInt8(other.normalSigns+i));
  assert.equal(main.positions,other.positions);assert.equal(main.normals,other.normals);assert.equal(main.uvs,other.uvs);
 }
 assert.deepEqual(readFileSync(proposedFile),sourceBytes);assert.deepEqual(readFileSync(masterFile),masterBytes);
});

test('Common origin requires complete exact left bank markers, one frame and intact shared buffers/rims',()=>{
 const id=roles.map(role=>`TA98ARTERIAL:male:${role}:left:finite-wall`)[0];
 const target=fixture.parts.find(p=>p.id===id);
 const changes=[
  {carotidFiniteWallRevision:undefined},{carotidFiniteWallBankSha256:undefined},
  {carotidFiniteWallRevision:null},{carotidFiniteWallBankSha256:''},
  {carotidFiniteWallRevision:'unreviewed'},{carotidFiniteWallBankSha256:'0'.repeat(64)},
  {geometryOrigin:[0,1.625,.000001]},{geometryOrigin:old.parts.lingual.geometryOrigin},
  {sectionAssembly:'left-carotid-lingual-coupled-finite-wall:male'},
  {positions:target.positions+4},{normals:target.normals+4},{uvs:undefined},
  {vertexCount:target.vertexCount-1},{chunk:target.chunk+1},
  {sourceOffset:[0,0,0]},{suppressed:true},
  {tissueDrawGroups:target.tissueDrawGroups.filter(g=>g.neighbor!==leftIds[0])},
  {tissueDrawGroups:target.tissueDrawGroups.concat(target.tissueDrawGroups.filter(g=>g.neighbor===leftIds[0]))},
 ];
 for(const change of changes){
  const parts=replace(fixture.parts,id,change);
  assert.equal(api.registeredLingualFamily(parts,'left'),undefined,JSON.stringify(change));
  assert.equal(api.registeredLingualFamily(parts,'right')?.size,6);
  assert.ok(!api.lingualArterialRelations(parts).get(id)?.some(r=>['before','after','supplies'].includes(r.kind)));
 }
 for(const id of leftIds){
  const p=fixture.parts.find(p=>p.id===id);
  assert.equal(api.registeredLingualFamily(replace(fixture.parts,id,{carotidFiniteWallBankSha256:'bad'}),'left'),undefined,p.id);
 }
 // Neither an old-origin fallback with partial markers nor a wholly copied
 // left marker on the original right family establishes a reviewed frame.
 const oldWithMarker=replace(legacy.parts,id,{carotidFiniteWallRevision:revision});
 assert.equal(api.registeredLingualFamily(oldWithMarker,'left'),undefined);
 const mixed=fixture.parts.map(p=>p.id===id?legacy.parts.find(q=>q.id===id):p);
 assert.equal(api.registeredLingualFamily(mixed,'left'),undefined);
 const rightIds=new Set(['ZA:External carotid artery.r',...roles.map(role=>`TA98ARTERIAL:male:${role}:right:finite-wall`)]);
 const copied=fixture.parts.map(p=>rightIds.has(p.id)?{...p,geometryOrigin:[0,1.625,0],carotidFiniteWallRevision:revision,carotidFiniteWallBankSha256:bankSha}:p);
 assert.equal(api.registeredLingualFamily(copied,'right'),undefined);
});

test('Actual bifurcation has only its reviewed native frame and ipsilateral common-carotid parent',()=>{
 const id='TA98ARTERIAL:male:carotid-bifurcation:left:finite-wall';
 const part=fixture.parts.find(p=>p.id===id),concept=fixture.concepts.find(c=>c.id===id);
 assert.ok(part&&concept);assert.equal(api.isInstalledCarotidBifurcationPart(part),true);
 assert.equal(api.nativePartFrame(part),'ZA:');assert.equal(concept.partOf,'ZA:Left common carotid artery');
 const parts=new Map(fixture.parts.map(p=>[p.id,p]));
 assert.equal(api.ta98ModeledParent(concept,fixture.concepts,parts)?.id,'ZA:Left common carotid artery');
 // Exercise the frame guard even when the latest whole-CCA concept also
 // explicitly owns this region: remove only that membership in memory.
 const direct=fixture.concepts.map(c=>c.id==='ZA:Left common carotid artery'?{...c,elements:c.elements.filter(p=>p!==id)}:c);
 assert.equal(api.ta98ModeledParent(concept,direct,parts)?.id,'ZA:Left common carotid artery');
 const noLeft=direct.filter(c=>c.id!=='ZA:Left common carotid artery');
 assert.equal(api.ta98ModeledParent(concept,noLeft,parts),undefined,'Opposite common carotid cannot replace its left parent');
 for(const change of [{carotidFiniteWallBankSha256:'bad'},{sourceId:'other'},{license:'unverified'},
  {sectionAssembly:'left-lingual-arterial-family:male'},{geometryOrigin:[0,0,0]}]){
  const bad={...part,...change},changed=new Map(parts);changed.set(id,bad);
  assert.equal(api.nativePartFrame(bad),undefined);
  assert.equal(api.ta98ModeledParent(concept,direct,changed),undefined);
 }
});
}
