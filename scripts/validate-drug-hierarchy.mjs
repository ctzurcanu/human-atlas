import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';

const json=path=>JSON.parse(readFileSync(path,'utf8'));
const raw=json('public/assets/drugs.json'),provenance=json('public/assets/drug-provenance.json');
const input=await loadChakraInputs(),hierarchy=input.parseGuestHierarchy(raw),byId=new Map(hierarchy.nodes.map(node=>[node.id,node]));
const genes=new Map(json('public/assets/genes.json').nodes.map(node=>[node.id,node]));
const cells=new Map(json('public/assets/cell-types.json').nodes.map(node=>[node.id,node]));
assert.deepEqual(hierarchy.roots,[...Array.from('ABCDGHJLMNPRSV',code=>'ATC:'+code),'UNCLASSIFIED']);
for(const row of json('work/drugs/atc.json').records){
 for(let level=1;level<=5;level++){
  const node=byId.get('ATC:'+row['level'+level]);assert(node,`Missing ATC branch ${row['level'+level]}`);
  if(level>1)assert(byId.get('ATC:'+row['level'+(level-1)]).childIds.includes(node.id));
 }
}
const molecules=new Map([...json('work/drugs/approved.json').records,...json('work/drugs/classified.json').records].map(row=>[row.molecule_chembl_id,row]));
for(const [id,molecule] of molecules){
 assert(byId.has(id),`Missing source drug ${id}`);
 if(molecule.atc_classifications.length)for(const code of molecule.atc_classifications)assert(byId.get('ATC:'+code).childIds.includes(id));
 else assert(byId.get('UNCLASSIFIED').childIds.includes(id));
}
assert.equal(molecules.size,provenance.statistics.drugs);
for(const node of hierarchy.nodes){
 assert(!node.source&&!node.description,'Keep credits and research notes out of the displayed hierarchy');
 if(node.id.startsWith('HGNC:'))assert.deepEqual(node.matches??[],genes.get(node.id).matches??[],`Unsubstantiated gene anatomy ${node.id}`);
 for(const link of node.links??[])assert((link.hierarchy==='genes'?genes:cells).has(link.id),`Broken hierarchy link ${link.id}`);
}
for(const row of provenance.mechanisms){
 const node=byId.get('MECHANISM:'+row.mec_id);assert(node);assert(byId.get(row.molecule_chembl_id).childIds.includes(node.id));
 if(row.humanGenes.length){assert.equal(row.target.tax_id,9606);assert.equal(row.target.organism,'Homo sapiens');assert.equal(row.molecular_mechanism,1);}
 else assert(!node.childIds?.some(id=>id.startsWith('HGNC:')),'Do not map a nonhuman target to human homologues');
 assert.deepEqual([...(node.childIds??[])].sort(),[...row.humanGenes].sort());
}
const bundle=await build({stdin:{contents:"export * from './app/hierarchical-explosion';export * from './app/hierarchy-navigation';export * from './app/viewer-state';export * from './app/guest-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
const flatten=roots=>{const byId=new Map(),visit=node=>{if(byId.has(node.id))return;byId.set(node.id,node);node.children.forEach(visit);};roots.forEach(visit);return byId;};
let mapped=false;
for(const {file,atlas} of input.catalogues){
 const start=performance.now(),roots=input.resolveGuestHierarchy(atlas,hierarchy),resolved=flatten(roots);
 const available=atlas.parts.filter(part=>!part.suppressed),ids=new Set(available.map(part=>part.id));
 for(const node of resolved.values()){
  assert.equal(node.parts.length,new Set(node.parts.map(part=>part.id)).size,`Duplicate geometry: ${node.id}`);
  for(const part of node.parts)assert(ids.has(part.id));
  if(node.id.startsWith('CHEMBL')&&node.parts.length){mapped=true;const choice=api.guestNodeChoice('drugs',node);assert(choice.children.length);assert.deepEqual(new Set(choice.elements),new Set(node.parts.map(part=>part.id)));}
 }
 const nav=api.hierarchyNavigation(atlas,'guest:drugs',hierarchy);assert(nav.children.length);
 const entries=input.hierarchyEntries(atlas),layout=api.createHierarchicalExplosionLayout(atlas,entries,ids,'guest',1.7,roots);
 for(const stage of layout.stages){assert.equal(stage.clusters.flat().length,available.length);assert.equal(new Set(stage.clusters.flat()).size,available.length);assert(stage.positions.every(Number.isFinite));}
 console.log(`${file}: ${[...resolved.values()].filter(node=>node.id.startsWith('CHEMBL')&&node.parts.length).length} drugs with anatomy; ${available.length} meshes assigned once; ${Math.round(performance.now()-start)} ms`);
}
assert(mapped,'At least one curated human target must connect to modeled anatomy');
const state={selected:[],visible:[],guestQuery:'ATC:A10BA02',explode:0,isolate:false,view:'front',rotate:false,reset:0};
const url=new URL(api.viewUrl('http://localhost:3016','local-reference',state));url.searchParams.set('tree','guest:drugs');
assert.equal(api.readViewUrl(url.search,input.catalogues[0].atlas,state).guestQuery,state.guestQuery);
console.log(`${provenance.statistics.atcSourceSubstances} ATC substance branches, ${molecules.size} source drugs and ${provenance.statistics.mechanisms} mechanisms validated.`);
