import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';
import {PHYSIOLOGY_MAP,buildPhysiologyHierarchy} from './build-physiology-hierarchy.mjs';
const read=path=>JSON.parse(readFileSync(path,'utf8'));
const input=await loadChakraInputs(),hierarchy=input.parseGuestHierarchy(read('public/assets/physiology.json')),evidence=read('public/assets/physiology-provenance.json');
assert.deepEqual(buildPhysiologyHierarchy(input).hierarchy,hierarchy,'Rebuild physiology after changing its mapping');
const ids=new Set(input.catalogues.flatMap(({atlas})=>atlas.parts.filter(part=>!part.suppressed).map(part=>part.id)));
const targets={genes:new Set(read('public/assets/genes.json').nodes.map(node=>node.id)),'cell-types':new Set(read('public/assets/cell-types.json').nodes.map(node=>node.id))};
const credits=readFileSync('ATTRIBUTION.md','utf8');
for(const source of Object.values(evidence.sources))assert(credits.includes(source.url),`Missing reference ${source.url}`);
assert.equal(readFileSync('public/ATTRIBUTION.md','utf8'),credits.replaceAll('(LICENSE)','(LICENSE.txt)').replaceAll('(LICENSES/MIT-upstream.txt)','(MIT-upstream.txt)'));
for(const node of hierarchy.nodes){
 assert(!node.description&&!node.source,'Keep scientific prose and credits out of the viewer');
 for(const id of node.partIds??[])assert(ids.has(id),`Unknown mesh ${id}`);
 for(const link of node.links??[])assert(targets[link.hierarchy]?.has(link.id),`Unknown cross-link ${link.hierarchy}:${link.id}`);
}
const bundle=await build({stdin:{contents:"export * from './app/hierarchical-explosion';export * from './app/guest-choice';export * from './app/hierarchy-navigation';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
const flatten=roots=>{const result=new Map(),visit=node=>{if(result.has(node.id))return;result.set(node.id,node);node.children.forEach(visit);};roots.forEach(visit);return result;};
for(const {file,atlas} of input.catalogues){
 const roots=input.resolveGuestHierarchy(atlas,hierarchy),resolved=flatten(roots),available=atlas.parts.filter(part=>!part.suppressed),availableIds=new Set(available.map(part=>part.id)),names=new Map(atlas.concepts.map(concept=>[concept.id,concept.name]));
 for(const row of PHYSIOLOGY_MAP.functions){
  const node=resolved.get(row.id);assert(node);
  assert.equal(new Set(node.parts.map(part=>part.id)).size,node.parts.length,'Do not duplicate meshes across overlapping function roles');
  for(const part of node.parts)assert(availableIds.has(part.id));
  if(node.parts.length){const choice=api.guestNodeChoice('physiology',node);assert.deepEqual(new Set(choice.elements),new Set(node.parts.map(part=>part.id)));}
 }
 const namesOf=id=>resolved.get('ANATOMY:'+id).parts.map(part=>names.get(part.conceptId)??part.name);
 assert(namesOf('retina').every(name=>/\bretina\b/i.test(name)&&!/retinaculum/i.test(name)),'Retinacula are not retina');
 assert(namesOf('thyroid').every(name=>!/(?:parathyroid|cartilage|ligament)/i.test(name)),'Thyroid cartilage is not the thyroid gland');
 assert(namesOf('systemic-arteries').every(name=>!/pulmonary|\blung\b|umbilical|placent/i.test(name)),'Separate pulmonary and systemic artery territories');
 for(const part of resolved.get('ANATOMY:skin').parts)assert(input.isSkinPart(part));
 if(file.includes('reference')){
  assert.equal(resolved.get('ANATOMY:kidneys').parts.length,2);
  assert.equal(resolved.get('ANATOMY:retina').parts.length,2);
  assert.equal(resolved.get('ANATOMY:parathyroid').parts.length,4);
  assert(resolved.get('ANATOMY:skin').parts.length&&resolved.get('ANATOMY:nodes').parts.length,'Map source names and display-system differences');
  assert(resolved.get('PHYS:quiet-inspiration').parts.some(part=>/diaphragm/i.test(names.get(part.conceptId)??part.name)));
  assert(!resolved.get('PHYS:quiet-inspiration').parts.some(part=>/internal intercostal/i.test(names.get(part.conceptId)??part.name)),'Quiet inspiration must not select forced expiratory muscles');
  assert(resolved.get('PHYS:insulin-secretion').parts.every(part=>/pancrea/i.test(names.get(part.conceptId)??part.name)));
  const choice=api.guestNodeChoice('physiology',resolved.get('PHYS:quiet-inspiration'));
  const ancestors=api.hierarchyAncestors(api.hierarchyNavigation(atlas,'guest:physiology',hierarchy),choice);
  assert(ancestors.some(node=>node.name==='Respiration'));
 }
 const layout=api.createHierarchicalExplosionLayout(atlas,input.hierarchyEntries(atlas),availableIds,'guest',1.7,roots);
 for(const stage of layout.stages){assert.equal(stage.clusters.flat().length,available.length);assert.equal(new Set(stage.clusters.flat()).size,available.length);assert(stage.positions.every(Number.isFinite));}
 console.log(`${file}: ${evidence.coverage[file].mappedFunctions}/${PHYSIOLOGY_MAP.functions.length} functions with geometry; ${available.length} meshes assigned once`);
}
const state={selected:[],visible:[],guestQuery:'PHYS:quiet-inspiration',explode:0,isolate:false,view:'front',rotate:false,reset:0};
const url=new URL(api.viewUrl('http://localhost:3016','local-reference',state));url.searchParams.set('tree','guest:physiology');
assert.equal(api.readViewUrl(url.search,input.catalogues[0].atlas,state).guestQuery,state.guestQuery);
const menu=readFileSync('app/use-guest-hierarchies.ts','utf8');assert(menu.includes("'genes','cell-types','physiology','dermatomes-myotomes'"),'Physiology follows Cell Types');
console.log('Physiology graph, references, cross-links, anatomy, selections, URL restoration and Explode validated.');
