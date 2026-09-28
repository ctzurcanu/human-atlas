import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';
import {EXERCISE_ROOTS,loadExerciseSources} from './build-exercise-hierarchy.mjs';
const json=path=>JSON.parse(readFileSync(path,'utf8'));
const input=await loadChakraInputs(),hierarchy=input.parseGuestHierarchy(json('public/assets/physical-exercise.json'));
const provenance=json('public/assets/exercise-provenance.json'),source=loadExerciseSources(),nodes=new Map(hierarchy.nodes.map(n=>[n.id,n]));
assert.deepEqual(hierarchy.roots,EXERCISE_ROOTS.map(([id])=>'TYPE:'+id));
const spinal=new Set(json('public/assets/dermatomes-myotomes.json').nodes.map(n=>n.id));
for(const row of source.data.exerciseinfo)assert(nodes.has('WGER:'+row.id),`Missing source exercise ${row.id}`);
assert.equal(provenance.exercises.length,provenance.statistics.sourceExercises+provenance.statistics.additionalExercises);
assert.equal(new Set(provenance.exercises.map(e=>e.id)).size,provenance.exercises.length);
const fullCredits=readFileSync('ATTRIBUTION.md','utf8');
for(const reference of Object.values(provenance.sources))assert(fullCredits.includes(reference.url),`Missing credit: ${reference.name}`);
for(const row of provenance.exercises.filter(e=>e.credit)){assert(fullCredits.includes(row.credit.source));assert(new URL(row.credit.license.url).hostname==='creativecommons.org');for(const author of row.credit.authors)assert(fullCredits.includes(author.replace(/[|\r\n]/g,' ').replace(/[[\]]/g,'')));}
assert.equal(readFileSync('public/ATTRIBUTION.md','utf8'),fullCredits.replaceAll('(LICENSE)','(LICENSE.txt)').replaceAll('(LICENSES/MIT-upstream.txt)','(MIT-upstream.txt)'));
for(const node of hierarchy.nodes){assert(!node.source&&!node.description,'Keep source credits and research prose out of the UI');for(const link of node.links??[])assert(link.hierarchy==='dermatomes-myotomes'&&spinal.has(link.id));}
const bundle=await build({stdin:{contents:"export * from './app/hierarchical-explosion';export * from './app/guest-choice';export * from './app/hierarchy-navigation';export * from './app/viewer-state';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
const flatten=roots=>{const result=new Map(),visit=node=>{if(result.has(node.id))return;result.set(node.id,node);node.children.forEach(visit);};roots.forEach(visit);return result;};
const availableIds=new Set(input.catalogues.flatMap(({atlas})=>atlas.parts.filter(p=>!p.suppressed).map(p=>p.id)));
for(const node of hierarchy.nodes)for(const id of node.partIds??[])assert(availableIds.has(id),`Unknown anatomy selector ${id}`);
for(const {file,atlas} of input.catalogues){
 const roots=input.resolveGuestHierarchy(atlas,hierarchy),resolved=flatten(roots),available=atlas.parts.filter(p=>!p.suppressed),ids=new Set(available.map(p=>p.id)),names=new Map(atlas.concepts.map(c=>[c.id,c.name]));
 for(const node of resolved.values()){assert.equal(node.parts.length,new Set(node.parts.map(p=>p.id)).size);for(const part of node.parts)assert(ids.has(part.id));}
 for(const row of provenance.exercises){
  const item=resolved.get(row.id);assert(item);if(item.parts.length){const choice=api.guestNodeChoice(hierarchy.id,item);assert.deepEqual(new Set(choice.elements),new Set(item.parts.map(p=>p.id)));assert(choice.children.length);}
  if(row.side)for(const p of item.parts){const name=names.get(p.conceptId)??p.name;assert(!new RegExp(`\\b${row.side==='left'?'right':'left'}\\b`,'i').test(name),`Crossed laterality in ${row.name}: ${name}`);}
  for(const role of item.children.filter(n=>/:primary$|:secondary$/.test(n.id)))for(const p of role.parts)assert(!/\btendon\b|\.[eo]\d*[lr]$/i.test(names.get(p.conceptId)??p.name),`Tendon/attachment marker in muscle role ${row.name}`);
 }
 if(file.includes('reference')){
  const muscles=id=>resolved.get(id).parts.map(p=>names.get(p.conceptId)??p.name).join(' ');
  assert(/vastus/i.test(muscles('ANATOMY:quads'))&&/rectus femoris/i.test(muscles('ANATOMY:quads')));
  assert(/adductor/i.test(muscles('ROLE:WGER:12:primary'))&&!/gluteus/i.test(muscles('ROLE:WGER:12:primary')));
  assert(/infraspinatus/i.test(muscles('ROLE:WGER:142:primary'))&&!/clavicular.*deltoid/i.test(muscles('ROLE:WGER:142:primary')));
  assert(/acromial.*deltoid/i.test(muscles('ROLE:WGER:348:primary'))&&!/clavicular.*deltoid/i.test(muscles('ROLE:WGER:348:primary')));
  assert(resolved.get('ANATOMY:heart').parts.length&&resolved.get('ANATOMY:lungs').parts.length,'Map subdivided organs, not only whole-organ names');
  assert(resolved.get('EX:pelvic-supine').parts.some(p=>p.system==='urinary'));
  assert(resolved.get('EX:plantar-roll').parts.length);
  const selectedRole=api.guestNodeChoice('physical-exercise',resolved.get('ROLE:WGER:203:primary'));
  const ancestors=api.hierarchyAncestors(api.hierarchyNavigation(atlas,'guest:physical-exercise',hierarchy),selectedRole);
  assert(ancestors.some(node=>node.id==='hierarchy:guest:physical-exercise:WGER:203'),'Shared muscle sets must preserve the selected exercise path');
  assert(resolved.get('EX:single-leg-balance').parts.some(p=>p.system==='sensory'));
 }
 assert.equal(!!api.hierarchyNavigation(atlas,'guest:physical-exercise',hierarchy).children.length,roots.some(root=>root.parts.length));
 const layout=api.createHierarchicalExplosionLayout(atlas,input.hierarchyEntries(atlas),ids,'guest',1.7,roots);
 for(const stage of layout.stages){assert.equal(stage.clusters.flat().length,available.length);assert.equal(new Set(stage.clusters.flat()).size,available.length);assert(stage.positions.every(Number.isFinite));}
 console.log(`${file}: ${provenance.exercises.filter(e=>resolved.get(e.id).parts.length).length} exercises with anatomy; ${available.length} meshes assigned once`);
}
const state={selected:[],visible:[],guestQuery:'WGER:203',explode:0,isolate:false,view:'front',rotate:false,reset:0};
const url=new URL(api.viewUrl('http://localhost:3016','local-reference',state));url.searchParams.set('tree','guest:physical-exercise');
assert.equal(api.readViewUrl(url.search,input.catalogues[0].atlas,state).guestQuery,state.guestQuery);
console.log(`${provenance.exercises.length} exercise records, roles, source credits, laterality, URL restoration and shared anatomy validated.`);
