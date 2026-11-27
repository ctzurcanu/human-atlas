import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {existsSync,readFileSync} from 'node:fs';
import test from 'node:test';
import {build} from 'esbuild';

// This is a metadata/selection regression. It neither loads nor qualifies meshes.
const base='.local-models/geometry-workbench/circulatory-completion/male-calcarine-source-partition-v1';
const candidatePath=process.env.ATLAS_CALCARINE_CANDIDATE??`${base}/proposed-atlas.json`;
const hierarchyPath=process.env.ATLAS_CALCARINE_HIERARCHY??`${base}/projected-hierarchy/public/assets/ta98.json`;
const navigationCandidatePath=process.env.ATLAS_CALCARINE_NAVIGATION_CANDIDATE??`${base}/proposed-navigation-fixed-atlas.json`;
const sides=['left','right'];
const ids=side=>({parent:`ZA:Medial occipital artery.${side[0]}`,calcarine:`TA98SOURCE:male:calcarine-branch.${side[0]}`,parieto:`ZA:Parieto-occipital artery.${side[0]}`});
const terms={parent:'A12.2.07.098',parieto:'A12.2.07.101',calcarine:'A12.2.07.102'};
const hash=bytes=>createHash('sha256').update(bytes).digest('hex');
const built=await build({stdin:{contents:`
 export {anatomicalComponentIds,storedAssemblyParent} from './app/anatomical-components.ts';
 export {hierarchyEntries,systemPathFor} from './app/anatomy-hierarchy.ts';
 export {buildSystemTree} from './app/system-tree.tsx';
 export {parseGuestHierarchy,resolveGuestHierarchy} from './app/guest-hierarchy.ts';
 export {guestNodeChoice,guestDirectAnatomyParts} from './app/guest-choice.ts';
 export {selectedAnatomyLabels,anatomyLabelChoice} from './app/anatomy-label-plan.ts';
 export {displayLaterality} from './app/laterality.ts';
 export {hierarchyNavigation,hierarchyAncestors} from './app/hierarchy-navigation.ts';
 export {localDescription} from './app/structure-description.ts';
`,resolveDir:process.cwd(),sourcefile:'calcarine-source-control-check.ts'},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import(`data:text/javascript;base64,${Buffer.from(built.outputFiles[0].contents).toString('base64')}`);

function fixture(){
 const parts=[],concepts=[];
 for(const [index,side]of sides.entries()){
  const key=ids(side),members=[key.parent,key.parieto];
  concepts.push({id:key.parent,name:`Medial occipital artery (${side})`,elements:members,componentElements:members,ta98Term:terms.parent,ta98Kind:'exact'},
   {id:key.calcarine,name:`Calcarine branch (${side})`,elements:[key.parent],partOf:key.parent,ta98Term:terms.calcarine,ta98Kind:'exact'},
   {id:key.parieto,name:`Parieto-occipital artery (${side})`,elements:[key.parieto],partOf:key.parent,ta98Term:terms.parieto,ta98Kind:'exact'});
  for(const [i,[id,conceptId,name,term]]of [[key.parent,key.calcarine,`Calcarine branch (${side})`,terms.calcarine],[key.parieto,key.parieto,`Parieto-occipital artery (${side})`,terms.parieto]].entries()){
   parts.push({id,conceptId,name,ta98Term:term,system:'arterial',groups:['Head'],chunk:0,positions:(index*2+i)*36,normals:0,indices:0,vertexCount:3,indexCount:3,bounds:[[index?-.01:.001,1.72,-.10],[index?-.001:.01,1.75,-.03]]});
  }
 }
 return {version:'calcarine-selection-fixture/v1',sex:'male',parts,concepts,chunks:[],triangles:4};
}

// The generic resolver needs explicit selectors for a newly introduced source
// concept. A ta98Term tag alone does not add it to a static hierarchy asset.
function fixtureHierarchy(){
 const source=JSON.parse(readFileSync('public/assets/ta98.json'));
 const nodes=source.nodes.filter(node=>Object.values(terms).includes(node.id)).map(node=>({...node,childIds:node.id===terms.parent?[terms.parieto,terms.calcarine]:undefined,
  matches:sides.map(side=>ids(side)[Object.entries(terms).find(([,term])=>term===node.id)[0]])}));
 return {schema:'human-atlas-hierarchy/v2',id:'ta98',name:'TA98 metadata fixture',roots:[terms.parent],nodes};
}

function focusedCatalogue(atlas){
 const partIds=new Set(sides.flatMap(side=>[ids(side).parent,ids(side).parieto]));
 const conceptIds=new Set(sides.flatMap(side=>Object.values(ids(side))));
 return {version:atlas.version,sex:atlas.sex,parts:atlas.parts.filter(part=>partIds.has(part.id)),concepts:atlas.concepts.filter(concept=>conceptIds.has(concept.id)),chunks:atlas.chunks,triangles:atlas.triangles};
}

function assertControls(atlas,hierarchy){
 const before=JSON.stringify(atlas),parts=new Map(atlas.parts.map(part=>[part.id,part])),concepts=new Map(atlas.concepts.map(concept=>[concept.id,concept]));
 assert.equal(parts.size,4,'There are four physical sources, with no new parent mesh');
 for(const side of sides){
  const key=ids(side),parent=concepts.get(key.parent),calcarine=concepts.get(key.calcarine),parieto=concepts.get(key.parieto),expected=[key.parent,key.parieto];
  assert.equal(parent?.ta98Term,terms.parent);assert.equal(parent.ta98Kind,'exact');
  assert.deepEqual(new Set(parent.elements),new Set(expected));assert.deepEqual(new Set(parent.componentElements),new Set(expected));
  assert.ok(atlas.parts.every(part=>part.conceptId!==parent.id),'The parent is a stored union without its own geometry');
  assert.deepEqual(new Set(api.anatomicalComponentIds(parent,[],parts,atlas.concepts)),new Set(expected));
  for(const [child,term,partId]of [[calcarine,terms.calcarine,key.parent],[parieto,terms.parieto,key.parieto]]){
   assert.equal(child?.ta98Term,term);assert.equal(child.ta98Kind,'exact');assert.equal(child.partOf,parent.id);assert.deepEqual(child.elements,[partId]);
   assert.equal(parts.get(partId)?.conceptId,child.id);assert.equal(api.displayLaterality(child.name).side,side);
   assert.equal(api.storedAssemblyParent(child,atlas.concepts,parts)?.id,parent.id);
   assert.deepEqual(api.anatomicalComponentIds(child,[],parts,atlas.concepts),[],'A leaf does not expand to its sibling');
  }
  assert.equal(calcarine.name,`Calcarine branch (${side})`);
  const expanded=api.selectedAnatomyLabels(atlas,parent.elements,parent,true);
  assert.deepEqual(new Set(expanded.map(label=>label.id)),new Set([calcarine.id,parieto.id]));
  for(const label of expanded)assert.deepEqual(api.anatomyLabelChoice(label).elements,concepts.get(label.id).elements,'A callout selects only its named source');
  const collapsed=api.selectedAnatomyLabels(atlas,parent.elements,parent,false);
  assert.equal(collapsed.length,1);assert.equal(collapsed[0].id,parent.id);assert.deepEqual(new Set(collapsed[0].components),new Set(expected));
 }
 const entries=api.hierarchyEntries(atlas);assert.equal(entries.length,4);
 for(const side of sides){
  const key=ids(side),calcarine=entries.find(entry=>entry.id===key.calcarine),parieto=entries.find(entry=>entry.id===key.parieto);
  assert.ok(!entries.some(entry=>entry.id===key.parent));
  assert.equal(calcarine.name,`Calcarine branch (${side})`);assert.equal(calcarine.terminology.ta98,terms.calcarine);assert.equal(parieto.terminology.ta98,terms.parieto);
  assert.ok(api.systemPathFor(calcarine).includes('Medial occipital artery'));assert.ok(api.systemPathFor(parieto).includes('Medial occipital artery'));
 }
 const tree=api.buildSystemTree(entries,'systems',atlas),leaves=[],groups=[];
 const walkSystems=nodes=>{for(const node of nodes){if(node.kind==='group'){groups.push(node);walkSystems(node.nodes);}else if(node.kind==='bilateral')leaves.push(...node.entries);else leaves.push(node.entry);}};
 tree.forEach(branch=>walkSystems(branch.nodes));
 assert.equal(tree.length,1);assert.equal(tree[0].id,'circulatory');
 assert.deepEqual(leaves.flatMap(entry=>entry.parts.map(part=>part.id)).sort(),[...parts.keys()].sort(),'Each physical source occurs exactly once in Systems');
 const parents=groups.filter(group=>group.terminology.ta98===terms.parent);
 assert.ok(parents.length,'Both branch families retain medial occipital ancestry');
 // General location heuristics may place the two source names on different
 // Head/Skull paths. Each physical leaf must still occur in exactly one path.
 assert.deepEqual(parents.flatMap(group=>group.parts.map(part=>part.id)).sort(),[...parts.keys()].sort());
 const roots=api.resolveGuestHierarchy(atlas,api.parseGuestHierarchy(hierarchy)),nodes=new Map();
 const walkGuest=node=>{if(nodes.has(node.id))return;nodes.set(node.id,node);node.children.forEach(walkGuest);};roots.forEach(walkGuest);
 for(const [role,term]of Object.entries(terms)){
  const expected=sides.flatMap(side=>role==='parent'?[ids(side).parent,ids(side).parieto]:[role==='calcarine'?ids(side).parent:ids(side).parieto]);
  const node=nodes.get(term);assert.ok(node,term);assert.deepEqual(new Set(node.parts.map(part=>part.id)),new Set(expected),term);
  assert.deepEqual(new Set(api.guestNodeChoice('ta98',node).elements),new Set(expected),term);
 }
 assert.deepEqual(api.guestDirectAnatomyParts('ta98',nodes.get(terms.parent)),[],'The TA98 parent does not repeat its child meshes as extra leaves');
 assert.deepEqual(new Set(api.anatomicalComponentIds(api.guestNodeChoice('ta98',nodes.get(terms.parent)),[],parts,atlas.concepts)),new Set(parts.keys()));
 for(const role of ['parieto','calcarine']){
  const node=nodes.get(terms[role]);assert.equal(node.directParts.length,2);
  const children=api.guestNodeChoice('ta98',node).children;
  assert.deepEqual(new Set(children.map(child=>child.id)),new Set(sides.map(side=>ids(side)[role])));
  for(const side of sides){
   const child=children.find(child=>child.id===ids(side)[role]);
   assert.equal(child.name,concepts.get(child.id).name);assert.equal(api.displayLaterality(child.name).side,side);assert.deepEqual(child.elements,concepts.get(child.id).elements);
  }
 }
 assert.equal(JSON.stringify(atlas),before,'Controls leave the metadata fixture unchanged');
}

test('meshless medial occipital parents expose two native source components with distinct leaf labels',()=>{
 assertControls(fixture(),fixtureHierarchy());
});

test('the local calcarine proposal and regenerated TA98 projection retain the same controls',{skip:!existsSync(candidatePath)},()=>{
 assert.ok(existsSync(hierarchyPath),'Run the read-only redirected hierarchy builder for the new source concept');
 const candidateBytes=readFileSync(candidatePath),hierarchyBytes=readFileSync(hierarchyPath);
 const candidate=JSON.parse(candidateBytes);
 assertControls(focusedCatalogue(candidate),JSON.parse(hierarchyBytes));
 // Also resolve against every unrelated stored concept: the exact-parent
 // index must not expand either side into a neighbouring source family.
 const parts=new Map(candidate.parts.map(part=>[part.id,part]));
 for(const side of sides){
  const key=ids(side),parent=candidate.concepts.find(concept=>concept.id===key.parent);
  assert.deepEqual(new Set(api.anatomicalComponentIds(parent,[],parts,candidate.concepts)),new Set([key.parent,key.parieto]));
  const labels=api.selectedAnatomyLabels(candidate,parent.elements,parent,true);
  assert.deepEqual(new Set(labels.map(label=>label.id)),new Set([key.calcarine,key.parieto]));
 }
 assert.equal(hash(readFileSync(candidatePath)),hash(candidateBytes));assert.equal(hash(readFileSync(hierarchyPath)),hash(hierarchyBytes));
});

test('cranial calcarine sources do not inherit false mediastinal membership or descriptions',{skip:!existsSync(navigationCandidatePath)},()=>{
 const raw=readFileSync(navigationCandidatePath),atlas=JSON.parse(raw);
 const hierarchy=api.parseGuestHierarchy(JSON.parse(readFileSync('public/assets/ta98.json')));
 const navigation=api.hierarchyNavigation(atlas,'guest:ta98',hierarchy),nodes=new Map(),pending=[navigation];
 while(pending.length){const node=pending.pop();if(nodes.has(node.id))continue;nodes.set(node.id,node);pending.push(...(node.children??[]));}
 for(const groupId of ['TA98GROUP:A07.1.02.101','TA98GROUP:A07.1.02.102']){
  const group=atlas.concepts.find(concept=>concept.id===groupId);assert.ok(group);
  assert.ok(sides.every(side=>!group.elements.includes(ids(side).parent)),'Mediastinal regions cannot own these cranial sources');
 }
 for(const side of sides){
  const key=ids(side),concept=atlas.concepts.find(record=>record.id===key.calcarine);
  // Label choices intentionally carry only identity/name/elements, as the UI does.
  const choice={id:concept.id,name:concept.name,elements:concept.elements};
  const ancestors=api.hierarchyAncestors(navigation,choice);
  assert.ok(ancestors.some(node=>node.id==='hierarchy:guest:ta98:A12.0.00.000'));
  assert.ok(ancestors.every(node=>!/^hierarchy:guest:ta98:A07\./.test(node.id)&&!/^trunk$|thorax|mediastinum/i.test(node.name)));
  const part=atlas.parts.find(record=>record.id===key.parent);
  const description=api.localDescription(choice.name,part,ancestors.map(node=>({label:node.name,search:node.name}))).text;
  assert.ok(!/trunk|thorax|thoracic|mediastinum/i.test(description),description);
  const canonical=nodes.get('hierarchy:guest:ta98:A12.2.07.102');assert.ok(canonical);
  assert.ok(canonical.children.some(child=>child.id===concept.id&&child.elements.length===1&&child.elements[0]===part.id),'Exact calcarine controls remain selectable');
 }
 assert.equal(hash(readFileSync(navigationCandidatePath)),hash(raw));
});
