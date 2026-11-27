import assert from 'node:assert/strict';
import {readFileSync,existsSync} from 'node:fs';
import test from 'node:test';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const bodyFile='.local-models/organ-models/spinal-cord/male-body.json';
const masterFile='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
if(!existsSync(bodyFile)||!existsSync(masterFile)){
 test('Lingual component controls require the installed private anatomy',{skip:'Local geometry is intentionally absent'},()=>{});
}else{
const bytes=readFileSync(bodyFile),atlas=JSON.parse(bytes),parts=new Map(atlas.parts.map(p=>[p.id,p]));
const masterBytes=readFileSync(masterFile),master=JSON.parse(masterBytes),source=new Map(master.parts.map(p=>[p.id,p]));
const sha=b=>createHash('sha256').update(b).digest('hex');
const built=await build({stdin:{contents:`
 export {nativePartFrame,ta98ModeledParent} from './app/ta98-modeled-parent.ts';
 export {anatomicalComponentIds} from './app/anatomical-components.ts';
 export {hierarchyEntries,majorSystemFor,systemPathFor,buildAnatomyNodes} from './app/anatomy-hierarchy.ts';
 export {buildSystemTree} from './app/system-tree.tsx';
 export {withLingualArterialSystem} from './app/lingual-arterial-system.ts';
 export {registeredLingualFamily} from './app/lingual-arterial-relations.ts';
 `,resolveDir:process.cwd(),sourcefile:'lingual-controls-test.ts'},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import(`data:text/javascript;base64,${Buffer.from(built.outputFiles[0].contents).toString('base64')}`);
const concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
const families=['left','right'].map(side=>api.registeredLingualFamily(atlas.parts,side));
assert.ok(families.every(Boolean),'Both installed seven-cell native banks must be present');
const owned=families.flatMap(family=>[...family.values()]);
assert.equal(owned.length,12);
const ownedIds=new Set(owned.map(p=>p.id));
for(const part of owned){
 assert.deepEqual({...part,chunk:atlas.sourceChunkIndices[part.chunk]},source.get(part.id),'Ordinary model must use the installed source records');
 assert.equal(api.nativePartFrame(part),'ZA:','Reviewed arterial cell is registered in the native male frame');
}
for(const side of ['left','right']){
 const family=families[side==='left'?0:1],trunk=family.get('lingual'),whole=concepts.get(trunk.conceptId);
 assert.equal(api.ta98ModeledParent(whole,atlas.concepts,parts)?.id,`ZA:External carotid artery.${side[0]}`);
 for(const part of family.values())if(part!==trunk){
  const child=concepts.get(part.conceptId);
  assert.equal(child.partOf,whole.id);
  assert.equal(api.ta98ModeledParent(child,atlas.concepts,parts)?.id,whole.id,'Every authored branch resolves to its ipsilateral parent');
 }
 assert.deepEqual(new Set(api.anatomicalComponentIds(whole,[],parts,atlas.concepts)),new Set([...family.values()].map(p=>p.id)),
  'Show components retains the owned trunk and all five child cells, excluding ECA and opposite side');
 assert.deepEqual(api.anatomicalComponentIds({...whole,id:'selection-set'},[],parts,atlas.concepts),[]);
}

const entries=api.hierarchyEntries(atlas),system=api.buildSystemTree(entries,'systems',atlas).find(b=>b.id==='circulatory');
const hits=[];
const visit=nodes=>{for(const node of nodes){if(node.name==='Lingual artery')hits.push(node);if(node.kind==='group')visit(node.nodes);}};
visit(system.nodes);
assert.equal(hits.length,1,'Systems exposes one complete Lingual artery control, not separate 10-cell and 2-cell controls');
const whole=hits[0];
assert.equal(whole.kind,'group');assert.equal(whole.itemCount,5);
assert.deepEqual(new Set(whole.parts.map(p=>p.id)),ownedIds);
assert.equal(whole.nodes.length,5);
const trunk=whole.nodes.find(n=>n.name==='Main trunk');
assert.equal(trunk.kind,'bilateral');assert.equal(trunk.parts.length,2);
assert.deepEqual(new Set(trunk.parts.map(p=>p.ta98Term)),new Set(['A12.2.05.015']));
assert.deepEqual(new Set(whole.nodes.flatMap(n=>n.parts.map(p=>p.id))),ownedIds);

const plain=api.buildAnatomyNodes(entries.filter(e=>api.majorSystemFor(e)==='circulatory'),e=>api.systemPathFor(e,atlas.scope));
const changed=api.withLingualArterialSystem(atlas,plain);
const unrelated=nodes=>nodes.flatMap(n=>{
 const ids=n.parts.filter(p=>!ownedIds.has(p.id)).map(p=>p.id).sort();
 if(!ids.length)return [];
 return [{id:n.id,name:n.name,kind:n.kind,parts:ids,...(n.kind==='group'?{nodes:unrelated(n.nodes)}:{})}];
});
assert.deepEqual(unrelated(changed),unrelated(plain),'Unrelated Systems branches and identities stay unchanged');

const sample=families[0].get('deep-lingual');
for(const [label,part] of [
 ['wrong term',{...sample,ta98Term:'A12.2.05.018'}],
 ['wrong side',{...sample,sectionAssembly:'right-lingual-arterial-family:male'}],
 ['wrong revision',{...sample,lingualArterialWallRevision:'male-right-lingual-finite-network-v2'}],
 ['offset donor',{...sample,sourceOffset:[0,0,0]}],
 ['instance donor',{...sample,geometryInstance:{}}],
 ['wrong provenance',{...sample,provenance:{url:'https://example.com/unregistered'}}],
 ['wrong source',{...sample,sourceId:'unregistered'}],
 ['unpaired',{...sample,pairedTissue:false}],
 ['suppressed',{...sample,suppressed:true}],
])assert.equal(api.nativePartFrame(part),undefined,label);
const incomplete={...atlas,parts:atlas.parts.filter(p=>p.id!==sample.id)};
assert.equal(api.withLingualArterialSystem(incomplete,plain),plain,'Incomplete native banks must not be promoted to whole assemblies');
const wrongParent={...atlas,concepts:atlas.concepts.map(c=>c.id===sample.conceptId?{...c,partOf:'ZA:External carotid artery.r'}:c)};
assert.equal(api.withLingualArterialSystem(wrongParent,plain),plain,'Incorrect declared parents must not be hidden by Systems aggregation');
assert.equal(sha(readFileSync(bodyFile)),sha(bytes));
assert.equal(sha(readFileSync(masterFile)),sha(masterBytes));
console.log(JSON.stringify({passed:true,source_sha256:sha(masterBytes),body_sha256:sha(bytes),registered_cells:12,whole_controls:1,component_cells_per_side:6,authority_terms:5,unrelated_systems_preserved:true,invalid_frame_cases_rejected:9}));

}
