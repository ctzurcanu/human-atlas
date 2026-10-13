import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const plan=JSON.parse(await fs.readFile('app/data/colon-wall-plan.json'));
const bundle=await build({stdin:{contents:"export * from './app/guest-hierarchy';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {parseGuestHierarchy,resolveGuestHierarchy}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const raw=JSON.parse(await fs.readFile('public/assets/ta98.json')),hierarchy=parseGuestHierarchy(raw),graph=new Map(raw.nodes.map(n=>[n.id,n]));
test('Colon wall preserves the supplied structure and canonical authority parentage',()=>{
 assert.deepEqual(graph.get(plan.owner).childIds.filter(id=>id.startsWith('CW:')),plan.innerToOuter.map(key=>'CW:colon:'+key));
 for(const [key,spec] of Object.entries(plan.terms)){const node=graph.get('CW:colon:'+key);assert.ok(node);assert.deepEqual(node.childIds??[],spec.children.map(child=>'CW:colon:'+child));assert.deepEqual(node.matches,['CW:colon:'+key,...(spec.sourcePartIds??[])]);if(spec.ta98)assert.equal(node.links[0].id,spec.ta98);}
 assert.ok(graph.get('CW:colon:subserosa').childIds.includes('CW:colon:epiploic-appendices'));
 assert.deepEqual(graph.get('CW:colon:longitudinal-muscle').childIds,['CW:colon:taeniae']);
 assert.ok(graph.get('A05.7.03.001').childIds.includes('A05.7.03.010'),'canonical appendix parent is retained');
});
for(const sex of ['female','male'])test(`${sex}: an outer colon mesh does not become a represented layer or gland`,async()=>{
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-${sex}.json`));const nodes=new Map();function walk(n){nodes.set(n.id,n);n.children.forEach(walk)}resolveGuestHierarchy(atlas,hierarchy).forEach(walk);
 const expected=key=>[...(plan.terms[key].sourcePartIds??[]),...plan.terms[key].children.flatMap(expected)];
 for(const key of Object.keys(plan.terms)){const node=nodes.get('CW:colon:'+key);assert.ok(node);assert.deepEqual(node.parts.map(p=>p.id).sort(),sex==='male'?expected(key).sort():[]);assert.deepEqual(node.directParts.map(p=>p.id),sex==='male'?(plan.terms[key].sourcePartIds??[]):[]);if(node.directParts.length)assert.match(node.description,/Open band surfaces/);else if(sex==='male'&&key==='taeniae')assert.match(node.description,/All three named native open band references/);else assert.match(node.description,/Missing separate representation/);}
 assert.ok(!nodes.get('CW:colon:longitudinal-muscle').directParts.length,'bands do not imply a complete longitudinal coat');
 assert.ok(nodes.get('A05.7.03.001').parts.length>0);
});
