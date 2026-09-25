import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/anatomy-path.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {anatomyPath}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);
const atlas=JSON.parse(readFileSync('public/models/atlas-male-complete.json'));
const concept=atlas.concepts.find(item=>item.name==='Pectoral region (left)');
assert.ok(concept);
const parts=atlas.parts.filter(part=>concept.elements.includes(part.id));
const labels=anatomyPath(concept.name,parts,atlas.parts,'male',atlas.scope).map(node=>node.label);
assert.deepEqual(labels,['Human Male','Body Surface','Trunk','Regions of trunk']);
console.log('Anatomical breadcrumbs contain ancestors without repeating the selected region.');
