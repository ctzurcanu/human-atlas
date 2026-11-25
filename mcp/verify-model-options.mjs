import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';

/** Compare an MCP response to the shipped package, independently of the MCP index. */
export async function verifyMaleTa98Options(options){
 const atlas=JSON.parse(await readFile(new URL('../public/models/male-ta98/atlas.json',import.meta.url),'utf8'));
 const hidden=new Set(atlas.defaultHidden);
 for(const id of atlas.selectableHidden??[])hidden.delete(id);
 const selectable=atlas.parts.filter(part=>!part.suppressed&&!hidden.has(part.id));
 assert.equal(options.model,'male-ta98');
 assert.deepEqual(new Set(options.systems.map(s=>s.id)),new Set(selectable.map(p=>p.system)));
 assert.equal(options.systems.reduce((total,s)=>total+s.pieces,0),selectable.length);
}
