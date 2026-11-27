import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';

const built=await build({stdin:{contents:"export * from './app/section-closure';export {Matrix4} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {filterSectionClosureCandidates:filter,Matrix4}=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const ids=['remainder','left-facet','right-facet'];
function fixture(){return ids.map((id,index)=>({part:{id,sectionAssembly:'shared-lamina',sectionClosureParts:[...ids],sectionSolid:false},index,cutIndex:0,matrix:new Matrix4()}));}
const selected=rows=>rows.map(row=>row.part.id);

test('legacy cap candidates retain identity/order without reading matrices',()=>{
 const legacy=[{part:{id:'old-open-assembly',sectionAssembly:'old'},cutIndex:0},{part:{id:'old-solid'},cutIndex:1}];
 const result=filter(legacy,()=>{throw Error('Legacy matrices must not be queried');});
 assert.deepEqual(result,legacy);result.forEach((row,index)=>assert.equal(row,legacy[index]));
});
test('complete reciprocal closure retains open patches without relabeling them as solids',()=>{
 const rows=fixture();rows[1].part.sectionClosureParts.reverse();
 const before=JSON.stringify(rows),result=filter(rows);
 assert.deepEqual(selected(result),ids);result.forEach((row,index)=>assert.equal(row,rows[index]));
 assert.equal(JSON.stringify(rows),before);assert.ok(result.every(row=>row.part.sectionSolid===false));
});
test('isolated or hidden required patches block all remaining declared candidates',()=>{
 const rows=fixture();
 for(let missing=0;missing<rows.length;missing++)assert.deepEqual(filter(rows.filter((_,index)=>index!==missing)),[]);
 assert.deepEqual(filter([rows[0]]),[]);
});
test('a peer visible for another cut cannot complete a closure',()=>{
 const rows=fixture();rows[2].cutIndex=1;assert.deepEqual(filter(rows),[]);
 const both=[...fixture(),...fixture().map(row=>({...row,cutIndex:1}))];
 assert.equal(filter(both).length,6);
 assert.deepEqual(filter(both.filter(row=>!(row.cutIndex===1&&row.part.id==='right-facet'))),both.slice(0,3));
});
test('mismatched declarations block every affected patch, including subsets',()=>{
 const rows=fixture();rows[1].part.sectionClosureParts=['remainder','left-facet'];
 assert.deepEqual(filter(rows),[]);
 rows[1].part.sectionClosureParts=[...ids,'missing-peer'];assert.deepEqual(filter(rows),[]);
});
test('mixed or missing assembly identity cannot certify a closure',()=>{
 for(const assembly of ['another-lamina',undefined,'']){
  const rows=fixture();rows[1].part.sectionAssembly=assembly;assert.deepEqual(filter(rows),[]);
 }
});
test('transforms must match exactly and are reread after an explosion change',()=>{
 const rows=fixture();assert.equal(filter(rows).length,3);
 rows[1].matrix.elements[12]=1e-12;assert.deepEqual(filter(rows),[],'No transform tolerance may join separately moved patches');
 rows[1].matrix.identity();assert.equal(filter(rows).length,3,'No stale transform eligibility cache');
 rows.forEach(row=>row.matrix.makeRotationX(.2));assert.equal(filter(rows).length,3);
 rows[2].matrix.elements[0]=Number.NaN;assert.deepEqual(filter(rows),[]);
});
test('raw scene candidates can supply picker world matrices before GPU or CPU work',()=>{
 const withMatrix=fixture(),rows=withMatrix.map(({matrix,...row})=>row);
 assert.deepEqual(filter(rows),[],'Declared closure requires actual transforms');
 assert.deepEqual(filter(rows,row=>withMatrix[row.index].matrix),rows);
 withMatrix[1].matrix.makeTranslation(0,.001,0);assert.deepEqual(filter(rows,row=>withMatrix[row.index].matrix),[]);
});
test('malformed closure lists and ambiguous duplicate visible peers are rejected',()=>{
 for(const list of [[],['remainder','remainder'],['left-facet','right-facet'],['remainder',''],null]){
  const rows=fixture();rows[0].part.sectionClosureParts=list;assert.deepEqual(filter(rows),[]);
 }
 const rows=fixture();assert.deepEqual(filter([...rows,rows[1]]),[]);
});
test('a missing declaration cannot bypass closure through the legacy cap path',()=>{
 const rows=fixture();delete rows[1].part.sectionClosureParts;
 assert.deepEqual(filter(rows),[],'An undeclared patch in the same assembly must not receive a standalone cap');
 const independent={part:{id:'legacy',sectionAssembly:'legacy'},cutIndex:0};
 assert.deepEqual(filter([rows[0],independent]),[independent]);
});
