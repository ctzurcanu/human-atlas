import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
const compiled=await build({stdin:{contents:"export * from './app/ta98-term-bindings';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {catalogueTa98TermParts}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
function fixture(){
 const owner={id:'whole',chunk:13,positions:100,normals:200,uvs:300,vertexCount:24,indices:400,indexCount:36};
 const part={...owner,id:'ascending',indices:544,indexCount:36,conceptId:'branch',ta98Term:'A12.2.16.028',sourceOwner:'whole',sourceOuterSurfaceOnly:true,sourceLumenPatencyUnverified:true};
 return {parts:[owner,part],concepts:[{id:'branch',elements:['ascending'],ta98Term:part.ta98Term,ta98Kind:'exact',coverageLimitation:'Original outside surface; wall and patency unverified.'}],ta98TermBindings:[{term:part.ta98Term,conceptId:'branch',partIds:['ascending'],representation:'source-surface-reference'}]};
}
test('explicit original source span is selectable without pretending to be finite tissue',()=>{const a=fixture();assert.equal(catalogueTa98TermParts(a).get('A12.2.16.028')[0].id,'ascending');delete a.ta98TermBindings[0].representation;assert.equal(catalogueTa98TermParts(a).size,0);});
test('source references require qualifiers, literal source storage and matching world transform',()=>{
 for(const edit of [a=>delete a.parts[1].sourceOuterSurfaceOnly,a=>delete a.parts[1].sourceLumenPatencyUnverified,a=>a.parts[1].indices++,a=>a.parts[1].chunk++,a=>a.parts[1].sourceOffset=[1,0,0],a=>a.concepts[0].coverageLimitation='',a=>a.parts[1].suppressed=true]){const a=fixture();edit(a);assert.equal(catalogueTa98TermParts(a).size,0);}
});
test('existing finite binding still requires an actual section solid',()=>{const a=fixture();delete a.ta98TermBindings[0].representation;a.parts[1].sectionSolid=true;assert.equal(catalogueTa98TermParts(a).size,1);});
