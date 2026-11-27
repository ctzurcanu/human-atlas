import test from 'node:test';
import assert from 'node:assert/strict';
import {diaphragmControlRegistry,diaphragmTermParts,diaphragmLabelReference} from './diaphragm-controls.ts';

const sha='a'.repeat(64),plan='b'.repeat(64);
const references=new Set([10,11,14,15,16]);
function fixture(){
 const terms=Array.from({length:17},(_,i)=>`A04.4.02.${String(i+1).padStart(3,'0')}`);
 const parts=terms.filter((_,i)=>i!==0&&!references.has(i+1)).map(term=>({id:term,conceptId:`concept:${term}`,ta98Term:term,
  pairedTissue:true,sectionSolid:true,normalSigns:12,chunk:0,indexCount:12,bounds:[[-1,-1,-1],[1,1,1]]}));
 const ids=parts.map(p=>p.id);
 const boundary=name=>({name,closed:true,points:[[.2,0,.1],[.3,0,.2],[.4,0,.1]],samplePartIds:Array(3).fill(terms[16]),
  faceIds:[0,1,2],barycentric:Array.from({length:3},()=>[1,0,0]),sourceBankSha256:sha});
 const groups=terms.map((term,i)=>({term,conceptId:`concept:${term}`,partIds:i===0?ids:references.has(i+1)?[terms[16]]:[term],
  representation:i===0?'whole-organ':references.has(i+1)?'boundary-reference':'finite-tissue',
  ...(references.has(i+1)?{boundaries:[15,16].includes(i+1)?[boundary('left'),boundary('right')]:[boundary('rim')]}:{})}));
 groups[1].representation='tissue-assembly';groups[1].partIds=parts.filter(p=>terms.slice(1,7).includes(p.ta98Term)).map(p=>p.id);
 return {sex:'male',parts,chunks:[{sha256:sha}],concepts:groups.map(g=>({id:g.conceptId,ta98Kind:'exact',ta98Term:g.term,elements:g.partIds})),
  finiteDiaphragmIntegration:{schema:'diaphragm-complete-local-v7',partIds:ids,bankSha256:sha,sourcePlanSha256:plan,termGroups:groups}};
}

test('all seventeen terms bind actual finite tissue or an explicit non-solid reference',()=>{
 const atlas=fixture(),terms=diaphragmTermParts(atlas);
 assert.equal(terms.size,17);
 assert.equal(terms.get('A04.4.02.001').length,atlas.parts.length);
 assert.equal(terms.get('A04.4.02.002').length,6);
 assert.equal(atlas.parts.some(p=>p.ta98Term==='A04.4.02.011'),false);
 assert.equal(terms.get('A04.4.02.011')[0].ta98Term,'A04.4.02.017');
});
test('opening labels anchor to their actual rim witness and retain ordinary selection members',()=>{
 const atlas=fixture(),g=atlas.finiteDiaphragmIntegration.termGroups[10];
 const anchor=diaphragmLabelReference(atlas,{id:g.conceptId,ta98Term:g.term,elements:g.partIds});
 assert.deepEqual(anchor,{partId:'A04.4.02.017',point:[.3,0,.2]});
 assert.equal(diaphragmLabelReference(atlas,{id:g.conceptId,ta98Term:g.term,elements:atlas.parts.map(p=>p.id)}),undefined);
});
test('stale banks, suppressed tissues, omitted terms and corrupt rim witnesses are rejected',()=>{
 const mutations=[a=>a.chunks[0].sha256=plan,a=>a.parts[0].suppressed=true,
  a=>a.finiteDiaphragmIntegration.termGroups.pop(),
  a=>a.finiteDiaphragmIntegration.termGroups[10].boundaries[0].faceIds[0]=99,
  a=>a.finiteDiaphragmIntegration.termGroups[10].boundaries[0].barycentric[0]=[-1,1,1],
  a=>a.finiteDiaphragmIntegration.termGroups[10].boundaries[0].sourceBankSha256=plan,
  a=>a.concepts[10].elements=['unrelated'],
  a=>a.finiteDiaphragmIntegration.termGroups[14].boundaries.pop()];
 for(const change of mutations){const atlas=fixture();change(atlas);assert.equal(diaphragmControlRegistry(atlas),null);}
});
test('existing catalogues without the new complete marker retain their ordinary controls',()=>{
 const atlas=fixture();atlas.finiteDiaphragmIntegration.schema='shared-finite-tissues-left-medial';
 assert.equal(diaphragmTermParts(atlas).size,0);
});
