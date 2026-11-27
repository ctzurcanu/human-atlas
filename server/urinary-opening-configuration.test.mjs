import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/urinary-wall';export * from './app/anatomy-label-plan';export * from './app/requested-selection-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const read=()=>JSON.parse(fs.readFileSync(process.env.ATLAS_URINARY_CANDIDATE??'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'));
const sex=read().sex,term=sex==='male'?'A09.4.02.002':'A09.2.03.002',id='URINARY:'+sex+':bladder:'+term;
function fixture(){
 const atlas=read(),m=atlas.finiteUrinaryIntegration,parent=m.termGroups.find(g=>g.term===(sex==='male'?'A08.3.01.027M':'A08.3.01.028F'));
 m.termGroups=m.termGroups.filter(g=>g.term!==term);atlas.concepts=atlas.concepts.filter(c=>c.id!==id);
 const group={term,conceptId:id,partIds:parent.partIds,openingConfiguration:{state:'voiding',openingTerm:parent.term,sourcePlanSha256:m.openingPlanSha256,sourceBankSha256:m.charts.find(c=>c.name==='bladder').bankSha256}};
 const choice={id,name:'Emptying internal urethral orifice',elements:parent.partIds,partOf:parent.conceptId,ta98Term:term,ta98Kind:'exact'};
 m.termGroups.push(group);atlas.concepts.push(choice);return{atlas,group,parent,choice};
}
test('Voiding uses the existing clinical neck rim with ordinary label and URL selection',()=>{
 const {atlas,group,parent,choice}=fixture();assert.ok(api.urinaryWallRegistry(atlas));
 const anchor=api.urinaryLabelReference(atlas,choice),b=parent.openingBoundaries[0];assert.ok(anchor);
 assert.ok(b.points.some((p,i)=>b.samplePartIds[i]===anchor.partId&&p.every((x,k)=>x===anchor.point[k])));
 const labels=api.selectedAnatomyLabels(atlas,choice.elements,choice,false);assert.equal(labels.length,1);
 assert.equal(labels[0].name,choice.name);assert.deepEqual(labels[0].anchor,anchor);
 assert.deepEqual(api.anatomyLabelChoice(labels[0]).elements,choice.elements);
 assert.equal(api.requestedSelectionChoice(atlas,[id],choice.elements).id,id);
 assert.ok(!api.layerLabelTargets(labels,atlas.parts).has(id));
});
test('A closed filling state cannot be claimed by the installed open funnel',()=>{
 const {atlas,group}=fixture();group.openingConfiguration.state='filling';assert.equal(api.urinaryWallRegistry(atlas),null);
 const other=fixture();other.group.term=sex==='male'?'A09.4.02.001':'A09.2.03.001';other.choice.ta98Term=other.group.term;assert.equal(api.urinaryWallRegistry(other.atlas),null);
});
test('Changing the source bank, source plan or opening owner rejects a configuration',()=>{
 for(const mutate of [g=>g.openingConfiguration.sourceBankSha256='f'.repeat(64),g=>g.openingConfiguration.sourcePlanSha256='e'.repeat(64),g=>g.partIds=['BLADDER:'+sex+':body:mucosa']]){
  const {atlas,group}=fixture();mutate(group);assert.equal(api.urinaryWallRegistry(atlas),null);
 }
});

test('A voiding configuration must use the model sex and its own neck opening',()=>{
 const opposite=sex==='male'?{term:'A09.2.03.002',opening:'A08.3.01.028F'}:{term:'A09.4.02.002',opening:'A08.3.01.027M'};
 const wrongTerm=fixture();wrongTerm.group.term=opposite.term;assert.equal(api.urinaryWallRegistry(wrongTerm.atlas),null);
 const wrongOpening=fixture();wrongOpening.group.openingConfiguration.openingTerm=opposite.opening;assert.equal(api.urinaryWallRegistry(wrongOpening.atlas),null);
});

test('Stored voiding control retains its explicit configuration identity',t=>{
 const atlas=read(),g=atlas.finiteUrinaryIntegration.termGroups.find(g=>g.term===term);
 if(!g){t.skip('Legacy atlas before configuration installation');return;}
 assert.equal(g.openingConfiguration.state,'voiding');assert.ok(api.urinaryWallRegistry(atlas));
 const c=atlas.concepts.find(c=>c.id===g.conceptId);assert.ok(api.urinaryLabelReference(atlas,c));
});

test('Female transferred coats peel completely and retain female source selection',t=>{
 const atlas=read();if(atlas.sex!=='female'){t.skip('Female transfer fixture');return;}
 const marker=atlas.finiteUrinaryIntegration;assert.equal(marker.partIds.length,97);
 for(const code of ['A08.3.01.001','A08.4.01.001F']){
  const group=marker.termGroups.find(g=>g.term===code),choice=atlas.concepts.find(c=>c.id===group.conceptId);
  const steps=api.urinaryTissuePeelSteps(atlas,choice);assert.ok(steps);
  assert.deepEqual(new Set(steps.flatMap(s=>s.ids)),new Set(choice.elements));
  assert.equal(steps.flatMap(s=>s.ids).length,choice.elements.length,'Every real tissue is peeled exactly once');
  if(code==='A08.4.01.001F')assert.ok(steps.some(s=>s.name==='Spongy layer'));
 }
 assert.ok(marker.partIds.every(id=>!id.includes(':male:')&&!id.includes('prostate')));
 for(const [source,chart] of [['LOCAL:female:bladder:187','bladder'],['LOCAL:female:urethra:188','urethra']]){
  const p=atlas.parts.find(p=>p.id===source);assert.ok(p.suppressed);
  assert.deepEqual(new Set(p.replacementParts),new Set(marker.charts.find(c=>c.name===chart).partIds));
 }
});
