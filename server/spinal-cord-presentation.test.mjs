import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createRequire} from 'node:module';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const bundle=await build({stdin:{contents:`
 export * from './app/spinal-cord-presentation';
 export * from './app/spinal-cord-controls';
 export * from './app/hierarchy-navigation';
 export * from './app/guest-hierarchy';
 export * from './app/anatomy-label-plan';
 export * from './app/preferred-source-representations';
 `,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const uiBundle=await build({stdin:{contents:`
 import React from 'react';
 import {renderToStaticMarkup} from 'react-dom/server';
 import Inspector from './app/selection-inspector';
 import {Sheet} from './components/ui/sheet';
 export function render(props){return renderToStaticMarkup(React.createElement(Sheet,{open:true},React.createElement(Inspector,props)));}
 `,resolveDir:process.cwd(),sourcefile:'spinal-inspector-test.tsx',loader:'tsx'},bundle:true,platform:'node',format:'cjs',write:false,define:{'import.meta.env.DEV':'false'}});
const ui={exports:{}};
new Function('require','module','exports',uiBundle.outputFiles[0].text)(createRequire(import.meta.url),ui,ui.exports);
const rawAtlas=JSON.parse(fs.readFileSync('.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'));
// Match page.tsx's loading pipeline: hidden old donors are suppressed before
// any hierarchy, label or details-panel controls are built.
const atlas=api.withHiddenReviewAlternatives(rawAtlas);
const before=createHash('sha256').update(JSON.stringify(atlas)).digest('hex');
const registry=api.spinalCordRegistry(atlas);
assert.ok(registry,'The installed female spinal assembly must be valid');
const whole=atlas.concepts.find(c=>c.id===registry.marker.wholeConceptId);
const byId=new Map(atlas.parts.map(p=>[p.id,p]));
const guest=api.parseGuestHierarchy(JSON.parse(fs.readFileSync('public/assets/ta98.json')));
const find=(choice,term)=>choice.terminology?.ta98===term||choice.ta98Term===term?choice:choice.children?.map(c=>find(c,term)).find(Boolean);
const noop=()=>{};

test('the webpage loading filter preserves verified canal controls while hiding the obsolete cord',()=>{
 assert.ok(registry.referenceParts.length>0);
 for(const p of registry.referenceParts){
  assert.ok(atlas.defaultHidden.includes(p.id),'Boundary controls start hidden');
  assert.ok(!p.suppressed,'Hidden boundary controls remain selectable');
 }
 assert.ok(byId.get(registry.marker.sourcePartId).suppressed,'Old whole-cord donor remains retired');
 const invalid=structuredClone(rawAtlas);
 invalid.parts.find(p=>p.id===registry.marker.referencePartIds[0]).spinalBoundaryReference.physicalOwnerRetained=false;
 assert.ok(api.withHiddenReviewAlternatives(invalid).parts.find(p=>p.id===registry.marker.referencePartIds[0]).suppressed,'Unverified references do not bypass duplicate hiding');
});

function inspect(choice){
 const componentChoices=api.spinalCordComponentChoices(atlas,choice);
 assert.ok(componentChoices,`${choice.name} must use the named component presentation`);
 assert.notEqual(api.spinalCordDisplayCount(choice.elements.map(id=>byId.get(id))),undefined);
 const html=ui.exports.render({titleRef:{current:null},choice,componentChoices,ancestors:[],
  selectedParts:choice.elements.map(id=>byId.get(id)),relationshipDepth:0,relationshipExhausted:false,
  partById:byId,concepts:atlas.concepts,state:{selected:choice.elements,view:'back',visible:['nervous'],explode:0,isolate:true,rotate:false,reset:0},
  viewUrl:'http://localhost:3016/',covering:[],depth:0,onDepth:noop,onCenter:noop,onIsolate:noop,
  onExpandRelationships:noop,onShowComponents:noop,onHide:noop,onHidePart:noop,onChoosePart:noop,onChooseChild:noop});
 assert.doesNotMatch(html,/class="member-list"|Included structures|[\d,]+ (?:modeled )?pieces/,'The actual details panel must not enumerate or count physical rendering pieces');
 assert.equal([...html.matchAll(/aria-label="Select /g)].length,componentChoices.length+(html.includes('aria-label="Anatomical parent"')?1:0));
 assert.ok(componentChoices.length<120,'Component rows must be named anatomical controls, not thousands of fragments');
 return {componentChoices,html};
}

for(const mode of ['direct','systems','guest:ta98'])test(`${mode}: rendered spinal details list named anatomy rather than thousands of fragments`,t=>{
 const choice=mode==='direct'?whole:find(api.hierarchyNavigation(atlas,mode,mode==='guest:ta98'?guest:undefined),'A14.1.02.001');
 assert.ok(choice);
 const {componentChoices,html}=inspect(choice);
 t.diagnostic(`${componentChoices.length} named component rows; ${choice.elements.length} underlying rendering pieces`);
 assert.ok(componentChoices.length>=75);
 assert.equal(new Set(componentChoices.map(c=>c.id)).size,componentChoices.length);
 assert.match(html,/Select Grey substance/);
 assert.match(html,/Select White substance/);
 const covered=new Set(componentChoices.flatMap(c=>c.elements));
 assert.ok(choice.elements.every(id=>covered.has(id)),'Every selected physical tissue remains reachable');
 const selected=new Set(choice.elements);
 for(const c of componentChoices)for(const id of c.elements){
  const p=byId.get(id);
  assert.ok(selected.has(id)||selected.has(p.spinalBoundaryReference?.ownerPartId),'A component cannot add another region or unrelated organ');
 }
});

test('a previously chosen hierarchy object cannot bring back a stale fragment list',()=>{
 const children=whole.elements.map(id=>{const p=byId.get(id);return {id:p.id,name:p.name,elements:[p.id]};});
 const {componentChoices}=inspect({...whole,id:'hierarchy:guest:ta98:A14.1.02.001',children});
 assert.equal(componentChoices.length,75);
});

test('regions, segments and individual nuclei remain bounded and independently selectable',()=>{
 const cervical=atlas.concepts.find(c=>c.ta98Term==='A14.1.02.014');
 const {componentChoices:segments}=inspect(cervical);
 assert.equal(segments.length,8);
 assert.deepEqual(new Set(segments.flatMap(c=>c.elements)),new Set(cervical.elements));
 for(const segment of registry.marker.segments){
  const c=atlas.concepts.find(c=>c.id===segment.conceptId);
  const {componentChoices}=inspect(c);
  assert.deepEqual(componentChoices.map(c=>c.name),[`White substance · ${segment.name}`,`Grey substance · ${segment.name}`,`Central canal · ${segment.name}`]);
  const labels=api.selectedAnatomyLabels(atlas,c.elements,c,true);
  assert.equal(labels.length,2,'Only the selected grey and white tissue controls need labels');
  for(const label of labels){
   const clicked=api.anatomyLabelChoice(label);
   assert.ok(clicked.elements.every(id=>segment.partIds.includes(id)));
   const nested=api.spinalCordComponentChoices(atlas,clicked);
   assert.ok(nested.length<60);
   assert.ok(nested.every(child=>child.elements.every(id=>segment.partIds.includes(id))));
  }
 }
 const nucleus=atlas.concepts.find(c=>c.id===registry.marker.termGroups.find(g=>g.term==='A14.1.02.121').conceptId);
 const {componentChoices}=inspect(nucleus);
 assert.equal(componentChoices.length,2);
 for(const child of componentChoices)assert.deepEqual(inspect(child).componentChoices,[],'Left and right named nuclei cannot recursively list themselves');
});

test('a single tissue such as pia has no self-repeating component branch',()=>{
 const pia=atlas.concepts.find(c=>c.ta98Term==='A14.1.01.309');
 assert.deepEqual(inspect(pia).componentChoices,[]);
});

test('covering relations merge fragment targets into independently selectable tissue names',()=>{
 const raw=registry.marker.partIds.map(id=>({kind:'covers',target:byId.get(id)}));
 const relations=api.spinalCordRelations(raw,atlas.concepts);
 assert.equal(relations.length,96);
 assert.deepEqual(new Set(relations.flatMap(r=>r.targetConcept.elements)),new Set(registry.marker.partIds));
 for(const relation of relations){
  assert.equal(relation.targetConcept.id,relation.target.conceptId);
  assert.ok(relation.targetConcept.elements.every(id=>byId.get(id).conceptId===relation.targetConcept.id));
 }
 assert.deepEqual(api.spinalCordRelations(raw.slice(0,2),atlas.concepts).flatMap(r=>r.targetConcept.elements),raw.slice(0,2).map(r=>r.target.id),'Partial relations must not select the rest of the cord');
});

test('other organs and the installed catalogue are unchanged by spinal presentation',()=>{
 const thyroid=atlas.concepts.find(c=>c.id===atlas.thyroidWallAssembly.wholeConcept);
 assert.equal(api.spinalCordComponentChoices(atlas,thyroid),undefined);
 assert.equal(api.spinalCordDisplayCount(thyroid.elements.map(id=>byId.get(id))),undefined);
 assert.equal(createHash('sha256').update(JSON.stringify(atlas)).digest('hex'),before);
});
