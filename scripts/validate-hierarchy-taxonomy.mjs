import assert from 'node:assert/strict';
import {readFileSync,existsSync} from 'node:fs';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/anatomy-hierarchy.ts','app/anatomical-terminology.ts','app/hierarchy-choice.ts','app/hierarchy-navigation.ts','app/depth-layers.ts','app/anatomy.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/hierarchy-taxonomy-validation'});
const load=async name=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(name)).text).toString('base64')}`);
const {MAJOR_SYSTEMS,REGION_ORDER,SKIN_DETAIL_NAMES,buildAnatomyNodes,depthPathFor,displaySystemFor,hierarchyEntries,majorSystemFor,regionPathFor,systemPathFor}=await load('anatomy-hierarchy.js');
const {DEPTH_LAYERS,depthLayerFor,isSkinPart}=await load('depth-layers.js');
const {SYSTEMS}=await load('anatomy.js');
const {anatomyNodeChoice}=await load('hierarchy-choice.js');
const {atlasIdentifier,terminologyTitle,terminologyForConcept,taHierarchyForConcept}=await load('anatomical-terminology.js');
const {hierarchyNavigation,hierarchyAncestors}=await load('hierarchy-navigation.js');
const ta98=JSON.parse(readFileSync('app/data/ta98-metadata.json','utf8'));
assert.equal(ta98.byConcept['LOCAL:female:cervix'].ta98,'A09.1.03.010','Pelvic cervix must not inherit dental-neck identity');
assert(ta98.byConceptPath['LOCAL:female:cervix'].includes('A09.1.03.001'),'Cervix must descend from uterus');
assert(!ta98.byConceptPath['LOCAL:female:cervix'].some(code=>code.startsWith('A05.')),'Pelvic cervix must not appear in digestive hierarchy');
for(const file of ['.local-models/female.json','.local-models/ta98-female-runtime.json']){
 if(!existsSync(file))continue;
 const atlas=JSON.parse(readFileSync(file,'utf8'));
 const cervix=hierarchyEntries(atlas).find(entry=>entry.id==='LOCAL:female:cervix');
 assert(cervix,`${file}: native cervix must be a selectable leaf`);
 assert.equal(majorSystemFor(cervix),'reproductive');
 assert.equal(systemPathFor(cervix)[0],'Female genital system');
 assert.equal(cervix.terminology.ta98,'A09.1.03.010');
}

// Palatine tonsil is cross-listed by TA98: oral identity, lymphoid system.
for(const file of ['.local-models/ta98-runtime.json','.local-models/ta98-female-runtime.json']){
 if(!existsSync(file))continue;
 const atlas=JSON.parse(readFileSync(file,'utf8'));
 const tonsils=hierarchyEntries(atlas).filter(entry=>entry.terminology.ta98==='A05.2.01.011');
 assert(tonsils.length>=2,`${file}: bilateral palatine tonsils must remain selectable`);
 for(const entry of tonsils){assert.equal(majorSystemFor(entry),'lymphatic');assert.equal(systemPathFor(entry)[0],'Lymphatic organs and vessels');}
}

assert.deepEqual(MAJOR_SYSTEMS.map(system=>system.name),[
 'Respiratory system','Digestive system','Circulatory system','Urinary system',
 'Integumentary system','Skeletal system','Muscular system','Endocrine system',
 'Exocrine system','Lymphatic system','Nervous system','Reproductive system',
]);
assert.deepEqual(ta98.byGroup.All,{ta98:'A01.0.00.000',tha:'THA:8719',fma:'FMA:20394',latin:'corpus humanum'});
assert.equal(ta98.byConcept['ZA:Stomach'].latin,'gaster');
assert.equal(ta98.byGroup['Pectoral region'].ta98,'A01.2.03.005');
assert.equal(ta98.byGroup['Pectoral region'].latin,'regio pectoralis');
assert.match(terminologyTitle('Stomach',ta98.byConcept['ZA:Stomach']),/TA98:A05\.\d\.\d{2}\.\d{3}/);
assert.match(terminologyTitle('Stomach',ta98.byConcept['ZA:Stomach']),/La:gaster/);
assert.equal(ta98.byConcept['HRA:VH_F_skin'].latin,'cutis');
assert.equal(ta98.byConcept['HRA:VH_F_skin'].ontology,'UBERON:0002097');
assert.equal(ta98.byConcept['HRA:VH_F_areolar_tubercles_L'].ontology,'UBERON:0011828');
assert.equal(ta98.byConcept['HRA:VH_F_areola_R'].fma,'FMA:223677');
assert.equal(ta98.byConcept['HRA:VH_F_amnion'].ontology,'UBERON:0000305');
assert.equal(ta98.byConcept['ZA:Abducens nerve (VI).l'].ta98,'A14.2.01.098');
assert.equal(ta98.byConcept['ZA:External abdominal oblique muscle.l'].ta98,'A04.5.01.008');
assert.equal(ta98.byConcept['ZA:External abdominal oblique muscle.ol'],undefined,'Attachment annotations must not inherit the muscle TA98 ID');
assert.equal(ta98.byConcept['ZA:Vestibulocochlear nerve (VIII).l'].ta98,'A14.2.01.121','Cranial nerve enumeration must survive bracket variants');
assert.equal(ta98.byConcept['ZA:Bucinator.l'].ta98,'A04.1.03.036','Reviewed source typo must resolve to buccinator');
assert.equal(ta98.byConcept['ZA:Sternocostal head of pectoralis major muscle.l'].ta98,'A04.4.01.004','Qualified head must use its own TA98 term');
assert.equal(ta98.byConcept['ZA:Abductor digiti minimi of hand.l'].ta98,'A04.6.02.062');
assert.equal(ta98.byConcept['ZA:Abductor digiti minimi of foot.l'].ta98,'A04.7.02.063','Homonymous hand and foot muscles must stay distinct');
assert.equal(ta98.byConceptMatch['ZA:Proximal phalanx of first finger of hand.l'].term,'A02.4.10.002');
assert.equal(ta98.byConceptMatch['ZA:Proximal phalanx of first finger of foot.l'].term,'A02.5.18.002');
assert(ta98.byConceptPath['ZA:Proximal phalanx of first finger of hand.l'].includes('A02.4.07.001'),'Hand phalanx must have the repaired hand-bone ancestry');
assert(!ta98.byConcept['ZA:Proximal phalanx of first finger of hand.l']?.ta98,'A numbered subdivision must not inherit the generic TA98 leaf ID');
assert(!ta98.byConcept['ZA:Proximal phalanx of first finger of hand.l']?.fma,'A numbered subdivision must not inherit its parent FMA ID');
assert.equal(ta98.byConceptMatch['O3M:1st lumbrical of hand.l'].kind,'parent');
assert.equal(ta98.byConcept['ZA:Left posterior lateral segment of liver (II)'].ta98,'A05.8.01.039');
assert.equal(ta98.byConcept['ZA:Posterior lateral segment of liver (VII)'].ta98,'A05.8.01.051','Left/right liver segment qualifiers must be preserved');
assert.equal(ta98.byConceptMatch['HRA:VH_F_suspensory_ligament_of_ovary_L'].term,'A09.1.01.018F');
assert(ta98.byConceptPath['HRA:VH_F_suspensory_ligament_of_ovary_L'].includes('A09.1.01.001'),'Gender-suffixed entities must retain their ovarian ancestry');
assert.equal(ta98.byConceptMatch['LOCAL:male:l_lumbrical1'],undefined,'Unqualified numbered lumbrical does not establish hand versus foot');
for(const term of Object.values(ta98.byConcept))if(term.ta98){
 assert.match(term.ta98,/^A\d{2}\.\d\.\d{2}\.\d{3}$/);
 assert.match(term.tha,/^THA:\d+$/);
 assert.ok(term.latin);
}
const canonical=new Set([
 ...REGION_ORDER,...SYSTEMS.map(system=>system.name),'Skin','Body surface','Named skin regions','Left arm','Right arm','Left leg','Right leg',
 'Hand & wrist','Forearm & elbow','Shoulder','Upper arm','Foot','Thigh & hip','Lower leg',
 'Head','Neck','Brain','Skull','Face & senses','Scalp & epicranium','Back & spine','Pelvis','Abdomen','Chest',
 'Mammary gland','Lacrimal glands','Subcutaneous tissue','Bones','Cartilage','Tendons','Muscles','Fascia','Heart','Arteries','Veins','Sensory organs',
 'Hair & nails','Joints and ligaments','Muscle attachments',
 'Digestive tract','Accessory digestive organs','Lungs','Airways','Kidneys','Urinary tract',
 'Lymph nodes','Lymphatic organs and vessels','Endocrine glands','Salivary glands','Skin glands',
 'Central nervous system','Peripheral nervous system','Placenta & pregnancy',
 'Female genital system','Male genital system',
 'Cell boundary','Nucleus','Cytoplasm',...SKIN_DETAIL_NAMES,
]);
function leaves(nodes){return nodes.flatMap(node=>node.kind==='group'?leaves(node.nodes):node.kind==='bilateral'?node.entries.flatMap(entry=>entry.parts):node.entry.parts);}
function inspect(nodes,knownParents){
 for(const node of nodes){
  if(node.kind!=='group')continue;
  assert(node.parts.length>0,`${node.name}: empty intermediate node`);
  assert(canonical.has(node.name)||knownParents.has(node.name)||node.terminology.ta98,`${node.name}: unresolvable middle node`);
  assert.deepEqual(node.parts.map(part=>part.id).sort(),leaves(node.nodes).map(part=>part.id).sort(),`${node.name}: parent and leaf membership differ`);
  const firstLeaf=node.nodes.findIndex(child=>child.kind==='entry'||child.kind==='bilateral');
  if(firstLeaf>=0)assert(node.nodes.slice(firstLeaf).every(child=>child.kind!=='group'),'Groups must precede leaves without an Other wrapper');
  const choice=anatomyNodeChoice(node,`test:${node.id}`);
  assert.deepEqual(choice.elements.sort(),node.parts.map(part=>part.id).sort(),`${node.name}: selectable membership differs`);
  assert.equal(choice.children.length,node.nodes.length,`${node.name}: details must list every direct child`);
  inspect(node.nodes,knownParents);
 }
}
function verify(atlas,name){
 const entries=hierarchyEntries(atlas),available=atlas.parts.filter(part=>!part.suppressed);
 const knownParents=new Set(entries.flatMap(entry=>entry.ancestry));
 const expected=available.map(part=>part.id).sort();
 assert.deepEqual(entries.flatMap(entry=>entry.parts).map(part=>part.id).sort(),expected,`${name}: entry coverage`);
 const modes={
  systems:()=>atlas.scope==='cell'
   ?[...new Set(entries.map(entry=>entry.system))].flatMap(system=>buildAnatomyNodes(entries.filter(entry=>entry.system===system),entry=>systemPathFor(entry,atlas.scope)))
   :MAJOR_SYSTEMS.flatMap(system=>buildAnatomyNodes(entries.filter(entry=>majorSystemFor(entry)===system.id),entry=>systemPathFor(entry,atlas.scope))),
  regions:()=>[...new Set(entries.map(entry=>entry.region))].flatMap(region=>buildAnatomyNodes(entries.filter(entry=>entry.region===region),regionPathFor)),
  depth:()=>DEPTH_LAYERS.flatMap(layer=>{
   const members=entries.flatMap(entry=>{const parts=entry.parts.filter(part=>depthLayerFor(part)===layer.id);return parts.length?[{...entry,parts}]:[];});
   return members.length?buildAnatomyNodes(members,entry=>depthPathFor(entry,layer.name)):[];
  }),
 };
 for(const [mode,make] of Object.entries(modes)){
  const nodes=make();inspect(nodes,knownParents);
  assert.deepEqual(leaves(nodes).map(part=>part.id).sort(),expected,`${name}: ${mode} has missing or duplicated items`);
  const root=hierarchyNavigation(atlas,mode),leafPaths=new Map();
  const walk=(node,path)=>{
   assert(atlasIdentifier(node.id),`${name}: ${mode} node ${node.name} has no ID`);
   if(!node.children){for(const id of node.elements)leafPaths.set(id,[...path,node]);return;}
   assert(node.children.length,`${name}: ${mode} has an empty group ${node.name}`);
   for(const child of node.children)walk(child,[...path,node]);
  };
  walk(root,[]);
  assert.deepEqual([...leafPaths.keys()].sort(),expected,`${name}: ${mode} navigation is missing leaves`);
  for(const id of expected.filter((_,index)=>index%Math.max(1,Math.floor(expected.length/50))===0)){
   const path=leafPaths.get(id),leaf=path.at(-1);
   assert.deepEqual(hierarchyAncestors(root,{id:leaf.id,name:leaf.name,elements:[id]}).map(item=>item.name),path.slice(0,-1).map(item=>item.name),`${name}: ${mode} breadcrumb differs for ${id}`);
  }
 }
 if(name==='male-detail'&&existsSync('public/assets/chakras.json')){
  const guest=JSON.parse(readFileSync('public/assets/chakras.json','utf8'));
  const root=hierarchyNavigation(atlas,`guest:${guest.id}`,guest),ids=new Set();
  const walk=node=>{assert(atlasIdentifier(node.id),`Chakras: ${node.name} has no ID`);if(node.children)node.children.forEach(walk);else node.elements.forEach(id=>ids.add(id));};
  walk(root);
  assert.deepEqual([...ids].sort(),expected,'Chakras hierarchy leaves do not cover the model');
 }
 if(name==='local-reference'){
  const skin=available.filter(part=>depthLayerFor(part)==='skin');
  assert.equal(skin.length,702,'All reference skin must be in the Skin depth layer');
  assert(skin.every(part=>displaySystemFor(part.system)==='integumentary'),'All reference skin must be under Body surface in Systems');
 }
 const expectedSkin={'male-detail':493,male:4,female:7,'local-male':2,'local-female':3,'local-reference':702};
 if(name in expectedSkin)assert.equal(available.filter(isSkinPart).length,expectedSkin[name],`${name}: skin classification changed`);
 for(const part of available.filter(part=>part.system==='integumentary'||part.system==='regions')){
  assert.equal(depthLayerFor(part)==='skin',isSkinPart(part),`${name}: ${part.name} has wrong skin depth`);
  const entry=entries.find(entry=>entry.parts.some(member=>member.id===part.id));
  assert(entry,`${name}: missing ${part.name}`);
  if(isSkinPart(part))assert(systemPathFor(entry,atlas.scope)[0]==='Skin',`${name}: ${part.name} is outside Skin`);
  else assert(systemPathFor(entry,atlas.scope)[0]!=='Skin',`${name}: ${part.name} is incorrectly under Skin`);
 }
 console.log(`${name}: ${available.length} pieces have one Systems, Regions, and Depth parent`);
}
const catalogues=[
 ['male-detail','public/models/atlas-male-complete.json'],['male','public/models/atlas.json'],
 ['female','public/models/atlas-hra-female.json'],['embryo','public/models/atlas-embryo.json'],
 ['cell','public/models/atlas-cell.json'],['local-male','.local-models/male.json'],
 ['local-female','.local-models/female.json'],['local-reference','.local-models/reference.json'],
 ['legacy-z-anatomy','public/models/atlas-z-anatomy.json'],
 ['local-ta98','.local-models/ta98-runtime.json'],['local-female-ta98','.local-models/ta98-female-runtime.json'],
];
const loaded=new Map();
for(const [name,file] of catalogues)if(existsSync(file)){
 const atlas=JSON.parse(readFileSync(file));loaded.set(name,atlas);verify(atlas,name);
}
function systemOf(model,label){
 const entry=hierarchyEntries(loaded.get(model)).find(entry=>entry.name.toLowerCase()===label.toLowerCase());
 assert(entry,`${model}: missing ${label}`);
 return majorSystemFor(entry);
}
assert.equal(systemOf('male','Third ventricle'),'nervous');
assert.equal(systemOf('male','Gingiva of upper jaw'),'digestive');
assert.equal(systemOf('male','Check ligament of left lateral rectus'),'nervous');
assert.equal(systemOf('female','Parotid gland (left)'),'exocrine');
assert.equal(systemOf('female','Mammary lobe'),'exocrine');
assert.equal(systemOf('female','Nipple (left)'),'integumentary');
assert.equal(systemOf('male-detail','Stomach'),'digestive');
const detailedEntries=hierarchyEntries(loaded.get('male-detail'));
const marker=detailedEntries.find(entry=>entry.id==='ZA:Bucinator.ol');
assert.equal(marker.taParentConcept,'ZA:Bucinator.l');
assert(marker.ancestry.includes('Facial muscles'),'Origin marker must follow the repaired parent muscle ancestry');
assert.equal(marker.terminology.ta98,null);
const phalanx=detailedEntries.find(entry=>entry.id==='ZA:Proximal phalanx of first finger of hand.l');
assert(systemPathFor(phalanx).includes('Bones of hand'));
assert(systemPathFor(phalanx).includes('Proximal phalanx'),'Parent-only placements must actually appear in Systems');
assert.equal(systemPathFor(detailedEntries.find(entry=>entry.id==='ZA:Vestibulocochlear nerve (VIII).l'))[0],'Peripheral nervous system','TA98 cranial nerves must override source sensory tags');
assert.equal(systemOf('male-detail','Left coronary artery'),'circulatory');
assert.equal(systemPathFor(hierarchyEntries(loaded.get('female')).find(entry=>entry.name==='Celiac trunk'))[0],'Arteries');
assert.equal(systemOf('local-reference','Body surface (derived)'),'integumentary');
if(loaded.has('male-detail')&&loaded.has('local-reference')){
 const detailed=loaded.get('male-detail'),reference=loaded.get('local-reference');
 const byName=new Map(hierarchyEntries(detailed).filter(entry=>systemPathFor(entry,detailed.scope)[0]==='Skin').map(entry=>[entry.name,entry]));
 let common=0;
 for(const entry of hierarchyEntries(reference).filter(entry=>systemPathFor(entry,reference.scope)[0]==='Skin')){
  const peer=byName.get(entry.name);if(!peer)continue;common++;
  assert.deepEqual(systemPathFor(entry,reference.scope),systemPathFor(peer,detailed.scope),`${entry.name}: system path differs between models`);
  assert.deepEqual(regionPathFor(entry),regionPathFor(peer),`${entry.name}: region path differs between models`);
 }
 assert(common>=250,`Only ${common} named skin structures share the canonical model paths`);
}

// Catalogue identity must override a contradictory global/source-name mapping.
const override={parts:[{id:'override-cervix',conceptId:'LOCAL:female:cervix',name:'Cervix',system:'reproductive',bounds:[[0,.9,0],[.02,.94,.02]]}],concepts:[{id:'LOCAL:female:cervix',name:'Cervix',elements:['override-cervix'],ta98Term:'A05.5.01.007',ta98Kind:'exact'}]};
const explicitEntry=hierarchyEntries(override)[0];
assert.equal(explicitEntry.terminology.ta98,'A05.5.01.007');
assert.equal(majorSystemFor(explicitEntry),'digestive');
assert(explicitEntry.ancestry.includes('Stomach'));
assert(!explicitEntry.ancestry.includes('Uterus'));
assert.equal(terminologyForConcept('LOCAL:female:cervix','A99.9.99.999','exact').ta98,null,'Unknown explicit identity must not resurrect shared fallback');
assert.deepEqual(taHierarchyForConcept('LOCAL:female:cervix','A99.9.99.999','exact'),[]);
assert.equal(terminologyForConcept('LOCAL:female:cervix','A05.5.01.007','parent').ta98,null,'Broader placement must not become an exact leaf badge');
assert(taHierarchyForConcept('LOCAL:female:cervix','A05.5.01.007','parent').some(parent=>parent.terminology.ta98==='A05.5.01.007'));
for(const name of ['local-ta98','local-female-ta98'])if(loaded.has(name)){
 const atlas=loaded.get(name),byId=new Map(atlas.concepts.map(c=>[c.id,c]));
 for(const entry of hierarchyEntries(atlas)){
  const source=byId.get(entry.id);
  if(!source?.ta98Term)continue;
  assert.equal(entry.ta98Term,source.ta98Term);
  if(source.ta98Kind==='parent')assert.equal(entry.terminology.ta98,null,`${name} ${entry.id}: placement identity falsely exact`);
  else assert.equal(entry.terminology.ta98,ta98.byCode[source.ta98Term]?.ta98??null,`${name} ${entry.id}: catalogue identity lost`);
  if(['A05.5.01.007','A05.5.01.009','A05.5.01.012','A05.5.01.014'].includes(source.ta98Term)){
   assert(entry.ancestry.includes('Stomach'),`${name} ${entry.id}: missing whole stomach ancestry`);
   assert.equal(majorSystemFor(entry),'digestive');
  }
 }
}
console.log('Explicit model identities, parent-placement guards and stomach ancestry verified for both runtime models.');

assert.equal(terminologyForConcept('HRA:VH_F_areola_R',ta98.byConceptMatch['HRA:VH_F_areola_R'].term,'exact').fma,'FMA:223677','Agreeing catalogue identities must preserve side-specific source FMA IDs');
assert.equal(terminologyForConcept('HRA:VH_F_aortic_arch','A12.2.04.001','exact').ontology,'UBERON:0001508','Agreeing source ontology must survive explicit catalogue identity');

// CS23 Coelom is an embryo cavity reference group, not a deep ligament.
if(existsSync('.local-models/embryo-cs23/atlas.json')){
 const embryo=JSON.parse(readFileSync('.local-models/embryo-cs23/atlas.json'));
 const cavities=hierarchyEntries(embryo).filter(e=>e.parts.some(p=>p.groups?.includes('Coelom')));
 assert.equal(cavities.length,3);
 for(const entry of cavities){
  assert.equal(depthLayerFor(entry.parts[0]),'anterior-organs');
  assert.deepEqual(depthPathFor(entry,'Anterior organs').slice(0,2),['Coelom','Torso & pelvis']);
 }
 console.log('CS23 Coelom: Anterior organs → Coelom → Torso & pelvis; all three cavities verified.');
}
