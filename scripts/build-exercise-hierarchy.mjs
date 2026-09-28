import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,existsSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';

export const EXERCISE_MAP=JSON.parse(readFileSync(new URL('./exercise-map.json',import.meta.url),'utf8'));
export const EXERCISE_ROOTS=[['strength','Strength'],['endurance','Aerobic endurance'],['mobility','Mobility and stretching'],['balance','Balance and coordination'],['breathing','Breathing'],['pelvic-floor','Pelvic floor'],['combined','Combined practices']];
const json=path=>JSON.parse(readFileSync(path,'utf8'));
const digest=content=>createHash('sha256').update(content).digest('hex');
const normalized=name=>name.toLowerCase().replace(/\((?:left|right)\)/g,'').replace(/\b(?:left|right)\b/g,'').replace(/\s+/g,' ').trim();
const sideOf=name=>/\bleft\b/i.test(name)?'left':/\bright\b/i.test(name)?'right':undefined;
const usableMuscle=part=>part.system==='muscular'&&!/\btendon\b|\.[eo]\d*[lr]$/i.test(part.name);
const unique=items=>[...new Set(items)].sort();
const slug=name=>name.toLowerCase().replace(/[^a-z0-9]+/g,'-').replace(/^-|-$/g,'');
const SOURCE_FILES=['muscle','exercisecategory','equipment','license','language','exerciseinfo'];
export function loadExerciseSources(){
 const manifest=json('scripts/exercise-sources.json');
 for(const source of manifest.files)assert.equal(digest(readFileSync(source.path)),source.sha256,`Source snapshot changed: ${source.path}`);
 const data=Object.fromEntries(SOURCE_FILES.map(key=>[key,json(`work/exercise/${key}.json`).results]));
 assert.equal(data.exerciseinfo.length,json('work/exercise/exerciseinfo.json').count,'Exercise source pagination is incomplete');
 return {manifest,data};
}

/** Human Atlas movement classification; original variation groups are retained. */
export function exerciseBranch(row,name){
 const n=name.toLowerCase();
 if(/kegel|pelvic (?:endurance|wave)/.test(n))return ['pelvic-floor',/reverse|relax/.test(n)?'Relaxation':'Contraction and relaxation'];
 if(/breathing|bobbing.*(?:exhale|recovery)|recovery bobbing/.test(n))return ['breathing',/diaphragm|belly|90\/90|wall baby/.test(n)?'Diaphragmatic breathing':'Breathing coordination'];
 if(/balance|wobble board|coordination/.test(n))return ['balance',/walk|coordination/.test(n)?'Dynamic balance':'Static balance'];
 if(/stretch|mobility|foam roll|foam roller|smr |cat.cow|cow.cat|child.?s pose|happy baby|pigeon|sphinx|hip circles|ankle roll|wrist circles|neck circle|neck cars|shoulder rotation|head turn|head tilt|chin tuck|arm circle|open book|hip switch|thoracic rotation|adductor rock back|pendular|shoulder dislocate|limber/.test(n)){
  const area=/neck|chin|head/.test(n)?'Neck':/ankle|calf|soleus|achilles|foot|tibialis/.test(n)?'Ankle and foot':/hip|glute|quad|hamstring|adductor|lunge|pigeon|split|clamshell|frog|pancake/.test(n)?'Hip and knee':/shoulder|arm|scapula|pectoral|lat |latissimus|sleeper|prayer/.test(n)?'Shoulder':/wrist|forearm/.test(n)?'Wrist and hand':'Spine and trunk';return ['mobility',area];
 }
 if(row.category?.name==='Cardio'||/cycling|stationary bike|treadmill|jogging|^run\b|running|swim|elliptical|jump rope|skipping|rowing erg|rower/.test(n)){
  const kind=/swim|kickboard|kick with board/.test(n)?'Swimming':/cycling|bike|bicycl/.test(n)?'Cycling':/walk|hike/.test(n)?'Walking':/run|jog/.test(n)?'Running':/elliptical/.test(n)?'Elliptical':/rope|skipping|jump/.test(n)?'Jumping':/row/.test(n)?'Rowing':'Other conditioning';return ['endurance',kind];
 }
 const patterns=[['Squat',/squat|v.squat/],['Lunge and step-up',/lunge|step.up/],['Hip hinge',/deadlift|\brdl\b|good morning|kettlebell swing|clean|snatch/],['Hip extension',/hip (?:thrust|raise|bridge)|glute (?:bridge|kick)|donkey|kickback.*(?:glute|hip)|kneeling kickback/],['Hip abduction',/hip abduction|clamshell|monster walk|fire hydrant/],['Hip adduction',/hip adduction|standing adduction|adductor/],['Knee flexion',/(?:leg|hamstring|nordic) curl|glute.ham raise/],['Knee extension',/leg extension|knee extension|quad sets/],['Ankle plantar flexion',/calf|heel raise/],['Horizontal push',/bench.?press|bench press|chest press|floor press|push.up|pushup|press.up|pec deck/],['Chest adduction',/fly|butterfly|cross.over|crossess/],['Vertical push',/shoulder press|overhead press|military press|arnold|landmine press|push press|thruster|handstand push|pike push/],['Vertical pull',/pull.up|pullup|chin.up|chin up|pulldown|pull.down|lat pull|muscle up/],['Horizontal pull',/row/],['Shoulder rotation',/external rotation|internal rotation|rotator cuff/],['Shoulder elevation',/lateral raise|side lateral|front raise|shoulder raise|delt raise|upright|high pull/],['Scapular control',/shrug|scapul|face.?pull|ytw|ywt|wall angel|pull.apart/],['Elbow flexion',/curl|bicep/],['Elbow extension',/tricep|skull.?crush|skullcrusher|dip|tate press/],['Wrist and grip',/wrist|forearm|grip|fingerboard|sloper|flexbar/],['Trunk flexion',/crunch|sit.up|situp|leg raise|knee raise|knee tuck|rollout|ab wheel/],['Trunk rotation',/twist|rotation|wood.?chop|windshield/],['Trunk stabilization',/plank|hollow|l.sit|l hold|abdominal|brac|deadbug|dead bug|bird dog/],['Back extension',/hyperextension|superman|back extension/],['Carry and hold',/carry|hold|hang|lever|planche/],['Combined movements',/burpee|bear|mountain climber|turkish get/]];
 const found=patterns.find(([,pattern])=>pattern.test(n));return ['strength',found?.[0]??`${row.category?.name??'Other'} movements`];
}

export async function buildExerciseHierarchy(input,source=loadExerciseSources()){
 const {manifest,data}=source,nodes=new Map(),add=node=>{assert(!nodes.has(node.id),`Duplicate ${node.id}`);nodes.set(node.id,node);return node;},node=id=>nodes.get(id);
 const roots=EXERCISE_ROOTS.map(([id,name])=>add({id:`TYPE:${id}`,name,childIds:[]}).id);
 const references={},relations=new Map(),anatomyEvidence={};
 const sideIds={left:new Set(),right:new Set()};
 for(const {atlas} of input.catalogues){const names=new Map(atlas.concepts.map(c=>[c.id,c.name]));for(const part of atlas.parts){const side=sideOf(names.get(part.conceptId)??part.name);for(const requested of ['left','right'])if(!side||side===requested)sideIds[requested].add(part.id);}}
 for(const [key,spec] of Object.entries(EXERCISE_MAP.anatomy)){
  const pattern=new RegExp(spec.pattern,'i'),fallback=spec.fallback?new RegExp(spec.fallback,'i'):null;
  const matches=[],bySide={left:[],right:[]},extras={nerves:[],tendons:[],bones:[]},levels=new Set();
  for(const {file,atlas} of input.catalogues){
   const conceptNames=new Map(atlas.concepts.map(c=>[c.id,c.name]));
   let context=relations.get(atlas);if(!context){const parts=atlas.parts.filter(p=>!p.suppressed).map(p=>({...p,name:conceptNames.get(p.conceptId)??p.name}));context={parts,cache:new Map()};relations.set(atlas,context);}
   const eligible=context.parts.filter(part=>spec.systems.includes(part.system)&&(part.system!=='muscular'||key==='achilles'||usableMuscle(part)));
   let found=eligible.filter(part=>pattern.test(normalized(part.name)));
   if(!found.length&&fallback)found=eligible.filter(part=>fallback.test(normalized(part.name)));
   const relevant=new Set(found.map(p=>p.conceptId));anatomyEvidence[`${key}:${file}`]=unique(found.map(p=>p.id));
   matches.push(...found.map(p=>p.id));
   for(const side of ['left','right'])bySide[side].push(...found.filter(p=>!sideOf(p.name)||sideOf(p.name)===side).map(p=>p.id));
   for(const id of relevant){
    const part=found.find(p=>p.conceptId===id);if(!usableMuscle(part))continue;
    let links=context.cache.get(id);if(!links){links=input.anatomicalRelations(part,context.parts);context.cache.set(id,links);}
    for(const link of links){
     if(link.kind==='innervation'&&link.target.system==='nervous')extras.nerves.push(link.target.id);
     if(link.kind==='tendons'&&/tendon|aponeurosis/i.test(link.target.name))extras.tendons.push(link.target.id);
     if(['origin','insertion','bones'].includes(link.kind)&&link.target.system==='skeletal')extras.bones.push(link.target.id);
    }
   }
  }
  references[key]={matches:unique(matches),bySide:Object.fromEntries(Object.entries(bySide).map(([side,ids])=>[side,unique(ids)])),extras:Object.fromEntries(Object.entries(extras).map(([role,ids])=>[role,unique(ids)]))};
 }
 const anatomy=(key,side)=>{
  const id=`ANATOMY:${key}${side?':'+side:''}`;if(nodes.has(id))return id;
  const ref=references[key],spec=EXERCISE_MAP.anatomy[key];assert(spec,`Unknown structure ${key}`);
  add({id,name:spec.name+(side?` (${side})`:''),partIds:side?ref.bySide[side]:ref.matches});return id;
 };
 const anatomyRole=(parent,role,name,keys,side)=>{
  if(!keys?.length)return;
  const id=`ROLE:${parent.id}:${role}`;add({id,name,childIds:unique(keys).map(key=>anatomy(key,side))});parent.childIds.push(id);
 };
 const partRole=(parent,role,name,ids)=>{
  if(!ids.length)return;
  const id=`ROLE:${parent.id}:${role}`;add({id,name,partIds:unique(ids)});parent.childIds.push(id);
 };
 const families=new Map();
 const family=(root,name)=>{const id=`FAMILY:${root}:${slug(name)}`;if(!families.has(id)){families.set(id,add({id,name,childIds:[]}));node(`TYPE:${root}`).childIds.push(id);}return families.get(id);};
 const provenance={sources:EXERCISE_MAP.sources,snapshot:manifest,exercises:[],anatomySelectors:anatomyEvidence,reviewRules:EXERCISE_MAP.reviews,statistics:{}};
 const groups=new Map(),sourceRows=[];
 for(const row of data.exerciseinfo){const en=row.translations.filter(t=>t.language===2).sort((a,b)=>a.id-b.id)[0];if(en)sourceRows.push({row,en});}
 for(const {row,en} of sourceRows){const [root,branch]=exerciseBranch(row,en.name),key=`${root}:${branch}:${row.variation_group??row.id}`;groups.set(key,[...(groups.get(key)??[]),{row,en}]);}
 const attach=(parent,record)=>{
  const primary=unique(record.primary??[]),secondary=unique(record.secondary??[]).filter(key=>!primary.includes(key));
  const isStretch=record.root==='mobility',side=record.side;
  anatomyRole(parent,'primary',isStretch?'Muscles involved':record.root==='pelvic-floor'?'Pelvic floor muscles':'Primary muscles',primary,side);
  anatomyRole(parent,'secondary','Secondary muscles',secondary,side);
  anatomyRole(parent,'stabilizers','Stabilizing muscles',record.stabilizers,side);
  anatomyRole(parent,'joints','Joints involved',record.joints,side);
  anatomyRole(parent,'tissues','Tendons and fascia',record.tendons,side);
  anatomyRole(parent,'organs',record.root==='pelvic-floor'?'Supported organs':'Respiratory involvement',record.organs);
  if(record.aerobic)anatomyRole(parent,'cardiopulmonary','Cardiopulmonary involvement',['heart','lungs']);
  if(record.balance)anatomyRole(parent,'balance','Vestibular involvement',['vestibular']);
  const keys=[...primary,...secondary,...(record.stabilizers??[])];
  for(const [role,label] of [['nerves','Motor innervation'],['tendons','Related tendons'],['bones','Attachment bones']]){
   const ids=keys.flatMap(key=>references[key].extras[role]);
   if(role==='nerves'&&primary.includes('diaphragm'))ids.push(...references.phrenic.matches);
   // Preserve explicit laterality when an exercise is left/right specific.
   const allowed=side?sideIds[side]:null;
   partRole(parent,role,label,allowed?ids.filter(id=>allowed.has(id)):ids);
  }
  if(primary.includes('diaphragm')){
   const ref=node(anatomy('diaphragm'));ref.links=['C3','C4','C5'].map(level=>({hierarchy:'dermatomes-myotomes',id:`ROOT:${level}`,name:level}));
  }
 };
 for(const [group,rows] of groups){
  const {row:first,en:firstName}=rows[0],[root,branch]=exerciseBranch(first,firstName.name),parent=family(root,branch);
  const common=rows.length>1?add({id:`VARIATIONS:${root}:${slug(branch)}:${first.variation_group}`,name:[...rows].map(({en})=>en.name).sort((a,b)=>a.length-b.length||a.localeCompare(b))[0]+' · variations',childIds:[]}):null;
  if(common)parent.childIds.push(common.id);
  for(const {row,en} of rows){
   const review=EXERCISE_MAP.reviews.find(rule=>new RegExp(rule.pattern,'i').test(en.name));
   const primary=review?.primary??row.muscles.map(m=>EXERCISE_MAP.wgerMuscles[m.id]);
   const secondary=review?.secondary??row.muscles_secondary.map(m=>EXERCISE_MAP.wgerMuscles[m.id]);
   assert(primary.every(Boolean)&&secondary.every(Boolean),'Unknown source muscle ID');
   const record={id:`WGER:${row.id}`,name:en.name,root,family:branch,primary,secondary,source:review?.source??'wger',review:review?.id,side:/\b(left|right)\s*$/i.exec(en.name)?.[1]?.toLowerCase(),organs:review?.organs,stabilizers:review?.stabilizers,joints:review?.joints,additionalSources:review?.additionalSources,aerobic:root==='endurance'&&!/drill|recovery|stretch|cool.down|bobbing/i.test(en.name),balance:root==='balance'};
   const synonyms=unique([...en.aliases.map(a=>a.alias),...row.equipment.map(e=>e.name),...row.muscles.map(m=>m.name_en)].filter(Boolean)).filter(s=>s.length<=200&&s.toLowerCase()!==en.name.toLowerCase()).slice(0,300);
   const item=add({id:record.id,name:record.name,...(synonyms.length?{synonyms}:{}),childIds:[]});(common??parent).childIds.push(item.id);attach(item,record);
   provenance.exercises.push({...record,original:{primary:row.muscles.map(m=>m.id),secondary:row.muscles_secondary.map(m=>m.id),category:row.category,equipment:row.equipment,variationGroup:row.variation_group},credit:{authors:unique([row.license_author,en.license_author,...(row.total_authors_history??[]),...(en.author_history??[])].filter(Boolean)),license:data.license.find(l=>l.id===en.license),source:`https://wger.de/api/v2/exerciseinfo/${row.id}/`}});
  }
 }
 for(const record of EXERCISE_MAP.curated){const item=add({id:record.id,name:record.name,childIds:[]});family(record.root,record.family).childIds.push(item.id);attach(item,record);provenance.exercises.push(record);}
 for(const practice of EXERCISE_MAP.combined){
  const matches=provenance.exercises.filter(row=>practice.patterns.some(pattern=>new RegExp(pattern,'i').test(row.name))).map(row=>row.id);
  const item=add({id:`PRACTICE:${practice.id}`,name:practice.name,childIds:unique(matches)});node('TYPE:combined').childIds.push(item.id);
 }
 for(const item of nodes.values())if(item.childIds){item.childIds=[...new Set(item.childIds)];if(!item.childIds.length)delete item.childIds;}
 // Keep the menu order and original source exercise order stable, sort only groups.
 for(const id of roots)if(node(id).childIds)node(id).childIds.sort((a,b)=>node(a).name.localeCompare(node(b).name));
 const hierarchy={schema:'human-atlas-hierarchy/v2',id:'physical-exercise',name:'Physical Exercise',version:`wger snapshot ${manifest.date} · anatomy map ${EXERCISE_MAP.version}`,roots,nodes:[...nodes.values()]};input.parseGuestHierarchy(hierarchy);
 provenance.statistics={sourceExercises:sourceRows.length,additionalExercises:EXERCISE_MAP.curated.length,reviewedSourceExercises:provenance.exercises.filter(e=>e.review).length,anatomicalGroups:Object.keys(references).length,nodes:hierarchy.nodes.length};
 return {hierarchy,provenance};
}

function report(hierarchy,provenance,input){
 const view=id=>'http://localhost:3016/?'+new URLSearchParams({model:'male-detail',tree:'guest:physical-exercise',bio:id,context:'0.18'});
 const lines=['# Physical Exercise hierarchy','','Physical Exercise appears immediately before Chakras, after Drugs, in Additional hierarchies. It loads only when selected.','','## Scope','',`Pinned wger snapshot: ${provenance.snapshot.date}. ${provenance.statistics.sourceExercises} English exercise records and ${provenance.statistics.additionalExercises} additional reference exercises. ${provenance.statistics.reviewedSourceExercises} source exercises have corrections or extra mappings against the named clinical references.`,'','The classification is an original Human Atlas grouping by exercise type and movement family. Source variation groups are preserved within each family. This is a catalogue of available exercises, not an exhaustive taxonomy of every human movement.','','## Types','','| Type | Open in viewer |','| --- | --- |'];
 for(const [id,name] of EXERCISE_ROOTS)lines.push(`| ${name} | [Open](${view('TYPE:'+id)}) |`);
 lines.push('','## Anatomical roles','','- Primary and secondary muscle associations retain the wger source roles, with the reviewed corrections recorded in the provenance. Secondary does not imply a quantified stabilizing role. Mobility nodes say “Muscles involved” rather than treating stretches as prime movers.','- Quadriceps expands to its four muscles. A source association with biceps femoris does not automatically include every hamstring. Anterior, middle and posterior deltoid portions remain separate where source geometry permits; models lacking separated portions use their whole-deltoid reference mesh.','- Named joint, tendon, plantar fascia, respiratory, supported pelvic organ and vestibular associations are independently sourced. Supported organs and physiological involvement do not indicate direct exercise targets.','- Motor innervation, attachment bones and related tendons reuse existing anatomical relationships of the mapped muscles. They are anatomical context, not a measured exercise activation result. Breathing links the phrenic nerves and the C3–C5 root entries.','- Exact selectors contain only existing, nonsuppressed meshes from the source catalogues. Muscle roles exclude separately labeled tendons and origin/insertion markers. Left/right variants select only that side plus source meshes with no side designation.','- Search includes exercise names, source IDs, aliases, equipment and muscle group names. Share links, slides, connections and embeds use `guest:physical-exercise`.','- Source instructions, training prescriptions, media, citations and research notes are not displayed in the hierarchy.','','## Geometry coverage','','| Model catalogue | Exercises with anatomy | Total exercise records |','| --- | ---: | ---: |');
 const mappedAnywhere=new Set();
 for(const {file,atlas} of input.catalogues){const resolved=new Map(),visit=n=>{if(resolved.has(n.id))return;resolved.set(n.id,n);n.children.forEach(visit);};input.resolveGuestHierarchy(atlas,hierarchy).forEach(visit);const mapped=provenance.exercises.filter(e=>resolved.get(e.id)?.parts.length).length;lines.push(`| ${file} | ${mapped} | ${provenance.exercises.length} |`);for(const e of provenance.exercises)if(resolved.get(e.id)?.parts.length)mappedAnywhere.add(e.id);}
 lines.push('','## Associations still needing review','','The following exercises have no mapped anatomy in any available catalogue. Links open their hierarchy entries; absent source associations are kept visible instead of filled with guessed structures. Geometry coverage differs by model as shown above; the placenta-only Embryo model has no exercise anatomy.','');
 for(const e of provenance.exercises.filter(e=>!mappedAnywhere.has(e.id)).sort((a,b)=>a.name.localeCompare(b.name)))lines.push(`- [${e.name.replace(/[[\]|]/g,'')}](${view(e.id)})`);
 for(const p of EXERCISE_MAP.combined)if(!hierarchy.nodes.find(n=>n.id==='PRACTICE:'+p.id)?.childIds?.length)lines.push(`- [${p.name}](${view('PRACTICE:'+p.id)}) — no individually mapped movements in the imported snapshot.`);
 lines.push('','## Sources and rebuilding','','All source credits and per-exercise authors/licenses are in [ATTRIBUTION.md](ATTRIBUTION.md#physical-exercise). Exact source hashes are in [scripts/exercise-sources.json](scripts/exercise-sources.json); original roles, corrections and geometry evidence are in [public/assets/exercise-provenance.json](public/assets/exercise-provenance.json).','','```sh','python3 scripts/download-exercise-sources.py  # fetch missing pinned snapshots','npm run build:exercise','npm run build:exercise -- --check','npm run test:exercise','npm run check','npm run build','```','');return lines.join('\n');
}
function credits(provenance){
 const lines=['<!-- exercise-credits:start -->','### Physical Exercise','','Exercise records and original primary/secondary muscle associations: [wger](https://wger.de/), snapshot '+provenance.snapshot.date+'. [Public API documentation](https://wger.readthedocs.io/en/stable/api/api.html). Each imported record retains its authors and license below. Derived exercise data preserves the applicable CC BY-SA 3.0 / 4.0 or CC0 terms independently of the application code license. Human Atlas adds the movement classification, reviewed corrections, exact anatomy selectors and role branches. No exercise descriptions, images or videos are imported.','','Classification and factual anatomical associations:',''];
 for(const [key,source] of Object.entries(EXERCISE_MAP.sources))if(key!=='wger')lines.push(`- [${source.name}](${source.url})${source.use?' — '+source.use:'. Original paraphrases and factual anatomical associations; instructions and illustrations are not reproduced.'}`);
 lines.push('- Existing muscle innervation, tendon and bone attachment references: the anatomical relation sources credited in this document.','','Pinned API snapshots:','');
 for(const file of provenance.snapshot.files)lines.push(`- [${file.key}](${file.url}) · SHA-256 \`${file.sha256}\``);
 lines.push('','#### wger exercise authors and licenses','','| Exercise | Authors retained from source history | License |','| --- | --- | --- |');
 const clean=s=>s.replace(/[|\r\n]/g,' ').replace(/[[\]]/g,'');
 for(const e of provenance.exercises.filter(e=>e.credit).sort((a,b)=>a.id.localeCompare(b.id,undefined,{numeric:true})))lines.push(`| [${clean(e.name)}](${e.credit.source}) | ${e.credit.authors.map(clean).join(', ')||'wger contributors'} | [${e.credit.license.short_name}](${e.credit.license.url}) |`);
 lines.push('','<!-- exercise-credits:end -->');return lines.join('\n');
}
export async function writeExerciseOutputs(result,input,check=false){
 const root=readFileSync('ATTRIBUTION.md','utf8'),block=credits(result.provenance);
 const attribution=root.includes('<!-- exercise-credits:start -->')?root.replace(/<!-- exercise-credits:start -->[\s\S]*?<!-- exercise-credits:end -->/,block):root.replace('### Chakras',block+'\n\n### Chakras');
 const publicAttribution=attribution.replaceAll('(LICENSE)','(LICENSE.txt)').replaceAll('(LICENSES/MIT-upstream.txt)','(MIT-upstream.txt)');
 const outputs=[['public/assets/physical-exercise.json',JSON.stringify(result.hierarchy)+'\n'],['public/assets/exercise-provenance.json',JSON.stringify(result.provenance)+'\n'],['PHYSICAL-EXERCISE-HIERARCHY.md',report(result.hierarchy,result.provenance,input)],['ATTRIBUTION.md',attribution],['public/ATTRIBUTION.md',publicAttribution]];
 for(const [path,content] of outputs)if(check)assert.equal(readFileSync(path,'utf8'),content,`${path} needs rebuilding`);else writeFileSync(path,content);
}
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){const input=await loadChakraInputs(),result=await buildExerciseHierarchy(input);await writeExerciseOutputs(result,input,process.argv.includes('--check'));console.log(result.provenance.statistics);}
