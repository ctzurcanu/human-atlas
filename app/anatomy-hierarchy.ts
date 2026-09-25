import {SYSTEMS,structureName,type Atlas,type Part,type SystemId} from './anatomy';

export type AnatomyEntry={id:string;name:string;system:SystemId;region:string;location:string[];parts:Part[]};
export type AnatomyNode=
 | {kind:'group';id:string;name:string;nodes:AnatomyNode[];parts:Part[]}
 | {kind:'bilateral';id:string;name:string;entries:AnatomyEntry[];parts:Part[]}
 | {kind:'entry';id:string;name:string;entry:AnatomyEntry;parts:Part[]};
export const REGION_ORDER=['Head & neck','Torso & pelvis','Upper limbs','Lower limbs','Whole body & spanning'];
const CELL_COMPARTMENTS=['Cell boundary','Nucleus','Cytoplasm'];
const normalized=(parts:Part[])=>structureName(parts[0].name).toLowerCase();
const groupsFor=(parts:Part[])=>new Set(parts.flatMap(part=>part.groups??[]).map(group=>group.toLowerCase()));
const has=(groups:Set<string>,pattern:RegExp)=>[...groups].some(group=>pattern.test(group));
const boundsFor=(parts:Part[])=>({
 minX:Math.min(...parts.map(part=>part.bounds[0][0])),
 maxX:Math.max(...parts.map(part=>part.bounds[1][0])),
 minY:Math.min(...parts.map(part=>part.bounds[0][1])),
 maxY:Math.max(...parts.map(part=>part.bounds[1][1])),
});

function regionFor(parts:Part[],scope?:string){
 if(scope==='cell')return parts[0].groups?.find(group=>CELL_COMPARTMENTS.includes(group))??'Cytoplasm';
 const {minX,maxX,minY,maxY}=boundsFor(parts);
 if(maxY-minY>1.05)return 'Whole body & spanning';
 const groups=groupsFor(parts);
 if(has(groups,/^(head|neck|brain|face|skull|cranium|cerebrum|cerebellum|left head|right head|head and neck)$/))return 'Head & neck';
 if(has(groups,/^(left |right )?(upper limb|arm|forearm|hand|wrist|shoulder)( region)?$/))return 'Upper limbs';
 if(has(groups,/^(left |right )?(lower limb|leg|thigh|foot|ankle|talocrural region)( region)?$/))return 'Lower limbs';
 if(has(groups,/^(trunk|torso|thorax|chest|abdomen|pelvis|pelvic region)$/))return 'Torso & pelvis';
 const name=normalized(parts);
 if(/^(skin of body|whole body|body surface)$/.test(name))return 'Whole body & spanning';
 const x=(minX+maxX)/2,y=(minY+maxY)/2;
 return y<.73?'Lower limbs':Math.abs(x)>.19&&y>.7&&y<1.55?'Upper limbs':y>1.4?'Head & neck':'Torso & pelvis';
}

function sideFor(parts:Part[],groups:Set<string>){
 const name=normalized(parts);
 const explicit=/\b(left|right)\b/.exec(name)?.[1];
 if(explicit)return explicit==='left'?'Left':'Right';
 const left=has(groups,/^left (?:upper limb|lower limb|arm|hand|foot|leg|head)$/);
 const right=has(groups,/^right (?:upper limb|lower limb|arm|hand|foot|leg|head)$/);
 if(left!==right)return left?'Left':'Right';
 const {minX,maxX}=boundsFor(parts);
 return minX>.025?'Left':maxX<-.025?'Right':'';
}

function locationFor(parts:Part[],region:string):string[]{
 const groups=groupsFor(parts),name=normalized(parts),side=sideFor(parts,groups);
 const {minY,maxY}=boundsFor(parts),y=(minY+maxY)/2;
 if(region==='Lower limbs'){
  const area=/\b(foot|toe|tars|metatars|calcane|talus|ankle|malleol|plantar|dorsum of foot)\b/.test(name)||has(groups,/^(left |right )?foot$/)?'Foot'
   :/\b(thigh|femur|patella|knee|inguinal|hip|glute|buttock)\b/.test(name)||has(groups,/^thigh$/)?'Thigh & hip'
   :/\b(leg|tibia|fibula|calf|popliteal)\b/.test(name)?'Lower leg':y<.22?'Foot':y<.44?'Lower leg':'Thigh & hip';
  return [side?`${side} leg`:'Central lower limbs',area];
 }
 if(region==='Upper limbs'){
  const area=/\b(hand|finger|thumb|digit|metacarp|carpal|wrist|palm)\b/.test(name)||has(groups,/^(left |right )?hand$/)?'Hand & wrist'
   :/\b(forearm|radius|ulna|radial|ulnar|elbow)\b/.test(name)?'Forearm & elbow'
   :/\b(shoulder|scapula|clavicle|pectoral girdle|axilla)\b/.test(name)?'Shoulder'
   :/\b(arm|humerus|brachial|biceps|triceps)\b/.test(name)?'Upper arm':y<.94?'Hand & wrist':y<1.12?'Forearm & elbow':y<1.4?'Upper arm':'Shoulder';
  return [side?`${side} arm`:'Central upper limbs',area];
 }
 if(region==='Head & neck'){
  if(/\b(neck|cervical|larynx|pharynx|thyroid|trachea)\b/.test(name)||has(groups,/^neck$/))return ['Neck'];
  const area=/\b(brain|cerebr|cerebell|telencephal|ventric|thalam|hypothalam|pons|medulla|pineal|pituitary)\b/.test(name)||has(groups,/^(brain|cerebrum|telencephalon)$/)?'Brain'
   :/\b(skull|crani|parietal|occipital|frontal bone|temporal bone|sphenoid|ethmoid|mandib|maxill|zygomat|vomer)\b/.test(name)||has(groups,/^cranium$/)?'Skull'
   :/\b(eye|ear|nose|nasal|mouth|tongue|tooth|teeth|facial|face|orbit|cheek|lip|chin)\b/.test(name)||has(groups,/^face$/)?'Face & senses':'Other head';
  return ['Head',area];
 }
 if(region==='Torso & pelvis'){
  const area=/\b(back|vertebr|spine|spinal|lumbar|dorsal|erector spinae)\b/.test(name)||has(groups,/^back$/)?'Back & spine'
   :/\b(pelvi|sacrum|coccyx|bladder|uterus|ovary|prostate|rectum|perine|genital|urethra|vagina|testis|penis)\b/.test(name)||has(groups,/^(pelvis|pelvic region)$/)?'Pelvis'
   :/\b(abdom|kidney|renal|liver|stomach|intestin|colon|spleen|pancreas|duodenum|jejunum|ileum)\b/.test(name)||has(groups,/^abdomen$/)?'Abdomen':'Chest';
  return [area];
 }
 return [];
}

export function hierarchyEntries(atlas:Atlas):AnatomyEntry[]{
 const names=new Map(atlas.concepts.map(concept=>[concept.id,concept.name]));
 const byConcept=new Map<string,Part[]>();
 for(const part of atlas.parts){const parts=byConcept.get(part.conceptId)??[];parts.push(part);byConcept.set(part.conceptId,parts);}
 return [...byConcept].map(([id,parts])=>{
  const region=regionFor(parts,atlas.scope);
  return {id,name:names.get(id)??structureName(parts[0].name),system:parts[0].system,region,location:atlas.scope==='cell'?[]:locationFor(parts,region),parts};
 });
}

export function entryLabel(entry:AnatomyEntry){
 const side=/^(left|right) deferent duct$/i.exec(entry.name);
 const name=side?`Ductus deferens (${side[1].toLowerCase()})`:entry.name;
 return name.replace(/^./,letter=>letter.toUpperCase());
}
function bilateral(name:string){
 const suffix=/^(.*) \((left|right)\)$/i.exec(name);
 if(suffix)return {base:suffix[1],side:suffix[2].toLowerCase()};
 const prefix=/^(left|right) (.+)$/i.exec(name);
 return prefix?{base:prefix[2].replace(/^./,letter=>letter.toUpperCase()),side:prefix[1].toLowerCase()}:null;
}
function leafNodes(entries:AnatomyEntry[],compare?:(a:AnatomyEntry,b:AnatomyEntry)=>number):AnatomyNode[]{
 const pairs=new Map<string,AnatomyEntry[]>();
 for(const entry of entries){const parsed=bilateral(entryLabel(entry));if(parsed){const list=pairs.get(parsed.base)??[];list.push(entry);pairs.set(parsed.base,list);}}
 const grouped=new Set<string>(),nodes:AnatomyNode[]=[];
 for(const entry of entries){
  if(grouped.has(entry.id))continue;
  const parsed=bilateral(entryLabel(entry)),siblings=parsed?pairs.get(parsed.base)??[]:[];
  if(parsed&&siblings.some(sibling=>bilateral(entryLabel(sibling))?.side==='left')&&siblings.some(sibling=>bilateral(entryLabel(sibling))?.side==='right')){
   siblings.forEach(sibling=>grouped.add(sibling.id));
   nodes.push({kind:'bilateral',id:`pair:${parsed.base}`,name:parsed.base,entries:siblings.sort((a,b)=>entryLabel(a).localeCompare(entryLabel(b))),parts:siblings.flatMap(sibling=>sibling.parts)});
  }else nodes.push({kind:'entry',id:`entry:${entry.id}`,name:entryLabel(entry),entry,parts:entry.parts});
 }
 return nodes.sort((a,b)=>Number(b.kind==='bilateral')-Number(a.kind==='bilateral')||(compare&&a.kind==='entry'&&b.kind==='entry'?compare(a.entry,b.entry):0)||a.name.localeCompare(b.name,undefined,{numeric:true,sensitivity:'base'}));
}

/** Each entry follows exactly one anatomical path; branches are added only when they contain multiple entries. */
export function buildAnatomyNodes(entries:AnatomyEntry[],path:(entry:AnatomyEntry)=>string[],compare?:(a:AnatomyEntry,b:AnatomyEntry)=>number):AnatomyNode[]{
 const descend=(items:AnatomyEntry[],depth:number):AnatomyNode[]=>{
  const byGroup=new Map<string,AnatomyEntry[]>(),leaves:AnatomyEntry[]=[];
  for(const entry of items){
   const name=path(entry)[depth];
   if(!name){leaves.push(entry);continue;}
   const group=byGroup.get(name)??[];group.push(entry);byGroup.set(name,group);
  }
  const groups:AnatomyNode[]=[];
  for(const [name,members] of byGroup){
   if(members.length===1){leaves.push(members[0]);continue;}
   groups.push({kind:'group',id:`group:${depth}:${name}`,name,nodes:descend(members,depth+1),parts:members.flatMap(entry=>entry.parts)});
  }
  groups.sort((a,b)=>{
   const ai=REGION_ORDER.indexOf(a.name),bi=REGION_ORDER.indexOf(b.name);
   return ai>=0&&bi>=0?ai-bi:a.name.localeCompare(b.name,undefined,{numeric:true,sensitivity:'base'});
  });
  return [...groups,...leafNodes(leaves,compare)];
 };
 return descend(entries,0);
}

export const systemName=(id:SystemId)=>SYSTEMS.find(system=>system.id===id)?.name??id;
