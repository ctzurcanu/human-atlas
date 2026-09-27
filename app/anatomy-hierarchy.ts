import {SYSTEMS,structureName,type Atlas,type Part,type SystemId} from './anatomy';
import {taHierarchyForConcept,terminologyForConcept,terminologyForGroup,type Terminology} from './anatomical-terminology';
import {isSkinPart} from './depth-layers';

export type AnatomyEntry={id:string;name:string;system:SystemId;region:string;location:string[];ancestry:string[];ancestryTerms?:Record<string,Terminology>;parts:Part[];terminology:Terminology};
export type AnatomyNode=
 | {kind:'group';id:string;name:string;nodes:AnatomyNode[];parts:Part[];terminology:Terminology}
 | {kind:'bilateral';id:string;name:string;entries:AnatomyEntry[];parts:Part[];terminology:Terminology}
 | {kind:'entry';id:string;name:string;entry:AnatomyEntry;parts:Part[]};
export const REGION_ORDER=['Head & neck','Torso & pelvis','Upper limbs','Lower limbs','Whole body & spanning'];
const CELL_COMPARTMENTS=['Cell boundary','Nucleus','Cytoplasm'];
export const displaySystemFor=(system:SystemId):SystemId=>system==='regions'?'integumentary':system;
export const MAJOR_SYSTEMS=[
 {id:'respiratory',name:'Respiratory system',source:'respiratory'},
 {id:'digestive',name:'Digestive system',source:'digestive'},
 {id:'circulatory',name:'Circulatory system',source:'cardiac'},
 {id:'urinary',name:'Urinary system',source:'urinary'},
 {id:'integumentary',name:'Integumentary system',source:'integumentary'},
 {id:'skeletal',name:'Skeletal system',source:'skeletal'},
 {id:'muscular',name:'Muscular system',source:'muscular'},
 {id:'endocrine',name:'Endocrine system',source:'endocrine'},
 {id:'exocrine',name:'Exocrine system',source:'digestive'},
 {id:'lymphatic',name:'Lymphatic system',source:'lymphatic'},
 {id:'nervous',name:'Nervous system',source:'nervous'},
 {id:'reproductive',name:'Reproductive system',source:'reproductive'},
] as const;
export type MajorSystemId=typeof MAJOR_SYSTEMS[number]['id'];
const isSkin=(entry:AnatomyEntry)=>entry.parts.every(isSkinPart);
const skinSourceRegions=(entry:AnatomyEntry)=>entry.system==='regions'||entry.parts.some(part=>(part.groups??[]).some(group=>/^9: regions of human body$/i.test(group)));
function integumentaryBranch(entry:AnatomyEntry){
 const name=entry.name.toLowerCase();
 return /mammar|breast|lactifer/.test(name)?'Mammary gland':/lacrimal/.test(name)?'Lacrimal glands':/adipose/.test(name)?'Subcutaneous tissue':'';
}
export function majorSystemFor(entry:Pick<AnatomyEntry,'name'|'system'|'parts'>&Partial<Pick<AnatomyEntry,'terminology'>>):MajorSystemId{
 const name=entry.name.toLowerCase();
 if(entry.parts.every(isSkinPart))return 'integumentary';
 if(/(?:salivary|parotid|submandibular|sublingual|lacrimal|sebaceous|sweat|ceruminous) gland|gland of (?:skin|eyelid)|mammary (?:gland|lobe)|lactiferous/.test(name))return 'exocrine';
 if(entry.system==='cardiac'&&/^(?:third|fourth|lateral) ventricle|interventricular foramen/.test(name))return 'nervous';
 if(entry.system==='skeletal'&&/\b(?:teeth|tooth|gingiva)\b/.test(name))return 'digestive';
 if(entry.system==='digestive'&&/articular disc|articular disk|acromioclavicular|sternoclavicular/.test(name))return 'skeletal';
 const chapter=entry.terminology?.ta98?.slice(0,3);
 const taSystem:Record<string,MajorSystemId>={A02:'skeletal',A03:'skeletal',A04:'muscular',A05:'digestive',A06:'respiratory',A08:'urinary',A09:'reproductive',A11:'endocrine',A12:'circulatory',A13:'lymphatic',A14:'nervous',A15:'nervous',A16:'integumentary'};
 if(chapter&&taSystem[chapter]&&entry.system!=='attachments')return taSystem[chapter];
 if(entry.system==='connective')return /eyeball|(?:lateral|medial) rectus|superior oblique|\blens\b/.test(name)?'nervous':/\bfascia\b/.test(name)?'muscular':'skeletal';
 if(entry.system==='fascia'||entry.system==='attachments'||entry.system==='muscular')return 'muscular';
 if(entry.system==='arterial'||entry.system==='venous'||entry.system==='cardiac')return 'circulatory';
 if(entry.system==='sensory'||entry.system==='nervous'||entry.system==='schematic')return 'nervous';
 if(entry.system==='regions'||entry.system==='integumentary')return 'integumentary';
 if(entry.system==='pregnancy')return 'reproductive';
 return entry.system as MajorSystemId;
}
function systemBranchFor(entry:AnatomyEntry):string{
 const name=entry.name.toLowerCase(),system=majorSystemFor(entry);
 if(system==='integumentary')return isSkin(entry)?'Skin':integumentaryBranch(entry)||(/adipose|fat/.test(name)?'Subcutaneous tissue':/\bhair|nail\b/.test(name)?'Hair & nails':'');
 if(system==='skeletal')return /cartilage|labrum|meniscus|\bdisc\b|\bdisk\b/.test(name)?'Cartilage':/\btendon\b/.test(name)?'Tendons':/ligament|capsule|joint|suture/.test(name)?'Joints and ligaments':'Bones';
 if(system==='muscular')return entry.system==='attachments'?'Muscle attachments':entry.system==='fascia'||/\bfascia\b/.test(name)?'Fascia':'Muscles';
 if(system==='circulatory')return entry.system==='arterial'||/\b(?:artery|aorta|celiac trunk)\b/.test(name)?'Arteries':entry.system==='venous'||/\b(?:vein|venous)\b/.test(name)?'Veins':'Heart';
 if(system==='nervous')return entry.system==='sensory'||entry.system==='connective'||/\b(?:eye|eyeball|retina|sclera|cornea|choroid|pupil|lens|ear|cochlea|macula lutea|conjunctiva|vitreous|lacrimal)\b/.test(name)?'Sensory organs':/^(?:allen )|\b(?:brain|cerebr|cerebell|spinal cord|spinocerebellar|vestibulospinal|tectospinal|rubrospinal|gracile fasciculus|white matter|ventricle|interventricular foramen|thalam|medulla|caudate|putamen|pallid|insula|cortex|accumbens)\b/.test(name)?'Central nervous system':'Peripheral nervous system';
 if(system==='exocrine')return /lacrimal/.test(name)?'Lacrimal glands':/mammary|lactiferous|breast/.test(name)?'Mammary gland':/sebaceous|sweat|ceruminous/.test(name)?'Skin glands':'Salivary glands';
 if(system==='digestive')return /liver|pancrea|gallbladder|biliar/.test(name)?'Accessory digestive organs':'Digestive tract';
 if(system==='respiratory')return /\blung|pulmonary|pleura/.test(name)?'Lungs':'Airways';
 if(system==='urinary')return /kidney|renal/.test(name)?'Kidneys':'Urinary tract';
 if(system==='lymphatic')return /\bnode|nodule\b/.test(name)?'Lymph nodes':'Lymphatic organs and vessels';
 if(system==='endocrine')return 'Endocrine glands';
 if(system==='reproductive')return entry.system==='pregnancy'?'Placenta & pregnancy':/uter|ovar|vagin|clitor|vulva|female/.test(name)?'Female genital system':'Male genital system';
 return '';
}
export function systemPathFor(entry:AnatomyEntry,scope?:string):string[]{
 if(scope==='cell')return [];
 const branch=systemBranchFor(entry);
 if(branch==='Skin')return ['Skin',skinSourceRegions(entry)?'Named skin regions':'Body surface',entry.region,...entry.location,...entry.ancestry];
 return [...(branch?[branch]:[]),entry.region,...entry.location,...entry.ancestry];
}
export function regionPathFor(entry:AnatomyEntry):string[]{return [...entry.location,isSkin(entry)?'Skin':integumentaryBranch(entry)||systemName(entry.system),...entry.ancestry];}
export function depthPathFor(entry:AnatomyEntry,layerName:string):string[]{return [entry.region,...entry.location.filter(name=>name.toLowerCase()!==layerName.toLowerCase()),...entry.ancestry];}
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
 const name=normalized(parts);
 if(/\b(?:acromioclavicular|sternoclavicular|clavicle|scapula|shoulder|deltoid)\b/.test(name)||has(groups,/pectoral girdle/))return 'Upper limbs';
 if(has(groups,/^(head|neck|brain|face|skull|cranium|cerebrum|cerebellum|left head|right head|head and neck)$/))return 'Head & neck';
 if(has(groups,/^(left |right )?(upper limb|arm|forearm|hand|wrist|shoulder)( region)?$/))return 'Upper limbs';
 if(has(groups,/^(left |right )?(lower limb|leg|thigh|foot|ankle|talocrural region)( region)?$/))return 'Lower limbs';
 if(has(groups,/^(trunk|torso|thorax|chest|abdomen|pelvis|pelvic region)$/))return 'Torso & pelvis';
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
  return side?[`${side} leg`,area]:[];
 }
 if(region==='Upper limbs'){
  const area=/\b(hand|finger|thumb|digit|metacarp|carpal|wrist|palm)\b/.test(name)||has(groups,/^(left |right )?hand$/)?'Hand & wrist'
   :/\b(forearm|radius|ulna|radial|ulnar|elbow)\b/.test(name)?'Forearm & elbow'
   :/\b(shoulder|scapula|clavicle|clavicular|acromioclavicular|pectoral girdle|axilla)\b/.test(name)?'Shoulder'
   :/\b(arm|humerus|brachial|biceps|triceps)\b/.test(name)?'Upper arm':y<.94?'Hand & wrist':y<1.12?'Forearm & elbow':y<1.4?'Upper arm':'Shoulder';
  return side?[`${side} arm`,area]:[];
 }
 if(region==='Head & neck'){
  if(/\b(neck|cervical|larynx|pharynx|thyroid|trachea)\b/.test(name)||has(groups,/^neck$/))return ['Neck'];
  const area=/\b(brain|cerebr|cerebell|telencephal|ventric|thalam|hypothalam|pons|medulla|pineal|pituitary|caudate|putamen|pallid|insula|cortex|nucleus accumbens)\b/.test(name)||/^allen /.test(name)||has(groups,/^(brain|cerebrum|telencephalon)$/)?'Brain'
   :/\b(skull|crani|parietal|occipital|frontal bone|temporal bone|sphenoid|ethmoid|mandib|maxill|zygomat|vomer)\b/.test(name)||has(groups,/^cranium$/)?'Skull'
   :/\b(scalp|epicran|mastoid|temporal region|hairs of head)\b/.test(name)||has(groups,/^regions of epicranium$/)?'Scalp & epicranium'
   :/\b(eye|eyebrow|eyelash|auricul\w*|antihelix|antitragus|tragus|helix|concha|conchae|scapha|fossa antihelica|ear|nose|nasal|mouth|oral|labial|philtrum|mentolabial|mental region|parotideomasseteric|buccal|zygomatic region|frontal region|tongue|tooth|teeth|facial|face|orbit|cheek|lip|chin)\b/.test(name)||has(groups,/^(face|regions of face|auricular region|oral region|orbital region|frontal region|hairs)$/)?'Face & senses':'';
  return area?['Head',area]:['Head'];
 }
 if(region==='Torso & pelvis'){
  const area=/\b(back|vertebr\w*|spine|spinal|lumbar|dorsal|scapular|infrascapular|interscapular|sacral region|triangle of auscultation|erector spinae)\b/.test(name)||has(groups,/^back$/)?'Back & spine'
   :/\b(pelvi\w*|sacrum|coccyx|bladder|uterus|ovary|prostate|rectum|perine\w*|genital\w*|urethra|vagina|testis|penis|pubic|anal region|urogenital)\b/.test(name)||has(groups,/^(pelvis|pelvic region)$/)?'Pelvis'
   :/\b(abdom\w*|epigastric|hypochondriac|hypogastric|umbilic\w*|inguinal|kidney|renal|liver|stomach|intestin\w*|colon|spleen|pancreas|duodenum|jejunum|ileum)\b/.test(name)||has(groups,/^(abdomen|regions of abdomen)$/)?'Abdomen'
   :/\b(chest|thorax|thorac\w*|rib|stern\w*|lung|heart|bronch\w*|mammar\w*|inframammary|pectoral|presternal|infraclavicular|deltopectoral|intercost\w*)\b/.test(name)||has(groups,/^(thorax|chest|regions of thorax|pectoral region|anterior region of thorax)$/)?'Chest':'';
  return area?[area]:[];
 }
 return [];
}

// The reference GLBs retain source region groups. Map their anatomical names
// into shared, optional skin branches; models without these groups stop earlier.
const SKIN_DETAIL_GROUPS:[RegExp,string[]][]=[
 [/^concha of auricle$/,['Auricular region','Concha of auricle']],
 [/^auricular region$/,['Auricular region']],
 [/^oral region$/,['Oral region']],
 [/^orbital region$/,['Orbital region']],
 [/^nasal region$/,['Nasal region']],
 [/^frontal region$/,['Frontal region']],
 [/^regions of face$/,['Facial regions']],
 [/^regions of epicranium$/,['Epicranial regions']],
 [/^hairs$/,['Hair']],
 [/^sternocleidomastoid region$/,['Sternocleidomastoid region']],
 [/^lateral region of neck$/,['Lateral neck']],
 [/^suprahyoid region$/,['Suprahyoid region']],
 [/^infrahyoid region$/,['Infrahyoid region']],
 [/^pectoral region$/,['Pectoral region']],
 [/^anterior region of thorax$/,['Anterior thorax']],
 [/^umbilical region$/,['Umbilical region']],
 [/^epigastric region$/,['Epigastric region']],
 [/^anal region$/,['Anal region']],
 [/^regions of perineum$/,['Perineal regions']],
 [/^gluteal region$/,['Gluteal region']],
 [/^anterior region of arm$/,['Anterior arm']],
 [/^brachial region$/,['Brachial region']],
 [/^anterior region of elbow$/,['Anterior elbow']],
 [/^cubital region$/,['Cubital region']],
 [/^antebrachial region$/,['Antebrachial region']],
 [/^carpal region$/,['Carpal region']],
 [/^region of digits of hand$/,['Digits of hand']],
 [/^palm$/,['Palm']],
 [/^region of hand$/,['Hand regions']],
 [/^anterior region of thigh$/,['Anterior thigh']],
 [/^femoral region$/,['Femoral region']],
 [/^posterior region of knee$/,['Posterior knee']],
 [/^knee region$/,['Knee region']],
 [/^posterior region of leg$/,['Posterior leg']],
 [/^leg region$/,['Leg regions']],
 [/^posterior region of ankle$/,['Posterior ankle']],
 [/^talocrural region$/,['Talocrural region']],
 [/^longitudinal arch of foot$/,['Longitudinal arch']],
 [/^metatarsal region$/,['Metatarsal region']],
 [/^region of digits of foot$/,['Digits of foot']],
 [/^sole$/,['Sole']],
 [/^regions of foot$/,['Foot regions']],
];
export const SKIN_DETAIL_NAMES=[...new Set(SKIN_DETAIL_GROUPS.flatMap(([,path])=>path))];
function skinDetailFor(parts:Part[]):string[]{
 const groups=groupsFor(parts);
 return SKIN_DETAIL_GROUPS.find(([pattern])=>has(groups,pattern))?.[1]??[];
}

const GENERIC_TA_PARENTS=new Set(['human body','the integument','alimentary system','muscles','bones','joints','arteries','veins','nerves','regions of human body','lymphatic system','nervous system','peripheral nervous system','central nervous system','respiratory system','digestive system','circulatory system','urinary system','integumentary system','skeletal system','muscular system','endocrine system','exocrine system','reproductive system']);
const GROUP_CHAPTERS:Record<MajorSystemId,string[]>={respiratory:['A06'],digestive:['A05'],circulatory:['A12'],urinary:['A08'],integumentary:['A01','A16'],skeletal:['A02','A03'],muscular:['A04'],endocrine:['A11'],exocrine:['A05','A15','A16'],lymphatic:['A13'],nervous:['A14','A15'],reproductive:['A09']};
const normalizedGroup=(name:string)=>name.toLowerCase().replace(/\bmuscle\b/g,'').replace(/[^a-z0-9]+/g,' ').trim();
function ancestryFor(entry:AnatomyEntry,groupCounts:Map<string,number>):string[]{
 const chapter=entry.terminology.ta98?.slice(0,3);
 const allowed=chapter?[chapter]:GROUP_CHAPTERS[majorSystemFor(entry)];
 const used=new Set([entry.name,...entry.location,entry.region,systemBranchFor(entry),systemName(entry.system),MAJOR_SYSTEMS.find(system=>system.id===majorSystemFor(entry))?.name??''].map(normalizedGroup));
 const usedCodes=new Set<string>();
 const result:string[]=[];
 const append=(name:string,code:string|null)=>{
  const key=normalizedGroup(name);
  if(!key||used.has(key)||code&&usedCodes.has(code))return;
  used.add(key);if(code)usedCodes.add(code);result.push(name);
 };
 for(const {name,terminology} of taHierarchyForConcept(entry.id)){
  if(GENERIC_TA_PARENTS.has(name.toLowerCase())||!allowed.includes((terminology.ta98??'').slice(0,3)))continue;
  append(name,terminology.ta98);
 }
 const source=[...new Set(entry.parts.flatMap(part=>part.groups??[]))]
  .filter(name=>!/^\d+:\s/.test(name)&&!GENERIC_TA_PARENTS.has(name.toLowerCase())&&!/^(?:left|right) (?:upper|lower) limb$/i.test(name))
  .map(name=>({name,term:terminologyForGroup(name)}))
  .filter(({term})=>!!term.ta98&&allowed.includes(term.ta98.slice(0,3))&&term.ta98!==entry.terminology.ta98)
  .sort((a,b)=>(groupCounts.get(b.name)??0)-(groupCounts.get(a.name)??0)||a.name.localeCompare(b.name));
 for(const {name,term} of source)append(name,term.ta98);
 return result;
}

export function hierarchyEntries(atlas:Atlas):AnatomyEntry[]{
 const names=new Map(atlas.concepts.map(concept=>[concept.id,concept.name]));
 const byConcept=new Map<string,Part[]>();
 const groupCounts=new Map<string,number>();
 for(const part of atlas.parts){if(part.suppressed)continue;const parts=byConcept.get(part.conceptId)??[];parts.push(part);byConcept.set(part.conceptId,parts);for(const group of new Set(part.groups??[]))groupCounts.set(group,(groupCounts.get(group)??0)+1);}
 const entries=[...byConcept].map(([id,parts])=>{
  const region=regionFor(parts,atlas.scope);
  const location=atlas.scope==='cell'?[]:locationFor(parts,region);
  if(parts.every(isSkinPart))location.push(...skinDetailFor(parts));
  const entry={id,name:names.get(id)??structureName(parts[0].name),system:parts[0].system,region,location,ancestry:[],parts,terminology:terminologyForConcept(id)} as AnatomyEntry;
  entry.ancestry=atlas.scope==='cell'?[]:ancestryFor(entry,groupCounts);
  return entry;
 });
 const muscleKey=(name:string)=>name.toLowerCase().replace(/\s*\((?:left|right)\)\s*$/,'').replace(/^\((.*)\)$/,'$1').trim();
 const muscles=new Map<string,AnatomyEntry>();
 for(const entry of entries)if(entry.system==='muscular'){
  const side=/\((left|right)\)\s*$/i.exec(entry.name)?.[1].toLowerCase();
  if(side)muscles.set(`${muscleKey(entry.name)}:${side}`,entry);
 }
 for(const entry of entries)if(entry.system==='attachments'){
  const marker=/^(.+)\.[eo]\d*([lr])$/i.exec(entry.name);
  if(!marker)continue;
  const parent=muscles.get(`${muscleKey(marker[1])}:${marker[2].toLowerCase()==='l'?'left':'right'}`);
  if(!parent)continue;
  const parentName=parent.name.replace(/\s*\((?:left|right)\)\s*$/i,'').replace(/^\((.*)\)$/,'$1').trim();
  entry.region=parent.region;entry.location=[...parent.location];
  entry.ancestry=[...parent.ancestry,parentName];
  entry.ancestryTerms={[parentName]:parent.terminology};
 }
 return entries;
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
   const term=siblings[0].terminology;
   nodes.push({kind:'bilateral',id:`pair:${parsed.base}`,name:parsed.base,entries:siblings.sort((a,b)=>entryLabel(a).localeCompare(entryLabel(b))),parts:siblings.flatMap(sibling=>sibling.parts),terminology:siblings.every(sibling=>sibling.terminology.ta98===term.ta98&&sibling.terminology.latin===term.latin)?term:terminologyForGroup(parsed.base)});
  }else nodes.push({kind:'entry',id:`entry:${entry.id}`,name:entryLabel(entry),entry,parts:entry.parts});
 }
 return nodes.sort((a,b)=>Number(b.kind==='bilateral')-Number(a.kind==='bilateral')||(compare&&a.kind==='entry'&&b.kind==='entry'?compare(a.entry,b.entry):0)||a.name.localeCompare(b.name,undefined,{numeric:true,sensitivity:'base'}));
}

/** Every entry has one path. Leaves without another named parent stay at this level. */
export function buildAnatomyNodes(entries:AnatomyEntry[],path:(entry:AnatomyEntry)=>string[],compare?:(a:AnatomyEntry,b:AnatomyEntry)=>number):AnatomyNode[]{
 const paths=new Map(entries.map(entry=>[entry,path(entry)]));
 const descend=(items:AnatomyEntry[],depth:number):AnatomyNode[]=>{
  const byGroup=new Map<string,AnatomyEntry[]>(),leaves:AnatomyEntry[]=[];
  for(const entry of items){
   const name=paths.get(entry)?.[depth];
   if(!name){leaves.push(entry);continue;}
   const group=byGroup.get(name)??[];group.push(entry);byGroup.set(name,group);
  }
  const groups:AnatomyNode[]=[];
  for(const [name,members] of byGroup){
   const parents=members.map(entry=>entry.ancestryTerms?.[name]??taHierarchyForConcept(entry.id).find(parent=>parent.name.toLowerCase()===name.toLowerCase())?.terminology);
   const common=parents[0]?.ta98&&parents.every(parent=>parent?.ta98===parents[0]?.ta98)?parents[0]:null;
   groups.push({kind:'group',id:`group:${depth}:${name}`,name,nodes:descend(members,depth+1),parts:members.flatMap(entry=>entry.parts),terminology:common??terminologyForGroup(name)});
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
