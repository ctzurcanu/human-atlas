import assert from 'node:assert/strict';
import {existsSync,readdirSync,readFileSync,writeFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import {build} from 'esbuild';

export async function loadChakraInputs(){
 const bundle=await build({stdin:{contents:`export {hierarchyEntries,majorSystemFor} from './app/anatomy-hierarchy'; export {anatomicalRelations} from './app/anatomical-relations'; export {isSkinPart} from './app/depth-layers'; export {guestAnatomyBase,parseGuestHierarchy,resolveGuestHierarchy} from './app/guest-hierarchy';`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
 const api=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
 const catalogues=[];
 for(const folder of ['public/models','.local-models'])if(existsSync(folder))for(const file of readdirSync(folder).filter(name=>name.endsWith('.json')&&!name.endsWith('.vr.json')).sort()){
  const atlas=JSON.parse(readFileSync(`${folder}/${file}`,'utf8'));
  if(Array.isArray(atlas.parts)&&Array.isArray(atlas.concepts)&&atlas.scope!=='cell')catalogues.push({file:`${folder}/${file}`,atlas});
 }
 return {...api,catalogues,metadata:JSON.parse(readFileSync('app/data/ta98-metadata.json','utf8'))};
}

// Chakra placements follow the user's outline. TA98 and the atlas relation graph
// identify anatomical branches; they do not establish chakra associations.
export function nerveGroups(entry,metadata){
 const name=entry.name.toLowerCase(),path=metadata.byConceptPath[entry.id]??[];
 const groups=entry.parts.flatMap(part=>part.groups??[]).join(' ').toLowerCase();
 const has=code=>path.includes(code);
 const knownSpinal=path.some(code=>/^A14\.2\.0[2-7]\./.test(code))||/spinal (?:nerve|ganglion)/.test(name);
 if(/cauda equina/.test(name))return ['sacral'];
 if(/vagus|vagal|laryngeal nerve|laryngeal branch/.test(name)||has('A14.2.01.153'))return ['vagus'];
 if(has('A14.2.01.121')||/vestibulocochlear|cochlear nerve|vestibular nerve/.test(name))return ['hearing'];
 // CNS, sensory organs, and unrelated cranial nerves are not spinal branches.
 if(path.some(code=>code.startsWith('A14.1.')||code.startsWith('A15.'))||has('A14.2.01.001')||/cranial nerves|central nervous system/.test(groups)&&!knownSpinal)return [];
 if(/olfactory|optic |oculomotor|trochlear|trigeminal|ophthalmic|abducens|facial nerve|glossopharyngeal|hypoglossal|accessory nerve|alveolar|mandibular|lingual nerve|buccal nerve|frontal nerve|auriculotemporal|infraorbital|supraorbital|nasociliary|zygomatic|ciliary gangli|chorda tympani|pterygopalatine|otic gangli|submandibular gangli/.test(name))return [];
 const level=/\b(c|t|l|s|co)\s*0?(\d{1,2})\b/i.exec(name);
 if(/spinal nerve/.test(name)&&level){
  const n=Number(level[2]);
  if(level[1]==='c')return n<=8?['cervical',...(n>=5?['brachial']:[])]:[];
  if(level[1]==='t')return [...(n<=6?['thoracicUpper']:[]),...(n>=6&&n<=9?['thoracicMiddle']:[]),...(n>=10&&n<=12?['thoracicLower']:[])];
  if(level[1]==='l')return ['lumbar'];
  if(level[1]==='s'||level[1]==='co')return ['sacral'];
 }
 const intercostal=/intercostal nerve\s*0?(\d{1,2})\b/.exec(name);
 if(intercostal){const n=Number(intercostal[1]);return [...(n<=6?['thoracicUpper']:[]),...(n>=6&&n<=9?['thoracicMiddle']:[]),...(n>=10&&n<=12?['thoracicLower']:[])];}
 const text=`${name} ${groups}`;
 if(path.some(code=>code.startsWith('A14.2.03.'))||/brachial plexus|median nerve|ulnar nerve|radial nerve|axillary nerve|musculocutan(?:eous|eus)|pectoral nerve|subscapular nerve|thoracodorsal|long thoracic|scapular nerve|brachial cutaneous|antebrachial cutaneous|interosseous nerve of forearm|palmar digital|(?:forearm|arm)\s*-\s*nerves/.test(text))return ['brachial'];
 if(has('A14.2.07.027')||has('A14.2.06.001')||/sacral|coccyg|sciatic|tibial nerve|fibular nerve|peroneal nerve|plantar|sural|pudendal|gluteal nerve|posterior femoral cutaneous|nerve to (?:piriformis|obturator internus|quadratus femoris|levator ani)|cauda equina/.test(text))return ['sacral'];
 if(has('A14.2.07.002')||has('A14.2.05.001')||/lumbar|femoral nerve|obturator nerve|saphenous|genitofemoral|iliohypogastric|ilio.?inguinal|femoral cutaneous/.test(text))return ['lumbar'];
 if(has('A14.2.02.002')||/cervical plexus|cervical nerve|phrenic|ansa cervicalis|occipital nerve|great auricular|supraclavicular nerve|transverse cervical nerve/.test(text))return ['cervical'];
 if(has('A14.2.04.001')||/intercostal|subcostal|thoracic nerve/.test(text))return ['thoracicUpper','thoracicMiddle','thoracicLower'];
 if(/(?:greater|lesser) splanchnic|celiac|coeliac|gastric plexus|hepatic plexus|splenic plexus/.test(text))return ['thoracicMiddle'];
 if(/least splanchnic|mesenteric|suprarenal plexus|aorticorenal|renal plexus/.test(text))return ['thoracicLower'];
 if(/hypogastric|pelvic splanchnic|gonadal|rectal gangli/.test(text))return ['lumbar','sacral'];
 if(/spinal (?:nerve|ganglion)|sympathetic|paravertebral|rami communicantes/.test(text)||path.some(code=>/^A14\.2\.(?:0[2-7])\./.test(code)))return ['unspecified'];
 return [];
}

export function buildChakraHierarchy(input){
 const {catalogues,metadata,hierarchyEntries,majorSystemFor,isSkinPart,guestAnatomyBase,anatomicalRelations,resolveGuestHierarchy,parseGuestHierarchy}=input;
 const hierarchy=JSON.parse(readFileSync('scripts/chakras-outline.json','utf8'));
 hierarchy.source='https://hackmd.io/21uN7y70SPOLKgK_mYZ-mA';
 const buckets=new Map();
 const node=(key,name,parent)=>{const item={name,matchBases:[]};parent.children??=[];parent.children.push(item);buckets.set(key,{item,names:new Set()});return item;};
 const [ajna,vish,ana,upper,lower,swa,mula]=hierarchy.nodes;
 node('cervical','Cervical spinal nerves and branches',vish);
 node('vagus','Vagus nerve and branches',vish);
 node('hearing','Inner ear and hearing',vish);
 node('brachial','Brachial plexus and branches',ana);
 node('thoracicUpper','Thoracic spinal nerves · T1–6 and unnumbered branches',ana);
 node('thoracicMiddle','Thoracic spinal nerves · T6–9 and unnumbered branches',upper);
 node('thoracicLower','Thoracic spinal nerves · T10–12 and unnumbered branches',lower);
 node('lumbar','Lumbar spinal nerves and branches',swa);
 node('sacral','Sacral and coccygeal spinal nerves and branches',mula);
 // Pelvic innervation is shared by these two roots in the source outline.
 node('pelvic','Sacral and pelvic nerves',swa);
 node('unspecified','Spinal and autonomic nerves · level unspecified',mula);
 node('skin','Skin · tactile',ana);
 node('respiratory','Lungs and airways',ana);
 node('circulatory','Heart and circulatory system',ana);
 node('digestive','Digestive system',upper);
 node('muscular','Muscles',upper);
 node('skeletal','Bones and skeleton',mula);
 node('lymphatic','Lymphatic system',swa);
 const rootKeys=new Map([[vish,['cervical','vagus','hearing']],[ana,['brachial','thoracicUpper']],[upper,['thoracicMiddle']],[lower,['thoracicLower']],[swa,['lumbar','pelvic']],[mula,['sacral','unspecified']]]);
 for(const [root] of rootKeys)node(`targets:${hierarchy.nodes.indexOf(root)}`,'Innervated organs and tissues',root);
 const add=(key,name)=>buckets.get(key).names.add(guestAnatomyBase(name));
 let nerveCount=0,targetCount=0;
 for(const {file,atlas} of catalogues){
  const entries=hierarchyEntries(atlas),byId=new Map(entries.map(entry=>[entry.id,entry]));
  const parts=atlas.parts.filter(part=>!part.suppressed).map(part=>({...part,name:byId.get(part.conceptId)?.name??part.name})),byPart=new Map(parts.map(part=>[part.id,byId.get(part.conceptId)]));
  const membership=new Map();
  for(const entry of entries){
   const sys=majorSystemFor(entry);
   if(buckets.has(sys)&&sys!=='nervous')add(sys,entry.name);
   if(entry.parts.some(isSkinPart))add('skin',entry.name);
   const path=metadata.byConceptPath[entry.id]??[];
   if(path.some(code=>code.startsWith('A15.3.03.'))||/cochle|organ of corti|spiral organ|inner ear|internal ear|vestibul|semicircular|endolymph|perilymph|spiral ganglion|auditory/.test(entry.name.toLowerCase()))add('hearing',entry.name);
   if(sys!=='nervous')continue;
   const keys=nerveGroups(entry,metadata);
   if(keys.length){membership.set(entry.id,new Set(keys));nerveCount++;}
  }
  const relations=new Map();
  const links=entry=>{
   if(!relations.has(entry.id))relations.set(entry.id,entry.parts.flatMap(part=>anatomicalRelations({...part,name:entry.name},parts)));
   return relations.get(entry.id);
  };
  // Follow only directed neural branches, never undirected 'related' anatomy.
  const queue=[...membership.keys()];
  for(let index=0;index<queue.length;index++){
   const entry=byId.get(queue[index]),keys=membership.get(entry.id);
   for(const link of links(entry)){
    if(!['after','contains'].includes(link.kind))continue;
    const target=byPart.get(link.target.id);
    if(!target||majorSystemFor(target)!=='nervous'||!(/nerve|plexus|gangli|ramus|rami|branch|root/i.test(target.name)))continue;
    const targetPath=metadata.byConceptPath[target.id]??[];
    if(targetPath.some(code=>code.startsWith('A14.1.')||code.startsWith('A15.')))continue;
    const targetKeys=membership.get(target.id)??new Set();
    let changed=false;for(const key of keys)if(!targetKeys.has(key)){targetKeys.add(key);changed=true;}
    if(changed){membership.set(target.id,targetKeys);queue.push(target.id);}
   }
  }
  for(const [id,keys] of membership)for(const key of keys){add(key,byId.get(id).name);if(key==='sacral')add('pelvic',byId.get(id).name);}
  // Include the original named plexuses, even when only their source label is known.
  const resolved=resolveGuestHierarchy(atlas,parseGuestHierarchy(hierarchy));
  const rootNerves=resolved.map(root=>new Set(root.parts.map(part=>part.conceptId).filter(id=>majorSystemFor(byId.get(id))==='nervous')));
  for(const [root,keys] of rootKeys){
   const rootIndex=hierarchy.nodes.indexOf(root),targets=new Set();
   const ids=new Set([...rootNerves[rootIndex],...[...membership].filter(([,found])=>keys.some(key=>found.has(key))||keys.includes('pelvic')&&found.has('sacral')).map(([id])=>id)]);
   for(const id of ids)for(const link of links(byId.get(id)))if(link.kind==='innervates'){
    const target=byPart.get(link.target.id);if(target){targets.add(target.id);targetCount++;}
   }
   // Named organ associations in the outline and the peripheral nerve reference
   // supplement the atlas's incomplete functional relation table.
   const organNames=root===vish
    ? /^(?:diaphragm|pharynx|larynx|trachea|(?:right |left )?main bronchus|lung|heart|[eo]sophagus|stomach|liver|pancreas|duodenum|jejunum|ileum|caecum|cecum|ascending colon|transverse colon)(?:\b|$)/i
    : null;
   for(const entry of entries){
    if(organNames?.test(guestAnatomyBase(entry.name)))targets.add(entry.id);
    if((root===swa||root===mula)&&['urinary','reproductive'].includes(majorSystemFor(entry)))targets.add(entry.id);
   }
   // Preserve entire organs: their TA98 descendants and material pieces follow them.
   const parentCodes=new Set([...targets].map(id=>metadata.byConceptMatch[id]?.term).filter(Boolean));
   for(const entry of entries)if(targets.has(entry.id)||(metadata.byConceptPath[entry.id]??[]).some(code=>parentCodes.has(code)))add(`targets:${rootIndex}`,entry.name);
  }
  console.log(`${file}: ${membership.size} spinal/vagus/hearing nerve concepts classified`);
 }
 // Ensure all named levels remain covered even when a catalogue omits a mesh.
 for(let i=1;i<=8;i++)add('cervical',`Cervical spinal nerve c${i}`);
 for(let i=1;i<=12;i++){const keys=i<=5?['thoracicUpper']:i===6?['thoracicUpper','thoracicMiddle']:i<=9?['thoracicMiddle']:['thoracicLower'];for(const key of keys){add(key,`Thoracic spinal nerve t${i}`);add(key,`Intercostal nerve${String(i).padStart(2,'0')}`);}}
 for(let i=1;i<=5;i++){add('lumbar',`Lumbar spinal nerve l${i}`);add('sacral',`Sacral spinal nerve s${i}`);add('pelvic',`Sacral spinal nerve s${i}`);}
 add('sacral','Coccygeal spinal nerve co1');add('pelvic','Coccygeal spinal nerve co1');
 for(const {item,names} of buckets.values())item.matchBases=[...names].sort();
 parseGuestHierarchy(hierarchy);
 const text=JSON.stringify(hierarchy,null,2)+'\n';
 assert(text.length<500_000,'Chakras must fit the existing guest hierarchy size limit');
 return {hierarchy,text,nerveCount,targetCount};
}

if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){
 const {text,nerveCount,targetCount}=buildChakraHierarchy(await loadChakraInputs());
 if(process.argv.includes('--check'))assert.equal(readFileSync('public/assets/chakras.json','utf8'),text,'Run npm run build:chakras to refresh the hierarchy');
 else writeFileSync('public/assets/chakras.json',text);
 console.log(`Chakras: ${text.length} bytes; ${nerveCount} nerve seeds and ${targetCount} modeled innervation links across catalogues.`);
}
