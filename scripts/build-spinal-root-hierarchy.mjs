import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';
import {createReportViewerLinks,markdownLabel} from './report-viewer-links.mjs';

export const SPINAL_REGIONS=[['C','Cervical',8],['T','Thoracic',12],['L','Lumbar',5],['S','Sacral',5],['Co','Coccygeal',1]];
export const SPINAL_LEVELS=SPINAL_REGIONS.flatMap(([prefix,,count])=>Array.from({length:count},(_,i)=>`${prefix}${i+1}`));
export const ROOT_MAP=JSON.parse(readFileSync(new URL('./spinal-root-map.json',import.meta.url),'utf8'));
const sideOf=name=>name.match(/\((left|right)\)\s*$/i)?.[1]?.toLowerCase()??name.match(/^(left|right)\s/i)?.[1]?.toLowerCase();
const base=name=>name.replace(/\s*\((?:left|right)\)\s*$/i,'').replace(/^(?:left|right)\s+/i,'').replace(/^\((.+)\)$/,'$1').trim().toLowerCase();
const canonicalLevel=(prefix,n)=>`${prefix.toLowerCase()==='co'?'Co':prefix.toUpperCase()}${Number(n)}`;
export function rootLevel(name){
 if(!/\b(?:root|spinal nerve|spinal ganglion|intercostal nerve|subcostal nerve)\b/i.test(name))return;
 const level=name.match(/\b(co|c|t|l|s)\s*0?(\d{1,2})\b/i);
 if(level){const id=canonicalLevel(level[1],level[2]);return SPINAL_LEVELS.includes(id)?id:undefined;}
 const intercostal=name.match(/intercostal nerve\s*0?(\d{1,2})\b/i);if(intercostal&&Number(intercostal[1])<=11)return `T${Number(intercostal[1])}`;
 if(/subcostal nerve/i.test(name))return 'T12';
}
const ordinals=['first','second','third','fourth','fifth','sixth','seventh','eighth','ninth','tenth','eleventh','twelfth'];
function segmentLevel(name){
 if(!/spinal cord segment|segment of .*spinal cord/i.test(name))return;
 const numbered=name.match(/\b(co|c|t|l|s)\s*(\d{1,2})\b/i);if(numbered)return canonicalLevel(numbered[1],numbered[2]);
 for(const [prefix,region] of SPINAL_REGIONS)if(name.toLowerCase().includes(region.toLowerCase())){const n=ordinals.findIndex(word=>name.toLowerCase().startsWith(word+' '))+1;if(n)return `${prefix}${n}`;}
}
const rootNames=(level,side)=>{
 const region=SPINAL_REGIONS.find(([prefix])=>prefix===level.replace(/\d+$/,''))[1];
 return [`${level} root`,`${region} spinal nerve ${level.toLowerCase()}`,`Spinal nerve ${level}`,`${level} spinal nerve`,`${level} anterior root`,`${level} posterior root`,`Anterior root of spinal nerve ${level}`,`Posterior root of spinal nerve ${level}`,`Anterior root of ${level} spinal nerve`,`Posterior root of ${level} spinal nerve`,`Dorsal root ${level}`,`Ventral root ${level}`].map(name=>`${name} (${side})`);
};

export function buildSpinalRootHierarchy(input){
 const entries=input.catalogues.flatMap(({atlas})=>input.hierarchyEntries(atlas).map(entry=>({...entry,parts:entry.parts.filter(part=>!part.suppressed)}))).filter(entry=>entry.parts.length);
 const nodes=[],byId=new Map(),add=node=>{assert(!byId.has(node.id));byId.set(node.id,node);nodes.push(node);return node;};
 const roots=[],unsegmented=new Set();
 for(const [prefix,region,count] of SPINAL_REGIONS){
  const group=add({id:`REGION:${prefix}`,name:`${region} · ${prefix}1–${prefix}${count}`,childIds:[]});roots.push(group.id);
  for(let n=1;n<=count;n++){
   const level=`${prefix}${n}`,parent=add({id:`ROOT:${level}`,name:`${level} · spinal nerve roots`,synonyms:[`${region} spinal nerve ${n}`,`${level} root`],description:`${level}: bilateral spinal nerve level. Root anatomy, dermatomal territory and representative myotomal targets are listed separately. Muscles and sensory territories overlap adjacent levels.`,source:ROOT_MAP.sources.roots.url,childIds:[]});group.childIds.push(parent.id);
   const segmentMatches=[...new Set(entries.filter(entry=>entry.parts.some(part=>part.system==='nervous')&&segmentLevel(entry.name)===level).map(entry=>entry.id))].sort();
   if(segmentMatches.length){const segment=add({id:`SEGMENT:${level}`,name:`${level} · spinal cord segment reference`,description:'Spinal cord segment, not a peripheral nerve root or a vertebra. Shown only as a central anatomical reference.',source:ROOT_MAP.sources.roots.url,matches:segmentMatches});parent.childIds.push(segment.id);}
   for(const side of ['left','right']){
    const branch=add({id:`ROOT:${level}:${side}`,name:`${level} (${side})`,childIds:[],description:`${level} ${side}. Unsided source muscle meshes may span both sides; sided meshes are selected only on this side.`});parent.childIds.push(branch.id);
    const found=entries.filter(entry=>entry.parts.some(part=>part.system==='nervous')&&rootLevel(entry.name)===level&&sideOf(entry.name)===side);
    const nerve=add({id:`NERVE:${level}:${side}`,name:`${level} · root and segmental nerve (${side})`,synonyms:[`${level} nerve root ${side}`],description:'Selects separately named roots, spinal nerves, ganglia and numbered segmental nerves. Unnumbered spinal root bundles are kept under “Unsegmented spinal roots” rather than assigned to every level.',source:ROOT_MAP.sources.roots.url,matches:[...new Set([...rootNames(level,side),...found.map(entry=>entry.id)])].sort()});branch.childIds.push(nerve.id);
    const sensorySource=level==='C1'?ROOT_MAP.sources.skin.url:level==='Co1'?ROOT_MAP.sources.coccygeal.url:ROOT_MAP.sources.isncsci.url;
    const sensory=add({id:`DERMATOME:${level}:${side}`,name:`Dermatome · ${level} (${side})`,synonyms:[`${level} sensory territory ${side}`],description:`${ROOT_MAP.sensory[level]} Territories are approximate and overlap. This atlas has no separately segmented dermatomal skin patches; nerve anatomy below is a reference, not a colored skin territory.`,source:sensorySource,childIds:[nerve.id]});branch.childIds.push(sensory.id);
    const key=ROOT_MAP.keyActions[level],myo=add({id:`MYOTOME:${level}:${side}`,name:`Myotome · ${level}${key?` · ${key.toLowerCase()}`:''} (${side})`,description:key?`ISNCSCI representative key action: ${key}. Additional listed muscle associations preserve their overlapping root contributions; they are not exclusive single-root maps.`:`No isolated ISNCSCI key muscle at ${level}. Listed muscles have documented shared segmental contributions; absence of a key muscle does not mean absence of motor fibers.`,source:ROOT_MAP.sources.isncsci.url,childIds:[]});branch.childIds.push(myo.id);
    const addMuscles=(mapping,index)=>{
     const aliases=new Set(mapping.aliases.map(base));
     const matches=[...new Set(entries.filter(entry=>(sideOf(entry.name)===side||!sideOf(entry.name))&&aliases.has(base(entry.name))).flatMap(entry=>entry.parts.filter(part=>part.system==='muscular'&&!/tendon|\.[eo]\d*[lr]$/i.test(part.name)).map(part=>part.id)))].sort();
     const item=add({id:`MUSCLE:${level}:${side}:${index}`,name:`${mapping.name} · ${mapping.roots.join('/')}`,description:`Documented root contributions: ${mapping.roots.join(', ')}. ${mapping.note??'Whole-muscle anatomy is shared by several roots; this is not an exclusive innervation territory.'} Source muscle surfaces may include tendon-coded materials.`,source:ROOT_MAP.sources[mapping.source].url,matches});myo.childIds.push(item.id);
    };
    ROOT_MAP.muscles.forEach((mapping,index)=>{if(mapping.roots.includes(level))addMuscles(mapping,index);});
    if(prefix==='T'){
     const number=n,labels=number<=11?['External intercostal muscles','Internal intercostal muscles','Innermost intercostal muscles']:['External abdominal oblique muscle','Internal abdominal oblique muscle','Transversus abdominis muscle','Rectus abdominis muscle'];
     if(number>=7&&number<=11)labels.push('External abdominal oblique muscle','Internal abdominal oblique muscle','Transversus abdominis muscle','Rectus abdominis muscle');
     const matches=[...new Set(entries.filter(entry=>(sideOf(entry.name)===side||!sideOf(entry.name))&&labels.map(base).includes(base(entry.name))).flatMap(entry=>entry.parts.filter(part=>part.system==='muscular'&&!/tendon/i.test(part.name)).map(part=>part.id)))].sort();
     const item=add({id:`AXIAL:${level}:${side}`,name:`${number<=11?'Intercostal':'Subcostal'} motor targets · shared source meshes`,description:`${level} supplies segmental body-wall muscles. The imported intercostal/abdominal muscle surfaces bundle several segments, so selecting this reference selects the whole available muscle surface, not just the ${level} motor strip.`,source:ROOT_MAP.sources.thoracic.url,matches});myo.childIds.push(item.id);
    }
    if(level==='Co1'){
     myo.description+=' Co1 participates in the S4–Co1 coccygeal plexus. An isolated Co1 muscle territory is not established here.';
     myo.source=ROOT_MAP.sources.coccygeal.url;
    }
   }
  }
 }
 for(const entry of entries)if(entry.parts.some(part=>part.system==='nervous')&&!rootLevel(entry.name)&&/anterior root of spinal nerve|posterior root of spinal nerve|spinal ganglion|cauda equina|^(?:cervical|thoracic|lumbar|sacral|coccygeal) spinal nerves?\b/i.test(entry.name))unsegmented.add(entry.id);
 const bundled=add({id:'UNSEGMENTED',name:'Unsegmented spinal roots',description:'These source meshes bundle several roots or omit their level. Their geometry cannot be assigned to a single C1–Co1 level without segmenting the source model.',source:ROOT_MAP.sources.roots.url,matches:[...unsegmented].sort()});roots.push(bundled.id);
 const hierarchy={schema:'human-atlas-hierarchy/v2',id:'dermatomes-myotomes',name:'Dermatomes and Myotomes',version:'Spinal root map 2026-09-28 · ISNCSCI 2019',source:ROOT_MAP.sources.isncsci.url,roots,nodes};
 input.parseGuestHierarchy(hierarchy);return hierarchy;
}

function report(hierarchy,input){
 const links=createReportViewerLinks(),view=(model,id)=>{const u=new URL('http://localhost:3016/');u.searchParams.set('model',model);u.searchParams.set('tree','guest:dermatomes-myotomes');u.searchParams.set('bio',id);return u.href;};
 const lines=['# Dermatomes and Myotomes','', 'The Additional hierarchies menu (`v`) lists Genes, Cell Types, Physiology, Dermatomes and Myotomes, Drugs, Physical Exercise, then Chakras. This tree contains all 31 spinal levels (C1–C8, T1–T12, L1–L5, S1–S5, Co1), with left/right branches. Each level contains root anatomy, a dermatomal description and representative myotomal muscles.','', '## How to read the mapping','', '- C1 normally has no cutaneous dermatome. Its motor associations remain present.','- S4 and S5 remain separate root entries; the standard sensory assessment groups them as S4–5.','- Co1 retains its coccygeal plexus / anococcygeal sensory association. No isolated Co1 muscle territory is asserted.','- Dermatomal areas and muscle root contributions overlap. The hierarchy is not a diagnostic classification tool.','- The atlas has no segmented dermatome skin meshes. Sensory nodes describe their skin territories and link root anatomy; they do not select the whole skin as if it were an individual dermatome.','- Myotomal muscles are representative documented associations, not a complete census of motor fibers. ISNCSCI key actions are distinguished from anatomical multi-root contributions.','- Muscles imported as a single surface spanning several thoracic segments or both sides remain whole surfaces, explicitly labeled. Standalone tendon concepts and attachment markers are excluded. Whole source muscle surfaces can include tendon-coded materials.','- Spinal cord segments are labeled as central references; they are not substituted for peripheral roots. Unnumbered root bundles have their own branch.','- The JSON is generated offline and loaded only when this hierarchy is selected. Explode assigns each mesh once even where root associations overlap.','', '## All spinal levels','', '| Level | Sensory territory / landmark | ISNCSCI key action |','| --- | --- | --- |'];
 for(const level of SPINAL_LEVELS)lines.push(`| [${level}](${view('male-detail',`ROOT:${level}`)}) | ${markdownLabel(ROOT_MAP.sensory[level])} | ${ROOT_MAP.keyActions[level]??'No isolated key muscle'} |`);
 lines.push('','## Source geometry coverage','','Counts below are sides with a separately named root/segmental nerve; the skin column is deliberately zero because no true dermatome patches have been imported.','','| Model | Named root / segmental nerve sides | Myotome sides with muscle anatomy | Segmented skin patches |','| --- | ---: | ---: | ---: |');
 const modelFor=new Map([...links.atlases].map(([model,atlas])=>[JSON.stringify(atlas.concepts.map(c=>c.id)),model]));
 for(const {file,atlas} of input.catalogues){
  const model=modelFor.get(JSON.stringify(atlas.concepts.map(c=>c.id)))??'male-detail',resolved=input.resolveGuestHierarchy(atlas,hierarchy),all=new Map();const visit=node=>{if(all.has(node.id))return;all.set(node.id,node);node.children.forEach(visit);};resolved.forEach(visit);
  const nerves=[...all.values()].filter(node=>node.id.startsWith('NERVE:')&&node.directParts.length).length,myotomes=[...all.values()].filter(node=>node.id.startsWith('MYOTOME:')&&node.parts.length).length;
  lines.push(`| [${markdownLabel(model)}](${view(model,'ROOT:C5')}) (${file}) | ${nerves} / 62 | ${myotomes} / 62 | 0 |`);
 }
 lines.push('','## Missing root geometry','','These links open each retained root entry. A root may have muscle mappings even when its own nerve mesh is absent.');
 for(const {file,atlas} of input.catalogues){
  const model=modelFor.get(JSON.stringify(atlas.concepts.map(c=>c.id)))??'male-detail',all=new Map(),visit=node=>{if(all.has(node.id))return;all.set(node.id,node);node.children.forEach(visit);};input.resolveGuestHierarchy(atlas,hierarchy).forEach(visit);
  const missing=SPINAL_LEVELS.flatMap(level=>['left','right'].filter(side=>!all.get(`NERVE:${level}:${side}`)?.directParts.length).map(side=>`[${level} ${side}](${view(model,`NERVE:${level}:${side}`)})`));
  if(missing.length)lines.push('',`### ${model} (${file})`,'',missing.join(', '));
 }
 lines.push('','## Sources','','Credits and references are collected in [ATTRIBUTION.md](ATTRIBUTION.md#dermatomes-and-myotomes).');
 lines.push('', 'Only factual root associations and original paraphrases are used; the source exam worksheet and published illustrations are not reproduced.','', '## Rebuild and validate','', '```sh','npm run build:spinal-roots','npm run test:spinal-roots','npm run check','npm run build','```','');return lines.join('\n');
}

if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){
 const input=await loadChakraInputs(),hierarchy=buildSpinalRootHierarchy(input),outputs=[['public/assets/dermatomes-myotomes.json',JSON.stringify(hierarchy,null,2)+'\n'],['DERMATOMES-MYOTOMES.md',report(hierarchy,input)]];
 for(const [file,content] of outputs){if(process.argv.includes('--check'))assert.equal(readFileSync(file,'utf8'),content,`${file} needs rebuilding`);else writeFileSync(file,content);}
 console.log(`${SPINAL_LEVELS.length} spinal levels, 62 side branches, ${hierarchy.nodes.length} nodes; ${input.catalogues.length} catalogues mapped.`);
}
