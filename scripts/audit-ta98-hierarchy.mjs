import assert from 'node:assert/strict';
import {existsSync,readdirSync,readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
import {createReportViewerLinks,markdownLabel,reportModel,reportViewUrl} from './report-viewer-links.mjs';

// Audit the same entries and navigation used by Systems, including local and
// legacy catalogues. Do not count a source folder or a spatial guess as a TA match.
const catalogues=[
 ['male-detail / male-full','public/models/atlas-male-complete.json'],
 ['male','public/models/atlas.json'],['female','public/models/atlas-hra-female.json'],
 ['embryo','public/models/atlas-embryo.json'],['cell','public/models/atlas-cell.json'],
 ['local-reference','.local-models/reference.json'],
 ['local-male','.local-models/male.json'],['local-female','.local-models/female.json'],
 ['legacy Z-Anatomy','public/models/atlas-z-anatomy.json'],
];
for(const folder of ['public/models','.local-models'])if(existsSync(folder)){
 for(const file of readdirSync(folder).filter(name=>name.endsWith('.json')).sort()){
  const path=`${folder}/${file}`;
  if(catalogues.some(([,known])=>known===path))continue;
  const data=JSON.parse(readFileSync(path,'utf8'));
  if(Array.isArray(data.parts)&&Array.isArray(data.concepts))catalogues.push([file.replace(/\.json$/,''),path]);
 }
}
const bundle=await build({entryPoints:['app/anatomy-hierarchy.ts','app/hierarchy-navigation.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/ta98-audit'});
const load=name=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(name)).text).toString('base64')}`);
const {hierarchyEntries}=await load('anatomy-hierarchy.js');
const {hierarchyNavigation}=await load('hierarchy-navigation.js');
const raw=readFileSync('app/data/ta98-metadata.json','utf8');
const metadata=JSON.parse(raw),matches=metadata.byConceptMatch;
const unresolved=new Map(),parentOnly=new Map(),outside=new Map(),suppressed=new Map();
const counts=[],unique=new Set(),verified=new Set(),vrCatalogues=[];
const root='A01.0.00.000';
const complete=id=>!!matches[id]&&metadata.byConceptPath[id]?.[0]===root;
const escape=value=>String(value??'').replaceAll('|','\\|').replaceAll('\n',' ');
const code=value=>'`'+escape(String(value).replaceAll('`',''))+'`';
const viewer=createReportViewerLinks();
const modelLink=model=>{const id=reportModel(model);return id?`[${markdownLabel(model)}](${reportViewUrl(id)})`:escape(model);};

for(const [id,path] of Object.entries(metadata.byConceptPath)){
 assert(matches[id],`${id}: path has no recorded mapping`);
 assert.equal(new Set(path).size,path.length,`${id}: duplicate ancestor`);
 for(const key of path)assert(metadata.byCode[key],`${id}: missing TA98 ancestor ${key}`);
 if(matches[id].kind==='parent'){
  assert.equal(path.at(-1),matches[id].term,`${id}: parent placement lost its containing term`);
  assert(!metadata.byConcept[id]?.ta98,`${id}: subpart has been assigned its parent's leaf ID`);
 }
}

function collect(target,entry,model,path,reason,parent){
 const key=JSON.stringify([entry.id,entry.name,entry.system,reason,parent]);
 const item=target.get(key)??{id:entry.id,name:entry.name,system:entry.system,models:new Set(),pieces:0,paths:new Set(),reason,parent};
 item.models.add(model);item.pieces+=entry.parts.length;item.paths.add(path);
 target.set(key,item);
}

for(const [model,file] of catalogues){
 if(!existsSync(file))continue;
 const atlas=JSON.parse(readFileSync(file,'utf8'));
 const available=atlas.parts.filter(part=>!part.suppressed);
 const vrFile=file.replace(/\.json$/,'.vr.json');
 if(existsSync(vrFile)){
  const vr=JSON.parse(readFileSync(vrFile,'utf8'));
  assert.deepEqual(vr.parts.map(part=>part.id).sort(),available.map(part=>part.id).sort(),`${model}: VR manifest must reuse exactly the available atlas pieces`);
  vrCatalogues.push(vrFile);
 }
 const entries=hierarchyEntries(atlas),byId=new Map(entries.map(entry=>[entry.id,entry]));
 const paths=new Map(),occurrences=new Map();
 const walk=(node,ancestors=[])=>{
  if(node.children){for(const child of node.children)walk(child,[...ancestors,node.name]);return;}
  for(const id of node.elements){paths.set(id,[...ancestors,node.name].join(' → '));occurrences.set(id,(occurrences.get(id)??0)+1);}
 };
 walk(hierarchyNavigation(atlas,'systems'));
 assert.equal(paths.size,available.length,`${model}: Systems does not account for every available piece`);
 for(const part of available)assert.equal(occurrences.get(part.id),1,`${model}: ${part.id} must occur exactly once in Systems`);
 const tally={model,file,pieces:available.length,concepts:entries.length,exact:0,parent:0,attachment:0,unresolved:0,outside:0,suppressed:atlas.parts.length-available.length};
 for(const entry of entries){
  const path=paths.get(entry.parts[0].id),n=entry.parts.length;
  unique.add(entry.id);
  if(atlas.scope==='cell'){
   tally.outside+=n;
   collect(outside,entry,model,path,'Cell-model component: TA98 macroscopic anatomy does not define this cellular hierarchy.');
  }else if(complete(entry.id)){
   verified.add(entry.id);
   if(matches[entry.id].kind==='exact')tally.exact+=n;
   else{
    tally.parent+=n;
    collect(parentOnly,entry,model,path,'Verified TA98 parent; the modeled member/subdivision has no asserted distinct TA98 leaf code.',matches[entry.id].term);
   }
  }else if(entry.taParentConcept&&complete(entry.taParentConcept)){
   tally.attachment+=n;verified.add(entry.id);
   assert.equal(entry.system,'attachments');
   assert(!entry.terminology.ta98,`${model}: attachment marker must not carry the muscle's TA98 ID`);
  }else{
   tally.unresolved+=n;
   const muscle=byId.get(entry.taParentConcept);
   const reason=matches[entry.id]
    ?'TA98 term is matched, but its reference-extract parent path is incomplete.'
    :entry.system==='attachments'
     ?muscle?`Attachment marker: parent muscle (${muscle.name}) still lacks a verified TA98 mapping.`:'Attachment marker: no matching parent muscle was found in this catalogue.'
     :/\?/.test(entry.name)?'Source label is unidentified; requires source/mesh review.'
     :atlas.scope==='embryo'?'Developmental structure: no confident TA98 match; review against embryological terminology.'
     :'No unique, reviewed TA98 term or specific parent; requires anatomical/source review.';
   collect(unresolved,entry,model,path,reason,matches[entry.id]?.term);
  }
 }
 assert.equal(tally.exact+tally.parent+tally.attachment+tally.unresolved+tally.outside,available.length,`${model}: audit accounting differs from Systems coverage`);
 for(const part of atlas.parts.filter(part=>part.suppressed)){
  collect(suppressed,{id:part.conceptId,name:part.name,system:part.system,parts:[part]},model,'Excluded by source catalogue (suppressed)','Source catalogue intentionally suppresses this mesh; it is not an available Systems leaf.');
 }
 counts.push(tally);
}

const total=key=>counts.reduce((sum,item)=>sum+item[key],0);
const lines=[
 '# TA98 Systems hierarchy audit','',
 'This report covers every piece in every catalogue present in this checkout, including the legacy catalogue and intentionally suppressed meshes. **Every available piece occurs exactly once in Systems.** Anatomical mapping remains incomplete; all unresolved cases are listed below.','',
 'Click a structure name to open it selected and focused in the local viewer. Links choose a supported model containing that concept; the Models column also opens the viewer. The legacy catalogue uses its selectable detailed-male counterpart. Suppressed source meshes link to an available equivalent with the same concept name when present; otherwise the link is marked as opening only the model.','',
 '## Scope and counting','',
 '- TA98 means **Terminologia Anatomica 1998**, as requested. TA2/TAH codes have not been substituted.',
 '- The existing Systems browsing roots and regional folders remain navigation aids. Verified anatomical containment follows TA98 ancestors; source collection folders are not appended below a known TA98 lineage.',
 '- An exact match has a TA98 term and a reference path to the human-body root. A parent placement has a verified containing TA98 term, while the modeled subdivision keeps its own atlas ID.',
 '- Origin/insertion annotation meshes inherit a verified parent muscle hierarchy and never receive the muscle’s TA98 leaf code.',
 '- Counts below are **mesh instances per catalogue**, not distinct human anatomical entities. The detailed and full male views share one catalogue, counted once. Identical concept IDs used across catalogues are grouped in the lists.',
 '- Multiple material/surface meshes sharing one catalogue concept inherit that concept’s mapping; their geometry roles are not claimed to be separate TA98 entities.',
 '- Cellular components are listed separately because TA98 does not provide the cell-model hierarchy. Suppressed meshes are listed separately because Systems deliberately excludes them.',
 `- ${vrCatalogues.length} VR geometry manifests were checked to reuse exactly their base catalogue's piece IDs; they introduce no separate anatomical identities.`,
 `- Concordance fingerprint: ${code(createHash('sha256').update(raw).digest('hex').slice(0,16))}.`,'',
 '| Model | Available pieces | Concepts | Exact TA98 | Verified parent | Muscle attachment | Needs review | Cell scope | Suppressed |',
 '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
 ...counts.map(item=>`| ${modelLink(item.model)} ([catalogue](${item.file})) | ${item.pieces} | ${item.concepts} | ${item.exact} | ${item.parent} | ${item.attachment} | ${item.unresolved} | ${item.outside} | ${item.suppressed} |`),
 `| **Total catalogue instances** | **${total('pieces')}** | **${total('concepts')}** | **${total('exact')}** | **${total('parent')}** | **${total('attachment')}** | **${total('unresolved')}** | **${total('outside')}** | **${total('suppressed')}** |`,'',
 `${unique.size} distinct available concept IDs were checked; ${verified.size} have a verified TA98 term, specific parent, or attachment-parent mapping in at least one catalogue.`,'',
 '## Repairs and sources','',
 'The concordance builder accepts conservative punctuation variants, including cranial-nerve brackets and hyphenation, and qualified names such as muscle heads. Reviewed synonym and subdivision rules are stored in [scripts/ta98-reviewed-mappings.json](scripts/ta98-reviewed-mappings.json). They distinguish hand/foot collectives, right/left bronchi and sex-specific entities. Gender-suffixed entity keys are retained for ancestry even when the published TA98 identifier is shared.','',
 'Missing hand-phalanx parents in the SQLite extract are repaired from the published [TA98 hand phalanges section](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/02.4.10%20TA98%20EN.htm).','',
 'Reference: [IFAA TA98 navigation](https://ifaa.unifr.ch/Public/EntryPage/ShowTA98.html), [reference-site scope](https://ifaa.unifr.ch/Public/EntryPage/AboutTA98New.html), and the [TA98 SQLite extract](https://github.com/mhalle/ta98-sqlite) used by the existing project. The extract reflects the IFAA reference site; reference parent gaps remain listed for review.','',
 '## Unresolved anatomical mappings','',
 `${unresolved.size} catalogue-concept rows below account for ${total('unresolved')} mesh instances. Each row includes its current Systems path so it can be found in the viewer. These structures remain available, with their atlas IDs; no anatomical identity has been guessed.`,'',
];

function table(items,includeParent=false){
 lines.push(`| Structure | Concept ID | Models | Pieces | ${includeParent?'TA98 parent/entity | ':''}Current Systems path | Review reason |`);
 lines.push(`| --- | --- | --- | ---: | ${includeParent?'--- | ':''}--- | --- |`);
 for(const item of items){
  const parent=item.parent?`${code(item.parent)} ${escape(metadata.byCode[item.parent]?.name??metadata.byConcept[item.id]?.latin??'')}`:'—';
  lines.push(`| ${viewer.link(item)} | ${code(item.id)} | ${[...item.models].map(modelLink).join(', ')} | ${item.pieces} | ${includeParent?parent+' | ':''}${[...item.paths].map(escape).join('<br>')} | ${escape(item.reason)} |`);
 }
 lines.push('');
}
const sorted=map=>[...map.values()].sort((a,b)=>a.system.localeCompare(b.system)||a.name.localeCompare(b.name)||a.id.localeCompare(b.id));
for(const system of [...new Set([...unresolved.values()].map(item=>item.system))].sort()){
 lines.push(`### ${system}`,'');table(sorted(unresolved).filter(item=>item.system===system),true);
}
lines.push('## Verified parent placements without a distinct leaf code','',
 `${parentOnly.size} catalogue-concept rows account for ${total('parent')} mesh instances. Their placement is resolved at the listed TA98 parent. A more specific leaf identity remains unasserted; use this list if a finer concordance is required.`,'');
table(sorted(parentOnly),true);
lines.push('## Cell-model components outside this TA98 hierarchy','');table(sorted(outside));
lines.push('## Suppressed source meshes','');table(sorted(suppressed));
lines.push('## Reproduce','',
 '```sh','npm run build:terminology','npm run audit:ta98','npm run test:hierarchies','```','',
 'The builder uses `/tmp/human-atlas-ta98.sqlite` when available, or downloads its configured extract. For an explicitly supplied reference database, run `python3 scripts/build-ta98-metadata.py /path/to/ta98.sqlite`. `npm run audit:ta98 -- --check` verifies that this report still matches the current catalogues and concordance.','');
const report=lines.join('\n'),file='TA98-HIERARCHY-AUDIT.md';
if(process.argv.includes('--check'))assert.equal(readFileSync(file,'utf8'),report,'TA98 audit report is stale; run npm run audit:ta98');
else writeFileSync(file,report);
console.log(`${total('pieces')} available mesh instances checked in ${counts.length} catalogues; every piece occurs exactly once in Systems.`);
console.log(`${total('exact')} exact TA98, ${total('parent')} parent placements, ${total('attachment')} verified muscle annotations, ${total('unresolved')} need review, ${total('outside')} cellular components, ${total('suppressed')} suppressed meshes.`);
console.log(`${file}: ${unresolved.size} unresolved rows, ${parentOnly.size} parent-placement rows.`);
