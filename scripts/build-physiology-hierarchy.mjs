import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {pathToFileURL} from 'node:url';
import {loadChakraInputs} from './build-chakra-hierarchy.mjs';

const read=path=>JSON.parse(readFileSync(path,'utf8'));
export const PHYSIOLOGY_MAP=read('scripts/physiology-map.json');
const unique=values=>[...new Set(values)].sort();
const slug=value=>value.toLowerCase().replace(/[^a-z0-9]+/g,'-').replace(/^-|-$/g,'');

export function buildPhysiologyHierarchy(input){
 const spec=PHYSIOLOGY_MAP,nodes=new Map(),add=node=>{assert(!nodes.has(node.id),`Duplicate ${node.id}`);nodes.set(node.id,node);return node;};
 const genes=read('public/assets/genes.json'),cells=read('public/assets/cell-types.json');
 const geneBySymbol=new Map(genes.nodes.filter(node=>node.id.startsWith('HGNC:')).map(node=>[node.name.split(' · ')[0],node]));
 const cellById=new Map(cells.nodes.map(node=>[node.id,node]));
 const evidence={version:spec.version,sources:spec.sources,mapSha256:createHash('sha256').update(readFileSync('scripts/physiology-map.json')).digest('hex'),anatomy:{},functions:spec.functions,coverage:{}};
 const used=new Set(spec.functions.flatMap(row=>row.anatomy));
 for(const key of used){
  const selector=spec.anatomy[key];assert(selector,`Unknown anatomy ${key}`);
  const pattern=new RegExp(selector.pattern,'i'),exclude=selector.exclude?new RegExp(selector.exclude,'i'):null,ids=[];
  evidence.anatomy[key]={name:selector.name,selector,models:{}};
  for(const {file,atlas} of input.catalogues){
   const names=new Map(atlas.concepts.map(concept=>[concept.id,concept.name]));
   const found=atlas.parts.filter(part=>{
    const name=names.get(part.conceptId)??part.name;
    if(key==='skin')return !part.suppressed&&input.isSkinPart(part);
    return !part.suppressed&&selector.systems.includes(part.system)&&pattern.test(name)&&!exclude?.test(name)&&!(part.system==='muscular'&&/\.[eo]\d*[lr]$/.test(name));
   }).map(part=>part.id);
   evidence.anatomy[key].models[file]=unique(found);ids.push(...found);
  }
  add({id:`ANATOMY:${key}`,name:selector.name,partIds:unique(ids)});
 }
 const roots=spec.roots.map(name=>add({id:`DOMAIN:${slug(name)}`,name,childIds:[]}).id);
 for(const row of spec.functions){
  assert(row.sources.length&&row.sources.every(key=>spec.sources[key]),`Missing evidence for ${row.id}`);
  let parent=nodes.get(`DOMAIN:${slug(row.path[0])}`);
  for(let depth=1;depth<row.path.length-1;depth++){
   const id=`GROUP:${row.path.slice(0,depth+1).map(slug).join(':')}`;
   if(!nodes.has(id)){add({id,name:row.path[depth],childIds:[]});parent.childIds.push(id);}parent=nodes.get(id);
  }
  const links=[];
  for(const id of row.cells??[]){const cell=cellById.get(id);assert(cell,`Unknown cell type ${id}`);links.push({hierarchy:'cell-types',id,name:cell.name});}
  for(const symbol of row.genes??[]){const gene=geneBySymbol.get(symbol);assert(gene,`Unknown gene ${symbol}`);links.push({hierarchy:'genes',id:gene.id,name:symbol});}
  add({id:row.id,name:row.path.at(-1),childIds:row.anatomy.map(key=>`ANATOMY:${key}`),...(links.length?{links}:{})});parent.childIds.push(row.id);
 }
 const hierarchy=input.parseGuestHierarchy({schema:'human-atlas-hierarchy/v2',id:'physiology',name:'Physiology',version:spec.version,roots,nodes:[...nodes.values()]});
 for(const {file,atlas} of input.catalogues){
  const resolved=new Map(),visit=node=>{if(resolved.has(node.id))return;resolved.set(node.id,node);node.children.forEach(visit);};input.resolveGuestHierarchy(atlas,hierarchy).forEach(visit);
  evidence.coverage[file]={mappedFunctions:spec.functions.filter(row=>resolved.get(row.id).parts.length).length,unmappedAnatomy:[...used].filter(key=>!resolved.get(`ANATOMY:${key}`).parts.length)};
 }
 evidence.statistics={domains:roots.length,functions:spec.functions.length,anatomyGroups:used.size,nodes:hierarchy.nodes.length,cellLinks:spec.functions.reduce((sum,row)=>sum+(row.cells?.length??0),0),geneLinks:spec.functions.reduce((sum,row)=>sum+(row.genes?.length??0),0)};
 return {hierarchy,evidence};
}

function report({hierarchy,evidence}){
 const viewer=id=>{const url=new URL('http://localhost:3016/');url.searchParams.set('model','local-reference');url.searchParams.set('tree','guest:physiology');url.searchParams.set('bio',id);url.searchParams.set('context','.18');return url.href;};
 const lines=['# Physiology hierarchy','','Physiology appears immediately after Cell Types. It loads only when selected.','','## Scope','',`Original Human Atlas function classification: ${evidence.statistics.domains} domains, ${evidence.statistics.functions} specific functions, and ${evidence.statistics.anatomyGroups} anatomical groups. This is a curated organ-level catalogue, not an import of the complete Gene Ontology.`,`Functions have ${evidence.statistics.cellLinks} links to Cell Types and ${evidence.statistics.geneLinks} links to Genes. These are selected mechanistic relationships; they are not an exhaustive list of every associated gene or cell.`, '', 'Microscopic processes such as glomerular filtration, surfactant secretion and hormone secretion select the corresponding modeled organ when the microscopic functional unit is absent. The organ is context, not a segmented representation of the functional unit. Names of unavailable structures remain visible without a selection checkbox; no geometry is invented. Cardiac conduction structures, bone marrow and sweat glands are among the groups whose mesh coverage varies or is missing.','','## Domains','','| Domain | Open in viewer |','| --- | --- |'];
 for(const id of hierarchy.roots)lines.push(`| ${hierarchy.nodes.find(node=>node.id===id).name} | [Open](${viewer(id)}) |`);
 lines.push('','## Geometry coverage','','| Model catalogue | Functions with geometry | Missing anatomical groups |','| --- | ---: | --- |');
 for(const [file,coverage] of Object.entries(evidence.coverage))lines.push(`| ${file} | ${coverage.mappedFunctions}/${evidence.statistics.functions} | ${coverage.unmappedAnatomy.map(key=>`[${evidence.anatomy[key].name}](${viewer('ANATOMY:'+key)})`).join(', ')||'None'} |`);
 lines.push('','## Function mappings','','| Function | Selected anatomy | Cell types | Genes |','| --- | --- | --- |');
 for(const row of evidence.functions)lines.push(`| [${row.path.join(' → ')}](${viewer(row.id)}) | ${row.anatomy.map(key=>evidence.anatomy[key].name).join(', ')} | ${(row.cells??[]).join(', ')} | ${(row.genes??[]).join(', ')} |`);
 lines.push('','## Rebuild','','`npm run build:physiology` regenerates the JSON, provenance, report and source credits. `npm run build:physiology -- --check` checks reproducibility. `npm run test:physiology` checks every model, cross-links, selection and Explode. All scientific references and licenses are in [ATTRIBUTION.md](ATTRIBUTION.md).','');
 return lines.join('\n');
}

export function writePhysiologyOutputs(result,check=false){
 const {evidence}=result;
 const credits=['<!-- physiology-credits:start -->','### Physiology','',`Human Atlas original function taxonomy and anatomy mappings, curated ${evidence.version}. ${evidence.statistics.functions} functions in ${evidence.statistics.domains} domains. Scientific references were consulted to verify structure–function facts. No textbook paragraphs, diagrams, images, tables or chapter taxonomy are reproduced in the viewer. The organ-level hierarchy and its selectors are original Human Atlas data, not the complete GO ontology or measured functional activation.`, '', 'Scientific reference: *Anatomy and Physiology 2e*, J. Gordon Betts, Peter DeSaix, Eddie Johnson, Jody E. Johnson, Oksana Korol, Dean H. Kruse, Brandon Poe, James A. Wise, Mark Womble and Kelly A. Young; OpenStax / Rice University. The current publisher [preface and terms](https://openstax.org/books/anatomy-and-physiology-2e/pages/preface) identify the textbook as [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). The book is used as a scientific reference for this independently curated classification.','',...Object.entries(evidence.sources).map(([key,source])=>`- [${source.name}](${source.url}) — reference key: ${key}.`),'','Selected Cell Ontology and HGNC cross-links point to the existing Cell Types and Genes hierarchies; their original sources and licenses are credited in the biological hierarchy section above. Node-level evidence, model coverage and the curated mapping checksum are stored in `public/assets/physiology-provenance.json`.','','<!-- physiology-credits:end -->'].join('\n');
 let attribution=readFileSync('ATTRIBUTION.md','utf8');
 attribution=attribution.includes('<!-- physiology-credits:start -->')?attribution.replace(/<!-- physiology-credits:start -->[\s\S]*?<!-- physiology-credits:end -->/,credits):attribution.replace('<!-- exercise-credits:start -->',credits+'\n\n<!-- exercise-credits:start -->');
 const outputs=[['public/assets/physiology.json',JSON.stringify(result.hierarchy)+'\n'],['public/assets/physiology-provenance.json',JSON.stringify(evidence)+'\n'],['PHYSIOLOGY-HIERARCHY.md',report(result)],['ATTRIBUTION.md',attribution],['public/ATTRIBUTION.md',attribution.replaceAll('(LICENSE)','(LICENSE.txt)').replaceAll('(LICENSES/MIT-upstream.txt)','(MIT-upstream.txt)')]];
 for(const [path,content] of outputs)if(check)assert.equal(readFileSync(path,'utf8'),content,`${path} needs rebuilding`);else writeFileSync(path,content);
}

if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href){const result=buildPhysiologyHierarchy(await loadChakraInputs());writePhysiologyOutputs(result,process.argv.includes('--check'));console.log(result.evidence.statistics);for(const [file,coverage] of Object.entries(result.evidence.coverage))console.log(file,coverage.mappedFunctions,coverage.unmappedAnatomy);}
