import {SYSTEMS,structureName,type Atlas,type Part,type SystemId} from './anatomy';
import {taHierarchyForConcept} from './anatomical-terminology';
import {anatomicalInnervationNames} from './anatomical-relations';
import {createDepthOrder} from './depth-sort';

export type GuestLink={hierarchy:string;id:string;name:string};
export type GuestNode={id?:string;name:string;synonyms?:string[];description?:string;source?:string;extension?:string;links?:GuestLink[];childIds?:string[];matchFrom?:string;matches?:string[];partIds?:string[];matchBases?:string[];systems?:SystemId[];children?:GuestNode[];unmapped?:boolean};
export type GuestHierarchy={schema:'human-atlas-hierarchy/v1'|'human-atlas-hierarchy/v2';id:string;name:string;source?:string;version?:string;revision?:number;roots?:string[];views?:{id:string;name:string;roots:string[]}[];nodes:GuestNode[]};
export type ResolvedGuestNode={id:string;name:string;synonyms?:string[];description?:string;source?:string;extension?:string;links?:GuestLink[];parts:Part[];directParts:Part[];children:ResolvedGuestNode[]};

const object=(value:unknown):value is Record<string,unknown>=>!!value&&typeof value==='object'&&!Array.isArray(value);
const strings=(value:unknown,max=100):value is string[]=>Array.isArray(value)&&value.length<=max&&value.every(item=>typeof item==='string'&&item.trim().length>0&&item.length<=200);

/** Exact anatomical names, with only the bilateral suffix removed. */
export const guestAnatomyBase=(name:string)=>name.trim().replace(/\s*\((?:left|right)\)\s*$/i,'').trim().toLocaleLowerCase();

/** Treat remote hierarchy files as data: bounded plain text and exact anatomy names. */
export function parseGuestHierarchy(value:unknown):GuestHierarchy{
 if(object(value)&&value.schema==='human-atlas-hierarchy/v2')return parseGraphHierarchy(value);
 if(!object(value)||value.schema!=='human-atlas-hierarchy/v1'||typeof value.id!=='string'||!/^[a-z0-9][a-z0-9-]{0,63}$/.test(value.id)||typeof value.name!=='string'||!value.name.trim()||value.name.length>80||!Array.isArray(value.nodes)||value.nodes.length<1||value.nodes.length>100)throw new Error('Invalid hierarchy file.');
 if(value.source!==undefined){if(typeof value.source!=='string'||value.source.length>500)throw new Error('Invalid hierarchy source.');try{if(!['http:','https:'].includes(new URL(value.source).protocol))throw new Error();}catch{throw new Error('Invalid hierarchy source.');}}
 let count=0;
 const parseNode=(input:unknown,depth:number):GuestNode=>{
  if(!object(input)||depth>8||++count>500||typeof input.name!=='string'||!input.name.trim()||input.name.length>160)throw new Error('Invalid hierarchy node.');
  if(input.matches!==undefined&&!strings(input.matches))throw new Error('Invalid anatomy matches.');
  if(input.matchBases!==undefined&&!strings(input.matchBases,2000))throw new Error('Invalid bilateral anatomy matches.');
  if(input.systems!==undefined&&!strings(input.systems,30))throw new Error('Invalid system matches.');
  if(input.children!==undefined&&(!Array.isArray(input.children)||input.children.length>100))throw new Error('Invalid hierarchy children.');
  if(input.unmapped!==undefined&&input.unmapped!==true)throw new Error('Invalid unmapped branch.');
  if(input.unmapped&&(depth!==0||input.matches||input.matchBases||input.systems||input.children))throw new Error('The unmapped branch must be a top-level empty node.');
  return {name:input.name.trim(),...(input.matches?{matches:input.matches as string[]}:{}),...(input.matchBases?{matchBases:input.matchBases as string[]}:{}),...(input.systems?{systems:input.systems as SystemId[]}:{}),...(input.children?{children:(input.children as unknown[]).map(child=>parseNode(child,depth+1))}:{}),...(input.unmapped?{unmapped:true}:{})};
 };
 const nodes=value.nodes.map(node=>parseNode(node,0));
 if(nodes.filter(node=>node.unmapped).length>1)throw new Error('Only one unmapped branch is allowed.');
 return {schema:'human-atlas-hierarchy/v1',id:value.id,name:value.name.trim(),...(value.source?{source:value.source as string}:{}),nodes};
}

const resolvedCache=new WeakMap<Atlas,WeakMap<GuestHierarchy,Map<string,ResolvedGuestNode[]>>>();
export function resolveGuestHierarchy(atlas:Atlas,hierarchy:GuestHierarchy,view=''):ResolvedGuestNode[]{
 let cache=resolvedCache.get(atlas);if(!cache){cache=new WeakMap();resolvedCache.set(atlas,cache);}
 let views=cache.get(hierarchy);if(!views){views=new Map();cache.set(hierarchy,views);}
 const cached=views.get(view);if(cached)return cached;
 const compareDepth=createDepthOrder(atlas);
 const availableParts=atlas.parts.filter(part=>!part.suppressed);
 const byName=new Map<string,Part[]>(),byBase=new Map<string,Part[]>(),byId=new Map(availableParts.map(part=>[part.id,part]));
 const concepts=new Map(atlas.concepts.map(concept=>[concept.id,concept]));
 for(const concept of atlas.concepts){const key=concept.name.toLocaleLowerCase();const members=concept.elements.map(id=>byId.get(id)).filter((part):part is Part=>!!part);byName.set(key,[...(byName.get(key)??[]),...members]);const base=guestAnatomyBase(concept.name);byBase.set(base,[...(byBase.get(base)??[]),...members]);}
 const graph=new Map(hierarchy.nodes.map(node=>[node.id,node])),memo=new Map<string,ResolvedGuestNode>(),directCache=new Map<GuestNode,Part[]>();
 const resolve=(node:GuestNode,path:string):ResolvedGuestNode=>{
  const cached=memo.get(path);if(cached)return cached;
  const children=hierarchy.schema==='human-atlas-hierarchy/v2'?(node.childIds??[]).map(id=>resolve(graph.get(id)!,id)):(node.children??[]).map((child,index)=>resolve(child,`${path}/${index}`));
  const selectors=node.matchFrom?graph.get(node.matchFrom)!:node;
  let directParts=directCache.get(selectors);
  if(!directParts){
   const matched=[...(selectors.partIds??[]).flatMap(id=>byId.has(id)?[byId.get(id)!]:[]),...(selectors.matchBases??[]).flatMap(term=>byBase.get(guestAnatomyBase(term))??[]),...(selectors.matches??[]).flatMap(term=>{
    const concept=concepts.get(term);return concept?concept.elements.map(id=>byId.get(id)).filter((part):part is Part=>!!part):byId.has(term)?[byId.get(term)!]:byName.get(term.toLocaleLowerCase())??[];
   })];
   const systems=new Set(selectors.systems??[]);
   directParts=[...new Map([...matched,...(systems.size?availableParts.filter(part=>systems.has(part.system)):[])].map(part=>[part.id,part])).values()];directCache.set(selectors,directParts);
  }
  const parts=children.every(child=>child.parts===directParts)?directParts:[...new Map([...directParts,...children.flatMap(child=>child.parts)].map(part=>[part.id,part])).values()];
  const result={id:path,name:node.name,synonyms:node.synonyms,description:node.description,source:node.source,extension:node.extension,links:node.links,parts,directParts,children};memo.set(path,result);return result;
 };
 if(hierarchy.schema==='human-atlas-hierarchy/v2'){
  const roots=hierarchy.views?.find(item=>item.id===view)?.roots??hierarchy.roots!;
  const result=roots.map(id=>resolve(graph.get(id)!,id));
  if(hierarchy.id==='dermatomes-myotomes'){
   const makeNode=(id:string,name:string,parts:Part[],children:ResolvedGuestNode[]=[],description?:string):ResolvedGuestNode=>({id,name,parts,directParts:children.length?[]:parts,children,description});
   const conceptNodes=(prefix:string,parts:Part[])=>{const grouped=new Map<string,Part[]>();for(const part of parts){const members=grouped.get(part.conceptId)??[];members.push(part);grouped.set(part.conceptId,members);}return [...grouped].map(([id,members])=>makeNode(`${prefix}:${id}`,structureName(concepts.get(id)?.name??members[0].name),members)).sort((a,b)=>a.name.localeCompare(b.name));};
   // Full catalogue access is distinct from a claim that every tissue has a known root assignment.
   const systems=SYSTEMS.flatMap(system=>{const parts=availableParts.filter(part=>part.system===system.id);return parts.length?[makeNode(`ANATOMY:${system.id}`,system.name,parts,conceptNodes(`ANATOMY:${system.id}`,parts))]:[];});
   result.push(makeNode('ANATOMY','All anatomy · by system',availableParts,systems,'Every available modeled structure. System membership does not assert a spinal root or innervation.'));
   const supply=new Map<string,Part[]>(),unresolved:Part[]=[];
   for(const part of availableParts){const names=anatomicalInnervationNames(part);if(!names.length)unresolved.push(part);for(const name of names){const members=supply.get(name)??[];members.push(part);supply.set(name,members);}}
   const nerves=[...supply].sort(([a],[b])=>a.localeCompare(b)).map(([name,parts])=>makeNode(`INNERVATION:${name}`,name,parts,conceptNodes(`INNERVATION:${name}`,parts),'Existing named innervation association. Several nerves may supply the same structure; no exclusive spinal level is inferred.'));
   const categorized=[['CRANIAL','Cranial nerve supply',/vagus|vagal|trigeminal|facial|oculomotor|trochlear|abducens|glossopharyngeal|hypoglossal|accessory|optic|olfactory|vestibulocochlear/i],['AUTONOMIC','Autonomic plexus supply',/sympathetic|splanchnic|(?:celiac|coeliac|mesenteric|hypogastric|uterovaginal|prostatic|vesical|cardiac|pulmonary|enteric|gastric) plexus/i]] as const;
   const used=new Set<string>(),groups:ResolvedGuestNode[]=[];
   for(const [id,name,pattern] of categorized){const children=nerves.filter(node=>!used.has(node.id)&&pattern.test(node.name));children.forEach(node=>used.add(node.id));if(children.length)groups.push(makeNode(`INNERVATION:${id}`,name,[...new Map(children.flatMap(node=>node.parts).map(part=>[part.id,part])).values()],children));}
   const peripheral=nerves.filter(node=>!used.has(node.id));if(peripheral.length)groups.push(makeNode('INNERVATION:PERIPHERAL','Peripheral nerve supply',[...new Map(peripheral.flatMap(node=>node.parts).map(part=>[part.id,part])).values()],peripheral));
   if(unresolved.length)groups.push(makeNode('INNERVATION:UNRESOLVED','Innervation mapping not established',unresolved,conceptNodes('INNERVATION:UNRESOLVED',unresolved),'No curated supply mapping is recorded for these source structures. This includes terms without direct tissue innervation; no spinal root is fabricated.'));
   result.push(makeNode('INNERVATION','Innervation · named nerves and unresolved anatomy',availableParts,groups));
   const nerveParts=availableParts.filter(part=>part.system==='nervous'),nerveText=new Map(nerveParts.map(part=>[part.id,[part.name,...(part.groups??[]),...taHierarchyForConcept(part.conceptId).map(entry=>entry.name)].join(' ')]));
   const cranialDefs:[string,string,RegExp][]=[['I','Olfactory',/olfactory nerve/i],['II','Optic',/optic nerve/i],['III','Oculomotor',/oculomotor nerve/i],['IV','Trochlear',/trochlear nerve/i],['V','Trigeminal',/trigeminal|ophthalmic nerve|mandibular nerve|maxillary nerve/i],['VI','Abducens',/abducens nerve/i],['VII','Facial',/facial nerve|chorda tympani/i],['VIII','Vestibulocochlear',/vestibulocochlear|cochlear nerve|vestibular nerve/i],['IX','Glossopharyngeal',/glossopharyngeal/i],['X','Vagus',/vagus|vagal|laryngeal nerve|laryngeal branch/i],['XI','Accessory',/accessory nerve/i],['XII','Hypoglossal',/hypoglossal nerve/i]];
   const cranialChildren=cranialDefs.map(([roman,name,pattern])=>{const anatomy=nerveParts.filter(part=>pattern.test(nerveText.get(part.id)!)),targets=nerves.filter(node=>pattern.test(node.name)),children=[...conceptNodes(`CRANIAL:${roman}:ANATOMY`,anatomy),...targets],members=[...new Map(children.flatMap(child=>child.parts).map(part=>[part.id,part])).values()];return makeNode(`CRANIAL:${roman}`,`${roman} · ${name} nerve`,members,children,'Named cranial nerve anatomy and existing supply relationships. Missing modeled branches remain absent; root-specific spinal innervation is not inferred.');});
   result.push(makeNode('CRANIAL','Cranial nerves',[...new Map(cranialChildren.flatMap(child=>child.parts).map(part=>[part.id,part])).values()],cranialChildren));
   const autonomousDefs:[string,string,RegExp][]=[['PARASYMPATHETIC','Parasympathetic pathways',/parasympathetic|pelvic splanchnic|vagus|vagal|oculomotor nerve|facial nerve|glossopharyngeal nerve|ciliary ganglion|pterygopalatine ganglion|submandibular ganglion|otic ganglion/i],['SYMPATHETIC','Sympathetic pathways',/\bsympathetic|(?:greater|lesser|least|thoracic|lumbar) splanchnic|sympathetic trunk|rami communicantes/i],['ENTERIC','Enteric pathways',/enteric|myenteric|submucosal (?:nerve )?plexus/i],['MIXED','Mixed autonomic plexuses',/(?:celiac|coeliac|mesenteric|hypogastric|uterovaginal|prostatic|vesical|cardiac|pulmonary|gastric|renal) (?:nerve |nervous )?plexus|hypogastric nerve|cavernous nerve/i]];
   const autonomousChildren=autonomousDefs.map(([id,name,pattern])=>{const anatomy=nerveParts.filter(part=>pattern.test(nerveText.get(part.id)!)),targets=nerves.filter(node=>pattern.test(node.name)),children=[...conceptNodes(`AUTONOMIC:${id}:ANATOMY`,anatomy),...targets],members=[...new Map(children.flatMap(child=>child.parts).map(part=>[part.id,part])).values()];return makeNode(`AUTONOMIC:${id}`,name,members,children,'Whole named source nerves and curated supply links. Mixed nerves can contain sensory and somatic fibers as well as autonomic fibers; this view does not isolate their autonomic fascicles.');});
   result.push(makeNode('AUTONOMIC','Autonomic nervous system',[...new Map(autonomousChildren.flatMap(child=>child.parts).map(part=>[part.id,part])).values()],autonomousChildren));

  }
  views.set(view,result);return result;
 }
 const mapped=hierarchy.nodes.flatMap((node,index)=>node.unmapped?[]:[resolve(node,String(index))]);
 const assigned=new Set(mapped.flatMap(node=>node.parts.map(part=>part.id)));
 const result=hierarchy.nodes.map((node,index)=>{
  if(!node.unmapped)return mapped.find(item=>item.id===String(index))!;
  const parts=availableParts.filter(part=>!assigned.has(part.id));
  const children=SYSTEMS.flatMap(system=>{
   const members=parts.filter(part=>part.system===system.id).sort((a,b)=>compareDepth({name:a.name,parts:[a]},{name:b.name,parts:[b]}));
   if(!members.length)return [];
   return [{id:`${index}/${system.id}`,name:system.name,parts:members,directParts:[],children:members.map((part,partIndex)=>({id:`${index}/${system.id}/${partIndex}`,name:structureName(part.name),parts:[part],directParts:[part],children:[]}))}];
  });
  return {id:String(index),name:node.name,parts,directParts:[],children};
 });
 views.set(view,result);return result;
}


/** Flat DAGs preserve all source parentage without expanding repeated paths. */
function parseGraphHierarchy(value:Record<string,unknown>):GuestHierarchy{
 if(typeof value.id!=='string'||!/^[a-z0-9][a-z0-9-]{0,63}$/.test(value.id)||typeof value.name!=='string'||!value.name.trim()||value.name.length>80||!Array.isArray(value.nodes)||!value.nodes.length||value.nodes.length>200_000||!strings(value.roots,10_000)||!value.roots.length)throw new Error('Invalid graph hierarchy.');
 const url=(input:unknown):string|undefined=>{if(input===undefined)return;try{if(typeof input!=='string'||input.length>500||!['http:','https:'].includes(new URL(input).protocol))throw new Error();return input;}catch{throw new Error('Invalid hierarchy source.');}};
 const nodes:GuestNode[]=value.nodes.map(input=>{
  if(!object(input)||typeof input.id!=='string'||!input.id||input.id.length>200||typeof input.name!=='string'||!input.name.trim()||input.name.length>500)throw new Error('Invalid graph node.');
  for(const key of ['childIds','matches','partIds','matchBases'])if(input[key]!==undefined&&!strings(input[key],20_000))throw new Error(`Invalid ${key}.`);
  if(input.synonyms!==undefined&&!strings(input.synonyms,300))throw new Error('Invalid synonyms.');
  if(input.matchFrom!==undefined&&(typeof input.matchFrom!=='string'||!input.matchFrom||input.matchFrom.length>200||input.matches||input.partIds||input.matchBases))throw new Error('Invalid inherited anatomy reference.');
  if(input.extension!==undefined&&(typeof input.extension!=='string'||!/^\/assets\/biology\/genes\/\d+\.json$/.test(input.extension)))throw new Error('Invalid gene extension.');
  if(input.description!==undefined&&(typeof input.description!=='string'||input.description.length>5000))throw new Error('Invalid node description.');
  if(input.links!==undefined&&(!Array.isArray(input.links)||input.links.length>5000||input.links.some(link=>!object(link)||typeof link.hierarchy!=='string'||!/^[a-z0-9][a-z0-9-]{0,63}$/.test(link.hierarchy)||!strings([link.id,link.name],2))))throw new Error('Invalid biology links.');
  return {id:input.id,name:input.name.trim(),...(input.synonyms?{synonyms:input.synonyms as string[]}:{}),...(input.matchFrom?{matchFrom:input.matchFrom as string}:{}),...(input.extension?{extension:input.extension as string}:{}),...(input.childIds?{childIds:input.childIds as string[]}:{}),...(input.matches?{matches:input.matches as string[]}:{}),...(input.partIds?{partIds:input.partIds as string[]}:{}),...(input.matchBases?{matchBases:input.matchBases as string[]}:{}),...(input.description?{description:input.description as string}:{}),...(input.links?{links:input.links as GuestLink[]}:{}),...(input.source?{source:url(input.source)}:{})};
 });
 const byId=new Map(nodes.map(node=>[node.id!,node]));if(byId.size!==nodes.length)throw new Error('Duplicate graph node ID.');
 const views:GuestHierarchy['views']=value.views===undefined?undefined:(()=>{
  if(!Array.isArray(value.views)||!value.views.length||value.views.length>10)throw new Error('Invalid hierarchy views.');
  return value.views.map(item=>{if(!object(item)||!strings([item.id,item.name],2)||!strings(item.roots,10_000)||!item.roots.length)throw new Error('Invalid hierarchy view.');return {id:item.id as string,name:item.name as string,roots:item.roots as string[]};});
 })();
 for(const node of nodes)if(node.matchFrom&&(!byId.has(node.matchFrom)||byId.get(node.matchFrom)!.matchFrom))throw new Error('Unknown or chained anatomy reference.');
 const state=new Map<string,number>();
 const visit=(id:string,depth:number)=>{if(!byId.has(id))throw new Error('Unknown graph node reference.');if(depth>64)throw new Error('Graph hierarchy is too deep.');if(state.get(id)===1)throw new Error('Hierarchy contains a cycle.');if(state.get(id)===2)return;state.set(id,1);for(const child of byId.get(id)!.childIds??[])visit(child,depth+1);state.set(id,2);};
 for(const id of [...value.roots,...(views??[]).flatMap(item=>item.roots)])visit(id,0);
 if(state.size!==nodes.length)throw new Error('Hierarchy has unreachable nodes.');
 return {schema:'human-atlas-hierarchy/v2',id:value.id,name:value.name,source:url(value.source),...(typeof value.version==='string'?{version:value.version.slice(0,200)}:{}),roots:value.roots,views,nodes};
}
