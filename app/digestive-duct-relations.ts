import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';
import {smallBowelWallSource} from './small-bowel-relations';
import {biliaryJunctionInterface} from './biliary-junction';

export function digestiveDuctFrame(p:Part):string|undefined{
 if(p.suppressed||p.system!=='digestive')return;
 const native=nativePartFrame(p);if(native)return native;
 if(p.sourceId==='Duodenum'&&/^TA98PART:ZA:Duodenum:A05\.6\.02\.005$/.test(p.id)&&p.name==='Descending part of duodenum')return 'ZA:';
 const trunkTerm=p.id==='TA98PART:FJ3079:A05.8.01.061'?'A05.8.01.061':p.id==='TA98PART:FJ3079:A05.8.02.013'?'A05.8.02.013':undefined;
 if(trunkTerm&&p.sourceId==='FJ3079'&&p.conceptId===`TA98BILIARY:${trunkTerm}`&&!p.sourceOffset)return 'reviewed-biliary-trunk:FJ3079';
 if(p.id==='FJ3080'&&p.sourceId==='FJ3080'&&p.conceptId==='FMA14539'&&p.name==='Cystic duct'&&!p.sourceOffset)return 'reviewed-cystic-source:FJ3080';
}
const stage=(p:Part)=>{
 if(!digestiveDuctFrame(p))return;
 if(p.sourceId==='bile_duct'&&p.conceptId===p.id&&p.pairedTissue&&p.sectionSolid&&p.tissueWallAssembly==='biliary:female'&&p.sectionAssembly==='gallbladder-cystic:female:folds'&&!p.sourceOffset){
  if(p.id==='LOCAL:female:BILIARY:A05.8.01.061'&&p.name==='Common hepatic duct')return 'hepatic';
  if(p.id==='LOCAL:female:BILIARY:hepatic-tributaries'&&p.name==='Hepatic tributary ducts')return 'hepatic-tributaries';
  if(p.id==='LOCAL:female:BILIARY:distal-assembly'&&p.name==='Distal biliary and pancreatic ducts')return 'combined-distal';
 }
 if(p.id==='LOCAL:female:bile_duct:1031'&&p.conceptId==='LOCAL:female:bile_duct'&&p.sourceId==='bile_duct'&&['Bile duct','Biliary and pancreatic ducts'].includes(p.name))return 'combined';
 if(p.id==='TA98PART:FJ3079:A05.8.01.061'&&p.name==='Common hepatic duct — provisional proximal trunk')return 'hepatic';
 if(p.id==='TA98PART:FJ3079:A05.8.02.013'&&p.name==='Bile duct — provisional distal trunk')return 'bile';
 if(p.id==='FJ3080'&&p.sourceId==='FJ3080'&&p.name==='Cystic duct')return 'cystic';
 if(p.id==='LOCAL:female:cystic_duct:authored'&&p.sourceId==='cystic_duct'&&p.conceptId==='LOCAL:female:cystic_duct'&&p.name==='Cystic duct'&&!p.sourceOffset&&p.pairedTissue&&p.sectionSolid&&p.tissueWallAssembly==='cystic:female'&&p.sectionAssembly==='gallbladder-cystic:female:folds')return 'cystic';
 if(p.id==='ZA:Pancreatic duct'&&p.sourceId==='Pancreatic duct'&&p.name==='Pancreatic duct')return 'main-pancreatic';
 if(p.id==='ZA:Accessory pancreatic duct'&&p.sourceId==='Accessory pancreatic duct'&&p.name==='Accessory pancreatic duct')return 'accessory-pancreatic';
};
export function digestiveDuctSourceLimitation(p:Part):string|undefined{
 const type=stage(p);if(!type)return;
 if(p.pairedTissue&&p.sectionSolid&&['cystic:male','biliary:male','cystic:female','biliary:female'].includes(p.tissueWallAssembly??''))return p.coverageLimitation;
 const scope=type==='combined'?'The combined native source includes hepatic branches, gallbladder communication and a branch following the pancreas. It is not a separately represented common bile duct. Regional boundaries remain unverified. ':type==='bile'?'The branched bile-duct source is not a verified partition of each hepatic, cystic or common-duct region. ':'';
 return `Named duct source exterior. ${scope}Patent lumens, duct walls, sphincters and papillary openings remain unverified. Related routes use text for missing intermediates and do not certify a physical junction. The accessory pancreatic duct and common ampullary termination vary between individuals.`;
}
export function digestiveDuctDescription(p:Part):{text:string;url:string}|undefined{
 if(stage(p)==='cystic'&&p.tissueDrawGroups?.some(g=>g.neighbor?.startsWith('CYSTICFOLD:')))return {text:'The cystic duct has an open gallbladder neck and an open hepatic/bile junction. Three raised crescentic mucosal folds project from its inner wall and are separately selectable under Spiral fold.',url:'https://www.thieme-connect.com/products/ejournals/html/10.1055/s-0035-1567879'};
 if(stage(p)==='cystic'&&p.tissueWallAssembly==='cystic:female')return {text:'The finite cystic wall shares one open gallbladder neck and one open confluence with the retained native compound biliary trunk. Gallbladder communication is bidirectional. The terminal neck and course are source-guided reconstructions; calibre, cystic spiral folds, named hepatic/common bile regions and upstream/downstream organ interfaces remain under review.',url:'https://www.ncbi.nlm.nih.gov/books/NBK459246/'};
 if(stage(p)==='cystic'&&p.pairedTissue&&p.sectionSolid&&p.sectionAssembly==='gallbladder-cystic:male:folds'&&p.tissueDrawGroups?.some(g=>g.neighbor?.startsWith('TA98PART:FJ3079:')))return {text:'The cystic duct allows gallbladder filling and emptying through one common hepatic/bile confluence. Its installed finite wall shares open paired boundaries with the gallbladder and both trunk regions. Diagnostic thickness, local contour, the spiral fold and upstream/downstream organ interfaces remain under review.',url:'https://www.ncbi.nlm.nih.gov/books/NBK459246/'};
 if(stage(p)==='cystic')return {text:'The cystic duct connects the gallbladder neck to the extrahepatic biliary junction, allowing gallbladder filling and emptying. Its junction with the hepatic/bile duct varies. '+(p.sectionAssembly==='gallbladder-cystic:male:folds'?'The installed finite wall shares the gallbladder neck opening. Its diagnostic thickness, terminal hepatic/bile junction and cystic spiral folds remain under review.':p.pairedTissue&&p.sectionSolid&&p.tissueWallAssembly==='cystic:male'?'The installed finite wall has two open diagnostic ends; end shape, thickness, the spiral fold and physical joins remain under review.':'The source currently represents the exterior only; a continuous lumen, the spiral fold and physical joins remain under review.')+'',url:'https://www.ncbi.nlm.nih.gov/sites/books/NBK459288/'};
 if(stage(p)==='combined')return {text:p.tissueWallAssembly==='biliary:female'?'This retained compound wall has an open shared interface with the distinct cystic duct. Its common hepatic, common bile and pancreatic branch regions are not yet independently partitioned; tributary and papillary interfaces remain under review.':'This native assembly contains biliary branches serving the liver and gallbladder, together with a branch following the pancreas towards the duodenum. The common hepatic, cystic, common bile and pancreatic duct regions have not yet been separately partitioned.',url:'https://www.ncbi.nlm.nih.gov/books/NBK459246/'};
 if(stage(p)==='accessory-pancreatic')return {text:'The accessory pancreatic duct (duct of Santorini) is distinct from the main pancreatic duct (duct of Wirsung). When present and patent, it may drain through the minor duodenal papilla. Its termination, patency and communication with the main duct vary between individuals.',url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC4413057/'};
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** NBK459246 and PMC4413057: gallbladder communication is bidirectional;
 * the accessory pancreatic duct may end at the minor papilla or be absent. */
export function digestiveDuctRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),available=parts.filter(p=>!!digestiveDuctFrame(p));
 const add=(a:Part,b:Part,kind:'after'|'adjacent'|'continuous',note:string,via?:string[])=>{
  for(const [from,to,k] of [[a,b,kind],[b,a,kind==='after'?'before':kind]] as const){const links=graph.get(from.id)??[];links.push({kind:k,target:to,note,...(via?{via}: {})});graph.set(from.id,links);}
 };
 for(const duct of available){
  const type=stage(duct),frame=digestiveDuctFrame(duct);if(!type)continue;
  if(type==='cystic'){
   for(const trunk of available.filter(p=>['combined','combined-distal'].includes(stage(p)??'')&&digestiveDuctFrame(p)===frame))
    if(biliaryJunctionInterface(duct,trunk))add(duct,trunk,'continuous','Open shared cystic/compound trunk confluence with literal paired tissue boundaries. Gallbladder filling and emptying are bidirectional. Clinical calibre, named trunk regions and papillary interfaces remain provisional.');
   for(const trunk of available.filter(p=>duct.id==='FJ3080'&&['hepatic','bile'].includes(stage(p)??'')&&p.sourceId==='FJ3079'&&!p.sourceOffset))
    if(biliaryJunctionInterface(duct,trunk))add(duct,trunk,'continuous','One open cystic/common hepatic/common bile confluence, with literal shared wall boundaries and reciprocal tissue-only annuli. Gallbladder communication is bidirectional. Thickness and local contour are diagnostic; hepatic tributary and papillary interfaces remain unresolved.');
    else add(duct,trunk,'adjacent','Reference to the single cystic/hepatic/bile junction. The proximal and distal subdivisions belong to one reviewed trunk, not two independent outlets. Gallbladder communication is bidirectional; the source ends and physical lumen junction remain unverified.');
  }
  if(type==='combined'){
   for(const [id,source,name] of [['LOCAL:female:liver:1032','liver','Liver'],['LOCAL:female:pancreas:1039','pancreas','Pancreas']] as const){
    const organ=available.find(p=>p.id===id&&p.sourceId===source&&p.name===name&&digestiveDuctFrame(p)===frame);
    if(organ)add(organ,duct,'after',`${name} secretions drain through the corresponding branch of the combined native duct source. Branch boundaries, finite duct walls and source lumen continuity remain unverified; this route does not make the entire duct assembly a tissue component of the organ.`);
   }
  }
  if(type==='hepatic-tributaries'){
   const liver=available.find(p=>p.id==='LOCAL:female:liver:1032'&&p.sourceId==='liver'&&digestiveDuctFrame(p)===frame);
   if(liver)add(liver,duct,'after','Liver secretion enters the collective native hepatic tributary tree. Individual branches and hepatic entry interfaces remain unverified.');
   for(const hepatic of available.filter(p=>stage(p)==='hepatic'&&digestiveDuctFrame(p)===frame))if(biliaryJunctionInterface(duct,hepatic)){
    add(duct,hepatic,'after','The represented hepatic tributary wall continues into the common hepatic segment; individual tributary identities and clinical transition limits remain under review.');
    add(duct,hepatic,'continuous','Reciprocal tissue-only annuli join the two source domains without closing the lumen.');
   }
  }
  if(type==='combined-distal'){
   const pancreas=available.find(p=>p.id==='LOCAL:female:pancreas:1039'&&p.sourceId==='pancreas'&&digestiveDuctFrame(p)===frame);
   if(pancreas)add(pancreas,duct,'after','Pancreatic secretion follows the pancreatic branch within this collective native distal assembly. Named pancreatic/bile boundaries and papillary interfaces remain unverified.',['Pancreatic duct']);
  }
  if(type==='hepatic'){
   for(const bile of available.filter(p=>['bile','combined-distal'].includes(stage(p)??'')&&digestiveDuctFrame(p)===frame)){
    const joined=biliaryJunctionInterface(duct,bile);
    add(duct,bile,'after',joined?'The common hepatic lumen continues through the shared cystic confluence into the bile duct. The region boundary follows that junction; clinical extent, tributary and papillary interfaces remain provisional.':'Common hepatic and bile duct regions are provisional partitions of the same reviewed trunk. The cystic junction boundary and continuous biological lumen remain unresolved; the connection describes the named route only.');
    if(joined)add(duct,bile,'continuous','Literal shared paired boundary of the two trunk regions; assembled tissue annuli hide reciprocally and leave the lumen open.');
   }
  }
  if(['combined','combined-distal','bile','main-pancreatic','accessory-pancreatic'].includes(type)){
   const desc=available.filter(p=>p.name==='Descending part of duodenum'&&digestiveDuctFrame(p)===frame);
   const targets=desc.length?desc:available.filter(p=>{
    const source=smallBowelWallSource(p,parts)??p;
    return source.id==='LOCAL:female:intestine_duodenum:183'&&source.sourceId==='intestine_duodenum'&&source.name==='Intestine duodenum'&&digestiveDuctFrame(p)===frame;
   });
   for(const target of targets){
    const accessory=type==='accessory-pancreatic';
    const inherited=target.bowelWallSource?' Organ-level route inherited from the native whole duodenum; the descending region and opening are not partitioned in this tissue interval.':'';
    add(duct,target,'after',(accessory?'When patent, the accessory pancreatic duct may discharge at the minor duodenal papilla in the descending duodenum. Its presence, termination and relation to the main duct vary; source walls, lumens and a functioning opening remain unverified.':'The named duct discharges into the descending duodenum through the major papillary region. A common hepatopancreatic ampulla is not assumed: separate duct openings and other variants occur. Papillary tissue, sphincters and source lumen continuity remain unverified.')+inherited,[accessory?'Minor duodenal papilla':'Major duodenal papilla']);
   }
  }
 }
 const gall=available.find(p=>p.id==='LOCAL:female:gallbladder:1030'&&p.sourceId==='gallbladder'&&p.name==='Gallbladder');
 if(gall){
  const cystic=available.find(p=>p.id==='LOCAL:female:cystic_duct:authored'&&stage(p)==='cystic'&&digestiveDuctFrame(p)===digestiveDuctFrame(gall));
  if(cystic)add(gall,cystic,'adjacent','The complete gallbladder communicates through the distinct finite cystic wall for filling and emptying. Its four wall intervals share reciprocal open-neck tissue interfaces; clinical calibre and microscopic allocation remain provisional.');
  for(const bile of available.filter(p=>['combined','combined-distal'].includes(stage(p)??'')&&digestiveDuctFrame(p)===digestiveDuctFrame(gall)))add(gall,bile,'adjacent',cystic&&biliaryJunctionInterface(cystic,bile)?'Gallbladder communication follows the represented cystic wall into its open shared compound-trunk confluence. The hepatic/common bile/pancreatic regions and upstream/downstream terminal interfaces remain unverified.':'Gallbladder communicates with the biliary tract through the cystic duct for filling and emptying. This bidirectional association is not a one-way serial route. The combined source also follows hepatic and pancreatic branches; a native cystic duct is not separately partitioned. Donor ducts are not substituted and source junctions remain unverified.',['Cystic duct']);
 }
 cache.set(parts,graph);return graph;
}
