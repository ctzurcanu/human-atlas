import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

export function digestiveDuctFrame(p:Part):string|undefined{
 if(p.suppressed||p.system!=='digestive')return;
 const native=nativePartFrame(p);if(native)return native;
 if(p.sourceId==='Duodenum'&&/^TA98PART:ZA:Duodenum:A05\.6\.02\.005$/.test(p.id)&&p.name==='Descending part of duodenum')return 'ZA:';
 if(p.sourceId==='FJ3079'&&/^TA98PART:FJ3079:A05\.8\.(?:01\.061|02\.013)$/.test(p.id))return 'reviewed-biliary-trunk:FJ3079';
}
const stage=(p:Part)=>{
 if(!digestiveDuctFrame(p))return;
 if(p.id==='LOCAL:female:bile_duct:1031'&&p.sourceId==='bile_duct'&&p.name==='Bile duct')return 'bile';
 if(p.id==='TA98PART:FJ3079:A05.8.01.061'&&p.name==='Common hepatic duct — provisional proximal trunk')return 'hepatic';
 if(p.id==='TA98PART:FJ3079:A05.8.02.013'&&p.name==='Bile duct — provisional distal trunk')return 'bile';
 if(p.id==='ZA:Pancreatic duct'&&p.sourceId==='Pancreatic duct'&&p.name==='Pancreatic duct')return 'main-pancreatic';
 if(p.id==='ZA:Accessory pancreatic duct'&&p.sourceId==='Accessory pancreatic duct'&&p.name==='Accessory pancreatic duct')return 'accessory-pancreatic';
};
export function digestiveDuctSourceLimitation(p:Part):string|undefined{
 const type=stage(p);if(!type)return;
 const scope=type==='bile'?'The branched bile-duct source is not a verified partition of each hepatic, cystic or common-duct region. ':'';
 return `Named duct source exterior. ${scope}Patent lumens, duct walls, sphincters and papillary openings remain unverified. Related routes use text for missing intermediates and do not certify a physical junction. The accessory pancreatic duct and common ampullary termination vary between individuals.`;
}
export function digestiveDuctDescription(p:Part):{text:string;url:string}|undefined{
 if(stage(p)==='accessory-pancreatic')return {text:'The accessory pancreatic duct (duct of Santorini) is distinct from the main pancreatic duct (duct of Wirsung). When present and patent, it may drain through the minor duodenal papilla. Its termination, patency and communication with the main duct vary between individuals.',url:'https://pmc.ncbi.nlm.nih.gov/articles/PMC4413057/'};
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** NBK459246 and PMC4413057: gallbladder communication is bidirectional;
 * the accessory pancreatic duct may end at the minor papilla or be absent. */
export function digestiveDuctRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),available=parts.filter(p=>!!digestiveDuctFrame(p));
 const add=(a:Part,b:Part,kind:'after'|'adjacent',note:string,via?:string[])=>{
  for(const [from,to,k] of [[a,b,kind],[b,a,kind==='after'?'before':'adjacent']] as const){const links=graph.get(from.id)??[];links.push({kind:k,target:to,note,...(via?{via}: {})});graph.set(from.id,links);}
 };
 for(const duct of available){
  const type=stage(duct),frame=digestiveDuctFrame(duct);if(!type)continue;
  if(type==='hepatic'){
   for(const bile of available.filter(p=>stage(p)==='bile'&&digestiveDuctFrame(p)===frame))add(duct,bile,'after','Common hepatic and bile duct regions are provisional partitions of the same reviewed trunk. The cystic junction boundary and continuous biological lumen remain unresolved; the connection describes the named route only.');
  }
  if(['bile','main-pancreatic','accessory-pancreatic'].includes(type)){
   const desc=available.filter(p=>p.name==='Descending part of duodenum'&&digestiveDuctFrame(p)===frame);
   const targets=desc.length?desc:available.filter(p=>p.id==='LOCAL:female:intestine_duodenum:183'&&p.sourceId==='intestine_duodenum'&&p.name==='Intestine duodenum'&&digestiveDuctFrame(p)===frame);
   for(const target of targets){
    const accessory=type==='accessory-pancreatic';
    add(duct,target,'after',accessory?'When patent, the accessory pancreatic duct may discharge at the minor duodenal papilla in the descending duodenum. Its presence, termination and relation to the main duct vary; source walls, lumens and a functioning opening remain unverified.':'The named duct discharges into the descending duodenum through the major papillary region. A common hepatopancreatic ampulla is not assumed: separate duct openings and other variants occur. Papillary tissue, sphincters and source lumen continuity remain unverified.',[accessory?'Minor duodenal papilla':'Major duodenal papilla']);
   }
  }
 }
 const gall=available.find(p=>p.id==='LOCAL:female:gallbladder:1030'&&p.sourceId==='gallbladder'&&p.name==='Gallbladder');
 if(gall)for(const bile of available.filter(p=>stage(p)==='bile'&&digestiveDuctFrame(p)===digestiveDuctFrame(gall)))add(gall,bile,'adjacent','Gallbladder communicates with the biliary tract through the cystic duct for filling and emptying. This bidirectional association is not a one-way serial route. A native cystic duct is not separately represented; donor ducts are not substituted and source junctions remain unverified.',['Cystic duct']);
 cache.set(parts,graph);return graph;
}
