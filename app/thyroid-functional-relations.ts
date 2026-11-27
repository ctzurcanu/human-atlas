import type {Concept,Part} from './anatomy';
import type {ResolvedRelation,RelationKind} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

type Sex='male'|'female';
interface Context {sex:Sex;members:Part[];whole:Concept;links:ResolvedRelation[]}
const cache=new WeakMap<Part[],Context[]>();
const name=(p:Part)=>p.name.replace(/\s*\((?:left|right)\)$/i,'').toLowerCase();
const independent=(p:Part)=>p.suppressed||!!p.sourceOffset||!!(p.provenance&&'mirrored' in p.provenance&&p.provenance.mirrored===true);
const native=(p:Part,frame:string)=>!independent(p)&&nativePartFrame(p)===frame&&!!p.sourceId
 &&(frame==='ZA:'?p.id==='ZA:'+p.sourceId:p.id.startsWith(frame+p.sourceId+':')&&/^\d+$/.test(p.id.slice((frame+p.sourceId+':').length)));
const wholeFemale=(p:Part)=>!independent(p)&&p.system==='endocrine'&&p.name==='Thyroid gland'&&p.id==='LOCAL:female:thyroid_gland:1081'
 &&p.conceptId==='LOCAL:female:thyroid_gland'&&p.sourceId==='thyroid_gland';
const finiteTissue=(p:Part,sex:Sex)=>!independent(p)&&p.system==='endocrine'&&p.id.startsWith(`THYROID:${sex}:`)
 &&p.tissueWallAssembly===`thyroid:${sex}`&&p.sectionAssembly===`TA98:thyroid-finite:${sex}`
 &&p.pairedTissue&&p.sectionSolid&&['capsule','parenchyma','stroma'].includes(p.tissueWallLayer??'');
const maleTissue=(p:Part)=>finiteTissue(p,'male');
const tissueSex=(p:Part):Sex|undefined=>maleTissue(p)?'male':finiteTissue(p,'female')||wholeFemale(p)?'female':undefined;

/** Female donor tissues activate whole-gland routes only as the complete,
 * same-bank eleven-cell assembly. Old native vessels remain separate sources;
 * this metadata guard establishes no vessel–gland mesh junction.
 */
function completeFemaleTissues(parts:Part[]):Part[]|undefined{
 const members=parts.filter(p=>p.id.startsWith('THYROID:female:')),first=members[0];
 if(members.length!==11||new Set(members.map(p=>p.id)).size!==11||!first)return;
 if(members.some(p=>!finiteTissue(p,'female')||!p.id.startsWith(`THYROID:female:${p.tissueWallLayer}:`)
  ||p.normalSigns===undefined||p.geometryOrigin?.length!==3||!p.geometryOrigin.every(Number.isFinite)
  ||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount
  ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return;
 const regions=['right-lobe','left-lobe','isthmus','pyramidal-lobe'];
 if(members.filter(p=>p.tissueWallLayer==='capsule').length!==4
  ||members.filter(p=>p.tissueWallLayer==='parenchyma').length!==6
  ||members.filter(p=>p.tissueWallLayer==='stroma').length!==1
  ||regions.some(region=>['capsule','parenchyma'].some(role=>!members.some(p=>p.id===`THYROID:female:${role}:${region}`))))return;
 const byId=new Map(members.map(p=>[p.id,p]));
 const contactCount=(p:Part,id:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id).reduce((n,g)=>n+g.count,0);
 for(const p of members){
  if(!p.tissueDrawGroups?.length)return;
  let end=0;
  for(const group of p.tissueDrawGroups){
   if(group.start!==end||!Number.isInteger(group.count)||group.count<=0||group.count%3)return;
   end+=group.count;
   if(group.neighbor){const peer=byId.get(group.neighbor);if(!peer||peer.id===p.id||contactCount(p,peer.id)!==contactCount(peer,p.id))return;}
  }
  if(end!==p.indexCount)return;
 }
 return members;
}
const cervicalOrigins=['Superior cervical ganglion','Middle cervical ganglion','Inferior cervical ganglion'];
const nerveNames=new Set(['ganglia of sympathetic trunk','sympathetic nerves','vagus nerve (x)','vagus nerve']);
const nodeNames=new Set(['thyroid nodes','pretracheal nodes','paratracheal cervical nodes',
 'inferior deep lateral cervical nodes','superior deep lateral cervical nodes',
 'paratracheal lymph nodes','superior deep lateral cervical lymph node']);

function joinedVessel(source:Part,context:Context,parts:Part[]):boolean{
 const assembly=`thyroid-vascular:${context.sex}`;
 const eligible=(p:Part)=>!independent(p)&&p.pairedTissue&&p.sectionSolid&&p.normalSigns!==undefined
  &&p.tissueWallAssembly===assembly&&['arterial','venous'].includes(p.system)
  &&p.chunk===source.chunk&&p.positions===source.positions&&p.normals===source.normals
  &&p.vertexCount===source.vertexCount&&JSON.stringify(p.geometryOrigin)===JSON.stringify(source.geometryOrigin);
 if(!eligible(source))return false;
 const byId=new Map(parts.map(p=>[p.id,p])),gland=new Set(context.members.map(p=>p.id));
 const count=(p:Part,id:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id).reduce((n,g)=>n+g.count,0);
 const queue=[source],seen=new Set([source.id]);
 for(let i=0;i<queue.length;i++)for(const group of queue[i].tissueDrawGroups??[]){
  if(!group.neighbor)continue;const peer=byId.get(group.neighbor);if(!peer||seen.has(peer.id))continue;
  const forward=count(queue[i],peer.id);if(!forward||forward!==count(peer,queue[i].id))continue;
  if(gland.has(peer.id))return peer.chunk===source.chunk&&peer.positions===source.positions&&peer.normals===source.normals
   &&peer.vertexCount===source.vertexCount&&JSON.stringify(peer.geometryOrigin)===JSON.stringify(source.geometryOrigin);
  if(eligible(peer)){seen.add(peer.id);queue.push(peer);}
 }
 return false;
}

/** These are sourced whole-gland functional paths, not mesh-junction proofs.
 * Gray's Anatomy (Elsevier Complete Anatomy thyroid page) describes cervical
 * sympathetic origins and central/deep cervical lymph drainage. Haviarová
 * et al., IJAR2025.245, human dissection figs1–6, also shows variable vagal
 * contributions. Cervical ganglia are parallel origins, never a serial chain.
 */
function contexts(parts:Part[]):Context[]{
 const prior=cache.get(parts);if(prior)return prior;
 const result:Context[]=[];
 const male=parts.filter(maleTissue),female=completeFemaleTissues(parts);
 if(male.filter(p=>p.tissueWallLayer==='capsule').length===4&&male.filter(p=>p.tissueWallLayer==='parenchyma').length>=4)
  result.push({sex:'male',members:male,whole:{id:'ZA:Thyroid gland',name:'Thyroid gland',elements:male.map(p=>p.id),ta98Term:'A11.3.00.001',ta98Kind:'exact'},links:[]});
 if(female)result.push({sex:'female',members:female,whole:{id:'LOCAL:female:thyroid_gland',name:'Thyroid gland',elements:female.map(p=>p.id),ta98Term:'A11.3.00.001',ta98Kind:'exact'},links:[]});
 else if(!parts.some(p=>p.id.startsWith('THYROID:female:'))){
  const nativeFemale=parts.find(wholeFemale);
  if(nativeFemale)result.push({sex:'female',members:[nativeFemale],whole:{id:nativeFemale.conceptId,name:'Thyroid gland',elements:[nativeFemale.id],ta98Term:'A11.3.00.001',ta98Kind:'exact'},links:[]});
 }
 for(const context of result){
  const frame=context.sex==='male'?'ZA:':'LOCAL:female:';
  const group=(system:string,names:string[],otherFrame=frame)=>{
   const selected=parts.filter(p=>p.system===system&&names.includes(name(p))&&native(p,otherFrame));
   const byConcept=new Map<string,Part[]>();for(const p of selected){const rows=byConcept.get(p.conceptId)??[];rows.push(p);byConcept.set(p.conceptId,rows);}return [...byConcept.values()];
  };
  const add=(rows:Part[],kind:RelationKind,note:string,extra:Partial<ResolvedRelation>={})=>{
   if(!rows.length)return;
   const target=rows[0];context.links.push({kind,target,note,...extra,
    ...(rows.length>1?{targetConcept:{id:target.conceptId,name:target.name,elements:rows.map(p=>p.id)}}:{})});
  };
  let sympathetic=group('nervous',['ganglia of sympathetic trunk']);
  if(!sympathetic.length)sympathetic=group('nervous',['sympathetic nerves']);
  for(const rows of sympathetic)add(rows,'innervation','Cervical sympathetic origins; the thyroid-specific postganglionic twigs are not separately modeled.',
   {viaRoutes:cervicalOrigins.map(origin=>({via:[origin,'Periarterial thyroid plexus'],modeled:false}))});
  for(const rows of group('nervous',['vagus nerve (x)','vagus nerve']))add(rows,'innervation',
   'Variable vagal/laryngeal autonomic contribution; thyroid-specific branches are not modeled.',
   {via:['Thyroid autonomic branches'],viaModeled:false});
  for(const rows of group('lymphatic',['thyroid nodes']))add(rows,'lymphaticDrainage',
   'Perithyroidal nodes, distinct from the prelaryngeal (Delphian) nodes. Connecting thyroid lymphatic vessels are not modeled.');
  for(const rows of group('lymphatic',['pretracheal nodes']))add(rows,'lymphaticDrainage',
   'Central thyroid drainage. A possible route passes through prelaryngeal nodes; those nodes and the connecting collectors are not separately modeled.',
   {via:['Prelaryngeal lymph nodes'],viaModeled:false});
  for(const rows of group('lymphatic',['paratracheal cervical nodes']))add(rows,'lymphaticDrainage',
   'Central drainage to cervical paratracheal nodes; connecting thyroid lymphatic collectors are not modeled.');
  for(const rows of group('lymphatic',['inferior deep lateral cervical nodes']))add(rows,'lymphaticDrainage',
   'Lower deep cervical drainage, directly and from central neck nodes; connecting lymphatic collectors are not modeled.');
  for(const rows of group('lymphatic',['superior deep lateral cervical nodes','superior deep lateral cervical lymph node']))add(rows,'lymphaticDrainage',
   context.sex==='female'?'Upper/lateral deep cervical drainage. Lower deep cervical nodes and the connecting collectors are not separately modeled.':
    'Upper/lateral deep cervical drainage; connecting thyroid lymphatic collectors are not modeled.');
  for(const rows of group('lymphatic',['paratracheal lymph nodes']))add(rows,'lymphaticDrainage',
   'Grouped paratracheal source. Prelaryngeal and pretracheal nodes, cervical/thoracic subdivisions and thyroid lymphatic collectors are not separately modeled.');
  // Keep the current named supply/drainage accessible. Its physical joins
  // remain explicitly unfinished until a coupled vascular assembly replaces it.
  for(const [system,kind,names] of [['arterial','arterial',['superior thyroid artery','inferior thyroid artery']],
   ['venous','venous',['superior thyroid vein','middle thyroid vein','inferior thyroid vein']]] as const){
   const nativeRows=group(system,[...names]);
   const foundNames=new Set(nativeRows.map(rows=>name(rows[0])));
   const alternatives=context.sex==='male'?group(system,names.filter(n=>!foundNames.has(n)),'LOCAL:male:'):[];
   for(const rows of [...nativeRows,...alternatives])add(rows,kind,rows.every(p=>joinedVessel(p,context,parts))?
    'Connected to the modeled thyroid vascular network; fine branches remain partial.':
    'The thyroid–vessel mesh junction is not yet installed.');
  }
 }
 cache.set(parts,result);return result;
}

export function thyroidWholeFunctionalRelations(concept:Concept,selected:Part[],parts:Part[]):ResolvedRelation[]|undefined{
 const thyroid=concept.ta98Term==='A11.3.00.001'||['ZA:Thyroid gland','LOCAL:female:thyroid_gland'].includes(concept.id);
 if(!thyroid)return;
 const requestedSex=concept.id==='LOCAL:female:thyroid_gland'?'female':concept.id==='ZA:Thyroid gland'?'male':undefined;
 const context=contexts(parts).find(c=>(!requestedSex||c.sex===requestedSex)&&c.members.length===selected.length&&c.members.every(p=>selected.includes(p))
  &&concept.elements.length===selected.length&&concept.elements.every(id=>selected.some(p=>p.id===id)));
 return context?.links??[];
}

const networks=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function thyroidFunctionalRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=networks.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const reverse:Partial<Record<RelationKind,RelationKind>>={innervation:'innervates',lymphaticDrainage:'lymphaticTributaries',arterial:'supplies',venous:'drains'};
 for(const context of contexts(parts)){
  if(context.members.length===1)graph.set(context.members[0].id,context.links);
  for(const link of context.links)for(const id of link.targetConcept?.elements??[link.target.id]){
   const rows=graph.get(id)??[];rows.push({...link,kind:reverse[link.kind]!,target:context.members[0],targetConcept:context.whole});graph.set(id,rows);
  }
 }
 networks.set(parts,graph);return graph;
}

/** Scope reciprocal endocrine innervation to the reviewed thyroid paths. */
export function isReviewedThyroidFunctionalRelation(source:Part,target:Part,kind:RelationKind):boolean{
 const gland=['innervates','lymphaticTributaries','supplies','drains'].includes(kind)?target:source;
 const other=gland===source?target:source,sex=tissueSex(gland);if(!sex)return false;
 const frame=sex==='male'?'ZA:':'LOCAL:female:';
 const otherNative=native(other,frame)||sex==='male'&&['arterial','venous'].includes(other.system)&&native(other,'LOCAL:male:');
 if(!otherNative)return false;
 if(['innervation','innervates'].includes(kind))return other.system==='nervous'&&nerveNames.has(name(other));
 if(['lymphaticDrainage','lymphaticTributaries'].includes(kind))return other.system==='lymphatic'&&nodeNames.has(name(other));
 if(['arterial','supplies'].includes(kind))return other.system==='arterial'&&['superior thyroid artery','inferior thyroid artery'].includes(name(other));
 if(['venous','drains'].includes(kind))return other.system==='venous'&&['superior thyroid vein','middle thyroid vein','inferior thyroid vein'].includes(name(other));
 return false;
}
