import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

// Source-qualified partial surfaces, not certified whole segment volumes.
// TA98 A06.5.02; NCBI NBK470197 and PMC6039811. Pulmonary arteries
// accompany segmental bronchi; veins run between segments, so they are not
// assigned one-to-one by matching a segment number.
const segments:[string,string,string[]][]=[
 ['003','Apical segmental bronchus of right lung (BI)',['Apical segmental artery of right lung']],
 ['004','Posterior segmental bronchus of right lung (BII)',['Posterior segmental artery of right lung']],
 ['005','Anterior segmental bronchus of right lung (BIII)',['Anterior segmental artery of right lung']],
 ['007','Lateral segmental bronchus of right lung (BIV)',['Lateral segmental artery of right lung']],
 ['008','Medial segmental bronchus of right lung (BV)',['Medial segmental artery of right lung']],
 ['010','Superior segmental bronchus of right lung (BVI)',['Superior segmental artery of right lung']],
 ['011','Medial basal segmental bronchus of right lung (BVII)',['Medial basal segmental artery of right lung']],
 ['012','Anterior basal segmental bronchus of right lung (BVIII)',['Anterior basal segmental artery of right lung']],
 ['013','Lateral basal segmental bronchus of right lung (BIX)',['Lateral basal segmental artery of right lung']],
 ['016','Apicoposterior segmental bronchus of left lung (BI+BII)',['Apical segmental artery of left lung','Posterior segmental artery of left lung']],
 ['017','Anterior segmental bronchus of left lung (BIII)',['Anterior segmental artery of left lung']],
 ['018','Superior lingular segmental bronchus of left lung (BIV)',['Superior lingular artery of left lung']],
 ['019','Inferior lingular segmental bronchus of left lung (BV)',['Inferior lingular artery of left lung']],
 ['021','Superior segmental bronchus of left lung (BVI)',['Superior segmental artery of left lung']],
 ['023','Anterior basal segmental bronchus of left lung (BVIII)',['Anterior basal segmental artery of left lung']],
 ['024','Lateral basal segmental bronchus of left lung (BIX)',['Lateral basal segmental artery of left lung']],
];
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Named covering association; never a certified mesh attachment or wall layer. */
export function isReviewedPulmonaryCovering(covering:Part,lung:Part):boolean{
 if(covering.suppressed||lung.suppressed||covering.system!=='respiratory'||lung.system!=='respiratory')return false;
 const frame=nativePartFrame(covering);if(!frame||nativePartFrame(lung)!==frame)return false;
 if(covering.id==='ZA:Pleura'&&covering.conceptId==='ZA:Pleura'&&covering.name==='Pleura')return /^ZA:(?:Inferior|Superior|Middle) lobe of (?:left|right) lung:Lung-(?:\d+|base)$/.test(lung.id);
 for(const [token,side] of [['l','left'],['r','right']] as const){
  if(covering.conceptId===`LOCAL:female:${token}_pleura`&&covering.id.startsWith(`${covering.conceptId}:`)&&covering.name===`Pleura (${side})`)
   return lung.conceptId===`LOCAL:female:${token}_lung`&&lung.id.startsWith(`${lung.conceptId}:`)&&lung.name===`Lung (${side})`;
 }
 return false;
}
export function pulmonarySourceLimitation(part:Part):string|undefined{
 if(part.system==='respiratory'&&/^LOCAL:female:[lr]_bronchus$/.test(part.conceptId)&&part.id.startsWith(`${part.conceptId}:`))return 'Native main-bronchus exterior from the tracheal bifurcation to the downstream branch tree. Separate wall layers, carinal anatomy and biological lumen continuity remain unverified; finer branches remain combined in an unpartitioned source.';
 if(part.system==='respiratory'&&part.conceptId==='LOCAL:female:bronchial_tree'&&part.id.startsWith(`${part.conceptId}:`))return 'Unpartitioned bilateral downstream bronchial source. Lobar and segmental branches are not separately selectable. The source retains two components, nonmanifold edges and degenerate faces; walls, lumens and fine branches remain unverified.';
 if(part.id==='ZA:Pleura'&&part.conceptId==='ZA:Pleura')return 'Undivided bilateral pleural source envelope. Separate visceral and parietal layers, fissural coverings, reflections, recesses and lung or mediastinal interfaces remain unverified.';
 if(/^LOCAL:female:[lr]_pleura$/.test(part.conceptId)&&part.id.startsWith(`${part.conceptId}:`))return 'Named pleural source surface. Separate visceral and parietal layers, fissural coverings, reflections, recesses and lung or mediastinal interfaces remain unverified.';
}
export function pulmonaryRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const previous=cache.get(parts);if(previous)return previous;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(source:Part,target:Part,kind:ResolvedRelation['kind'],reverse:ResolvedRelation['kind'],note:string)=>{
  for(const [a,b,k] of [[source,target,kind],[target,source,reverse]] as const){
   const links=graph.get(a.id)??[];if(!links.some(r=>r.kind===k&&r.target.id===b.id))links.push({kind:k,target:b,note});graph.set(a.id,links);
  }
 };
 const available=parts.filter(p=>!p.suppressed&&nativePartFrame(p));
 const coverings=available.filter(p=>p.id==='ZA:Pleura'||/^LOCAL:female:[lr]_pleura$/.test(p.conceptId));
 for(const covering of coverings)for(const lung of available){
  if(isReviewedPulmonaryCovering(covering,lung))add(covering,lung,'covers','coveredBy','Pleural covering association for the named native lung source. Separate visceral/parietal layers, reflections, fissural coverage and physical tissue interfaces remain unverified; this is not a certified mesh junction.');
 }
 const named=(source:Part,system:Part['system'],names:string[])=>available.filter(p=>p.system===system&&nativePartFrame(p)===nativePartFrame(source)&&names.some(n=>p.name.toLowerCase()===n.toLowerCase()));
 // Native female source names omit "main"; do not resolve them through
 // independently registered HRA lobar/segmental donors.
 for(const trachea of available.filter(p=>p.system==='respiratory'&&p.conceptId==='LOCAL:female:trachea'&&p.id.startsWith(`${p.conceptId}:`)&&p.name==='Trachea')){
  for(const [token,side] of [['l','left'],['r','right']] as const){
   for(const bronchus of named(trachea,'respiratory',[`Bronchus (${side})`])){
    if(bronchus.conceptId!==`LOCAL:female:${token}_bronchus`||!bronchus.id.startsWith(`${bronchus.conceptId}:`))continue;
    add(trachea,bronchus,'after','before','Tracheal bifurcation into the native ipsilateral main bronchus. Source surfaces do not verify carinal shape, wall thickness or lumen continuity; separate native female lobar and segmental representations remain missing.');
   }
  }
 }
 for(const bronchus of available.filter(p=>p.system==='respiratory'&&/^LOCAL:female:[lr]_bronchus$/.test(p.conceptId)&&p.id.startsWith(`${p.conceptId}:`)&&/^Bronchus \((left|right)\)$/.test(p.name))){
  for(const tree of named(bronchus,'respiratory',['Bronchial tree'])){
   if(tree.conceptId!=='LOCAL:female:bronchial_tree'||!tree.id.startsWith(`${tree.conceptId}:`))continue;
   add(bronchus,tree,'after','before','Downstream native bronchial source shares boundary coordinates with this main bronchus. It combines both sides: showing it includes the bilateral tree, not a separately verified ipsilateral lobar or segmental branch. Source defects, walls and lumen continuity remain unverified.');
  }
 }
 // Lobar return uses the source assembly, never a one-to-one segment vein.
 // PMC6039811; PMID17236990 documents middle-lobe return variants.
 for(const source of available){
  if(source.system!=='respiratory')continue;
  const lobe=/^ZA:(Inferior|Superior|Middle) lobe of (left|right) lung:Lung-(?:\d+|base)$/.exec(source.id);
  if(!lobe||lobe[1]==='Middle'&&lobe[2]!=='right')continue;
  const level=lobe[1]==='Inferior'?'inferior':'superior';
  const veinName=`${lobe[2]} ${level} pulmonary vein`;
  for(const vein of named(source,'venous',[veinName])){
   if(vein.id!==`ZA:${lobe[2][0].toUpperCase()+lobe[2].slice(1)} ${level} pulmonary vein`)continue;
   add(source,vein,'venous','drains',
    `Typical lobar pulmonary venous return for this native source assembly; not a one-to-one segmental vein assignment. ${lobe[1]==='Middle'?'Middle-lobe return usually enters the right superior pulmonary vein; separate atrial or inferior-vein routes also occur. ':''}Intersegmental tributaries, subject-specific drainage, ostia and physical lumen junctions remain unverified.`);
  }
 }
 for(const [suffix,bronchus,arteries] of segments){
  for(const source of available.filter(p=>p.system==='respiratory'&&p.id.startsWith('ZA:')&&p.conceptId===`TA98LUNG:A06.5.02.${suffix}`)){
   for(const target of named(source,'respiratory',[bronchus]))add(source,target,'adjacent','adjacent','Matching named segmental airway and partial source lung surface. Segment volume, airway-wall contact and bronchiolar/alveolar continuity remain unverified.');
   for(const target of named(source,'arterial',arteries))add(source,target,'arterial','supplies','Pulmonary arterial branch associated with the named segment territory; distinct from nutritive bronchial arterial supply. The partial source surface does not verify a complete perfusion volume or capillary/lumen continuity.');
  }
 }
 // The female model has whole lungs and named main bronchi, not individual
 // segment volumes. Never borrow male patches to fill that missing anatomy.
 for(const side of ['left','right']){
  for(const source of available.filter(p=>p.system==='respiratory'&&/^LOCAL:(male|female):/.test(p.id)&&p.name.toLowerCase()===`lung (${side})`)){
   for(const target of named(source,'respiratory',[`Bronchus (${side})`,`${side} main bronchus`]))add(source,target,'adjacent','adjacent','Ipsilateral main bronchus enters the lung root. The source surfaces do not establish complete hilar attachments, lobar/segmental subdivisions or lumen continuity.');
  }
 }
 cache.set(parts,graph);return graph;
}
