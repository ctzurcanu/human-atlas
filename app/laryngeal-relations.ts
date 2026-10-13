import type {Atlas,Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';

type Region='thyroid'|'cricoid'|'arytenoid-left'|'arytenoid-right'|'corniculate-left'|'corniculate-right'|'median-ligament';
// Literal source identities; a shared name does not register an HRA donor.
const male:Record<string,[Region,string,string]>={
 'Thyroid cartilage':['thyroid','Thyroid cartilage','A06.2.02.002'],
 'Cricoid cartilage':['cricoid','Cricoid cartilage','A06.2.03.001'],
 'Arytenoid cartilage.l':['arytenoid-left','Arytenoid cartilage (left)','A06.2.04.001'],
 'Arytenoid cartilage.r':['arytenoid-right','Arytenoid cartilage (right)','A06.2.04.001'],
 'Corniculate cartilage.l':['corniculate-left','Corniculate cartilage (left)','A06.2.05.001'],
 'Corniculate cartilage.r':['corniculate-right','Corniculate cartilage (right)','A06.2.05.001'],
 'Median cricothyroid ligament':['median-ligament','Median cricothyroid ligament','A06.2.03.009'],
};
const female:Record<string,[Region,string,string]>={
 'thyroid_cartilage:1079':['thyroid','Thyroid cartilage','A06.2.02.002'],
 'cricoid_cartilage:1072':['cricoid','Cricoid cartilage','A06.2.03.001'],
 'l_arytenoid_cartilage:1071':['arytenoid-left','Arytenoid cartilage (left)','A06.2.04.001'],
 'r_arytenoid_cartilage:1080':['arytenoid-right','Arytenoid cartilage (right)','A06.2.04.001'],
 'l_corniculate:1070':['corniculate-left','Corniculate (left)','A06.2.05.001'],
 'r_corniculate:1076':['corniculate-right','Corniculate (right)','A06.2.05.001'],
 'cricothyroid_ligament:1078':['median-ligament','Cricothyroid ligament','A06.2.03.009'],
};
function identity(p:Part):{frame:string;region:Region;term:string}|undefined{
 if(p.suppressed||p.sourceOffset)return;
 if(p.id.startsWith('ZA:')){
  const source=p.id.slice(3),row=male[source];
  if(row&&p.sourceId===source&&p.conceptId===p.id&&p.name===row[1])return {frame:'ZA:',region:row[0],term:row[2]};
 }
 if(p.id.startsWith('LOCAL:female:')){
  const key=p.id.slice(13),row=female[key],source=key.split(':')[0];
  if(row&&p.sourceId===source&&p.conceptId===`LOCAL:female:${source}`&&p.name===row[1])return {frame:'LOCAL:female:',region:row[0],term:row[2]};
 }
}
/** Display classification only; borrow every existing geometry/material bank. */
export function withLaryngealSystems(atlas:Atlas):Atlas{
 let changed=false;
 const parts=atlas.parts.map(p=>{
  if(!identity(p)||!['skeletal','connective','respiratory'].includes(p.system)||p.system==='respiratory')return p;
  changed=true;return {...p,system:'respiratory' as const};
 });
 return changed?{...atlas,parts}:atlas;
}
export function laryngealDescription(p:Part):{text:string;url:string}|undefined{
 const id=identity(p);if(!id||p.system!=='respiratory')return;
 const url='https://www.ncbi.nlm.nih.gov/books/NBK538202/';
 if(id.region==='thyroid')return {text:'The thyroid cartilage forms the anterior and lateral framework of the larynx. Its inferior horns articulate with the cricoid cartilage. It is cartilage, distinct from the thyroid gland.',url};
 if(id.region==='cricoid')return {text:'The cricoid cartilage forms a complete ring below the thyroid cartilage, with a narrow anterior arch and broad posterior lamina. It articulates with the thyroid and both arytenoid cartilages.',url};
 if(id.region.startsWith('arytenoid'))return {text:'The paired arytenoid cartilages articulate with the posterior cricoid. Their vocal and muscular processes participate in vocal-fold positioning; the ipsilateral corniculate lies at the apex.',url:'https://www.ncbi.nlm.nih.gov/books/NBK513252/'};
 if(id.region.startsWith('corniculate'))return {text:'The paired corniculate cartilages lie on the apices of the corresponding arytenoid cartilages, within the aryepiglottic folds.',url};
 return {text:'The median cricothyroid ligament joins the thyroid cartilage and anterior cricoid arch. This named source does not separately establish the lateral conus elasticus or its tissue attachments.',url};
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Allow these cartilage joints through the general bone/ligament filter,
 * while preserving exact frame, side, identity and availability guards. */
export function isReviewedLaryngealRelation(source:Part,target:Part,kind:ResolvedRelation['kind']):boolean{
 const a=identity(source),b=identity(target);
 if(!a||!b||a.frame!==b.frame||source.system!=='respiratory'||target.system!=='respiratory')return false;
 if(kind==='articulates')return a.region==='cricoid'&&['thyroid','arytenoid-left','arytenoid-right'].includes(b.region)||b.region==='cricoid'&&['thyroid','arytenoid-left','arytenoid-right'].includes(a.region);
 if(kind==='connects')return a.region==='median-ligament'&&['thyroid','cricoid'].includes(b.region);
 if(kind==='connectedBy')return b.region==='median-ligament'&&['thyroid','cricoid'].includes(a.region);
 return false;
}
/** NCBI NBK538202 / NBK513252: cartilage joints and median ligament.
 * Named anatomical relations, not a claim of verified mesh contact or motion. */
export function laryngealRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),byFrame=new Map<string,Map<Region,Part>>();
 for(const p of parts){const id=identity(p);if(!id||p.system!=='respiratory')continue;const rows=byFrame.get(id.frame)??new Map<Region,Part>();rows.set(id.region,p);byFrame.set(id.frame,rows);}
 const add=(rows:Map<Region,Part>,a:Region,b:Region,kind:ResolvedRelation['kind'],inverse:ResolvedRelation['kind'],note:string)=>{
  const source=rows.get(a),target=rows.get(b);if(!source||!target)return;
  for(const [from,to,k] of [[source,target,kind],[target,source,inverse]] as const){const links=graph.get(from.id)??[];links.push({kind:k,target:to,note});graph.set(from.id,links);}
 };
 for(const rows of byFrame.values()){
  add(rows,'cricoid','thyroid','articulates','articulates','Cricothyroid articulation between the thyroid inferior horns and cricoid. Joint facets, capsule and permitted motion are not separately represented by these whole cartilage sources.');
  for(const side of ['left','right'] as const){
   add(rows,'cricoid',`arytenoid-${side}`,'articulates','articulates','Ipsilateral cricoarytenoid articulation. Sliding and rotation change the vocal processes; joint surfaces, capsule and movement remain unverified in the source meshes.');
   add(rows,`corniculate-${side}`,`arytenoid-${side}`,'adjacent','adjacent','The corniculate cartilage sits on the apex of the ipsilateral arytenoid. This topographic reference does not certify their source tissue interface.');
  }
  for(const cartilage of ['thyroid','cricoid'] as const)add(rows,'median-ligament',cartilage,'connects','connectedBy','Median cricothyroid ligament attachment between thyroid and cricoid in the anterior midline. Tissue extent, entheses and the lateral conus elasticus are not established by this named source.');
 }
 cache.set(parts,graph);return graph;
}
