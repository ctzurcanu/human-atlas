import type {Atlas,Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';

type Region='thyroid'|'cricoid'|'arytenoid-left'|'arytenoid-right'|'corniculate-left'|'corniculate-right'|'median-ligament';
type CricoidPatch=Part&{cricoidSourcePartId?:string;cricoidRegion?:'arch'|'lamina';cricoidSurface?:string};
const cricoidSuffixes=['arch','lamina-remainder','arytenoid-articular:left','arytenoid-articular:right','thyroid-articular:left','thyroid-articular:right'] as const;

function cricoidPatchIdentity(p:Part):{frame:string;region:Region;term:string}|undefined{
 const patch=p as CricoidPatch;
 if(p.suppressed||p.sourceOffset||p.system!=='respiratory'||p.conceptId!==p.id)return;
 const sex=p.id.startsWith('CRICOID:female:')?'female':p.id.startsWith('CRICOID:male:')?'male':undefined;
 if(!sex)return;
 const suffix=p.id.slice(`CRICOID:${sex}:`.length);
 if(!cricoidSuffixes.includes(suffix as typeof cricoidSuffixes[number]))return;
 const source=sex==='female'?'LOCAL:female:cricoid_cartilage:1072':'ZA:Cricoid cartilage';
 if(patch.cricoidSourcePartId!==source||p.sourceId!==(sex==='female'?'cricoid_cartilage':'Cricoid cartilage'))return;
 if(patch.cricoidRegion!==(suffix==='arch'?'arch':'lamina'))return;
 if(suffix!=='arch'&&patch.cricoidSurface!==(suffix==='lamina-remainder'?'remainder':suffix.replace(':','-')))return;
 return {frame:sex==='female'?'LOCAL:female:':'ZA:',region:'cricoid',term:'A06.2.03.001'};
}

/** Require the complete stored replacement before inheriting native joint references. */
function cricoidPatches(parts:Part[]):Part[]{
 const result:Part[]=[];
 for(const sex of ['female','male'] as const){
  const sourceId=sex==='female'?'LOCAL:female:cricoid_cartilage:1072':'ZA:Cricoid cartilage';
  const source=parts.find(p=>p.id===sourceId);
  const ids=cricoidSuffixes.map(suffix=>`CRICOID:${sex}:${suffix}`),required=new Set(ids);
  if(!source?.suppressed||source.sourceOffset||source.replacementParts?.length!==ids.length
   ||new Set(source.replacementParts).size!==ids.length||source.replacementParts.some(id=>!required.has(id)))continue;
  if(!identity({...source,suppressed:false}))continue;
  const rows=ids.map(id=>parts.filter(p=>p.id===id));
  if(rows.some(matches=>matches.length!==1))continue;
  const members=rows.map(matches=>matches[0]),first=members[0];
  if(members.some(p=>!cricoidPatchIdentity(p)||p.chunk!==first.chunk||p.positions!==first.positions
   ||p.normals!==first.normals||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))continue;
  if(parts.some(p=>!p.suppressed&&(p as CricoidPatch).cricoidSourcePartId===sourceId&&!required.has(p.id)))continue;
  const lamina=members.slice(1),laminaIds=new Set(lamina.map(p=>p.id));
  if(lamina.some(p=>p.sectionAssembly!==`TA98:cricoid:${sex}:lamina`||p.sectionSolid!==false
   ||p.sectionClosureParts?.length!==laminaIds.size||new Set(p.sectionClosureParts).size!==laminaIds.size
   ||p.sectionClosureParts.some(id=>!laminaIds.has(id))))continue;
  if(first.sectionSolid!==true||first.sectionAssembly!==`TA98:cricoid:${sex}:arch`)continue;
  result.push(...members);
 }
 return result;
}
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
 const patch=cricoidPatchIdentity(p);if(patch)return patch;
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
 const patch=p as CricoidPatch;
 if(cricoidPatchIdentity(p)){
  if(patch.cricoidRegion==='arch')return {text:'The narrow anterior arch is part of the complete cricoid cartilage ring. Its lateral partition boundaries are provisional.',url};
  if(patch.cricoidSurface?.startsWith('arytenoid-articular-'))return {text:'A superior cricoid surface beneath the ipsilateral arytenoid base. The source-based footprint remains provisional; it is a surface region, not a separate solid or joint capsule.',url};
  if(patch.cricoidSurface?.startsWith('thyroid-articular-'))return {text:'A provisional lateral/posterolateral cricoid surface region facing the thyroid inferior horn. Its present source apposition is unverified; it is not a separate solid or a joint capsule.',url};
  return {text:'The broad posterior lamina is part of the cricoid cartilage ring. Its articular surface patches belong to the same cartilage; lateral partition boundaries remain provisional.',url};
 }
 if(id.region==='thyroid')return {text:'The thyroid cartilage forms the anterior and lateral framework of the larynx. Its inferior horns articulate with the cricoid cartilage. It is cartilage, distinct from the thyroid gland.',url};
 if(id.region==='cricoid')return {text:'The cricoid cartilage forms a complete ring below the thyroid cartilage, with a narrow anterior arch and broad posterior lamina. It articulates with the thyroid and both arytenoid cartilages.',url};
 if(id.region.startsWith('arytenoid'))return {text:'The paired arytenoid cartilages articulate with the posterior cricoid. Their vocal and muscular processes participate in vocal-fold positioning; the ipsilateral corniculate lies at the apex.',url:'https://www.ncbi.nlm.nih.gov/books/NBK513252/'};
 if(id.region.startsWith('corniculate'))return {text:'The paired corniculate cartilages lie on the apices of the corresponding arytenoid cartilages, within the aryepiglottic folds.',url};
 return {text:'The median cricothyroid ligament joins the thyroid cartilage and anterior cricoid arch. This named source does not separately establish the lateral conus elasticus or its tissue attachments.',url};
}

/** Describe a complete stored selection as a whole, rather than its first patch. */
export function laryngealSelectionDescription(selected:Part[],allParts:Part[]):{text:string;url:string}|undefined{
 if(selected.length<2||selected.some(p=>!cricoidPatchIdentity(p)))return;
 const registered=cricoidPatches(allParts),ids=new Set(selected.map(p=>p.id));
 if(ids.size!==selected.length||selected.some(p=>!registered.includes(p)))return;
 const sex=selected[0].id.startsWith('CRICOID:female:')?'female':'male';
 const members=cricoidSuffixes.map(s=>`CRICOID:${sex}:${s}`);
 const same=(expected:string[])=>expected.length===ids.size&&expected.every(id=>ids.has(id));
 const url='https://www.ncbi.nlm.nih.gov/books/NBK538202/';
 if(same(members))return {text:'The cricoid cartilage forms a complete ring below the thyroid cartilage, with a narrow anterior arch and broad posterior lamina. It articulates with the thyroid and both arytenoid cartilages.',url};
 if(same(members.slice(1)))return {text:'The broad posterior lamina is part of the cricoid cartilage ring. This selection includes its remainder and all four articular surface patches. Regional boundaries and joint apposition remain provisional.',url};
 if(same(members.slice(2,4)))return {text:'The paired superior cricoid surface regions lie beneath the arytenoid bases. They are open surface patches, not separate solids or joint capsules; their footprints and joint apposition remain provisional.',url};
 if(same(members.slice(4,6)))return {text:'The paired lateral cricoid surface regions face the inferior thyroid horns. They are open surface patches, not separate solids or joint capsules; their footprints and joint apposition remain provisional.',url};
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Allow these cartilage joints through the general bone/ligament filter,
 * while preserving exact frame, side, identity and availability guards. */
export function isReviewedLaryngealRelation(source:Part,target:Part,kind:ResolvedRelation['kind']):boolean{
 const a=identity(source),b=identity(target);
 if(!a||!b||a.frame!==b.frame||source.system!=='respiratory'||target.system!=='respiratory')return false;
 const permits=(cartilage:Part,other:Region)=>{
  if(!cricoidPatchIdentity(cartilage))return true;
  const p=cartilage as CricoidPatch;
  return other==='median-ligament'?p.cricoidRegion==='arch'
   :other==='thyroid'?p.cricoidSurface?.startsWith('thyroid-articular-')===true
   :other==='arytenoid-left'?p.cricoidSurface==='arytenoid-articular-left'
   :other==='arytenoid-right'?p.cricoidSurface==='arytenoid-articular-right':false;
 };
 if(kind==='articulates')return a.region==='cricoid'&&['thyroid','arytenoid-left','arytenoid-right'].includes(b.region)&&permits(source,b.region)
  ||b.region==='cricoid'&&['thyroid','arytenoid-left','arytenoid-right'].includes(a.region)&&permits(target,a.region);
 if(kind==='connects')return a.region==='median-ligament'&&['thyroid','cricoid'].includes(b.region)&&permits(target,a.region);
 if(kind==='connectedBy')return b.region==='median-ligament'&&['thyroid','cricoid'].includes(a.region)&&permits(source,b.region);
 return false;
}
/** NCBI NBK538202 / NBK513252: cartilage joints and median ligament.
 * Named anatomical relations, not a claim of verified mesh contact or motion. */
export function laryngealRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),byFrame=new Map<string,Map<Region,Part>>();
 for(const p of parts){const id=identity(p);if(!id||p.system!=='respiratory'||cricoidPatchIdentity(p))continue;const rows=byFrame.get(id.frame)??new Map<Region,Part>();rows.set(id.region,p);byFrame.set(id.frame,rows);}
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
 for(const patch of cricoidPatches(parts)){
  const id=cricoidPatchIdentity(patch)!,rows=byFrame.get(id.frame);if(!rows)continue;
  for(const [region,peer] of rows){
   const kind:ResolvedRelation['kind']=region==='median-ligament'?'connectedBy':'articulates';
   if(!isReviewedLaryngealRelation(patch,peer,kind))continue;
   const note=kind==='connectedBy'?'Median cricothyroid ligament reference at the anterior arch; tissue attachment remains under review.'
    :'Source-qualified articular-surface reference; the footprint is provisional and does not certify a closed joint capsule or permitted motion.';
   const forward=graph.get(patch.id)??[];forward.push({kind,target:peer,note});graph.set(patch.id,forward);
   const back=graph.get(peer.id)??[];back.push({kind:kind==='connectedBy'?'connects':kind,target:patch,note});graph.set(peer.id,back);
  }
 }
 cache.set(parts,graph);return graph;
}
