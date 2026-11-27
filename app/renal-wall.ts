import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';
import type {ResolvedRelation} from './anatomical-relations';
import {urinaryExternalPeer} from './urinary-wall';

type Source={sourcePartId?:string;wholeConceptId:string;partIds:string[]};
type Assembly={side:'left'|'right';partIds:string[];sources:Source[]};
type Marker={schema:string;sex:string;partIds:string[];assemblies:Assembly[]};
type Registry={marker:Marker;assemblies:{definition:Assembly;parts:Part[]}[]};
type RenalPart=Part&{renalSourcePartId?:string;renalRegion?:string;renalRepresentation?:string;renalParentKidneyPartId?:string;renalWholeConceptId?:string};
const registries=new WeakMap<Atlas,Registry|null>();
const same=(a:readonly string[],b:readonly string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));
const regions=['capsule','renal-tissue','pelvis','ureter-abdominal','ureter-pelvic','cortex','pyramid','papilla','upj','ureter'];
const coats=['adventitia','muscularis','mucosa'];
const detailed=(p:RenalPart)=>['renal-detail-v416','renal-layered-v1'].includes(p.renalRepresentation??'');
const relationCaches=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();

/** Only stored reciprocal wall interfaces establish these continuations.
 * Neither spatial proximity nor an authority parent creates a tube join. */
export function renalTissueRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=relationCaches.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),byId=new Map(parts.map(p=>[p.id,p]));
 for(const p of parts as RenalPart[]){
  const order=detailed(p)?['pelvis','upj','ureter']:['pelvis','ureter-abdominal','ureter-pelvic'];
  const from=order.indexOf(p.renalRegion??'');
  if(from<0||p.suppressed||p.system!=='urinary'||!p.pairedTissue||!p.tissueWallAssembly?.startsWith('renal:'))continue;
  for(const group of p.tissueDrawGroups??[]){
   const peer=group.neighbor?byId.get(group.neighbor)as RenalPart|undefined:undefined,to=order.indexOf(peer?.renalRegion??'');
   if(!peer||to<0||Math.abs(to-from)!==1||peer.suppressed||peer.system!=='urinary'||!peer.pairedTissue
    ||peer.tissueWallAssembly!==p.tissueWallAssembly||peer.tissueWallLayer!==p.tissueWallLayer
    ||peer.chunk!==p.chunk||peer.positions!==p.positions||peer.normals!==p.normals||peer.vertexCount!==p.vertexCount
    ||JSON.stringify(peer.geometryOrigin)!==JSON.stringify(p.geometryOrigin)
    ||peer.tissueDrawGroups?.filter(g=>g.neighbor===p.id).reduce((n,g)=>n+g.count,0)!==group.count)continue;
   const links=graph.get(p.id)??[];
   if(!links.some(r=>r.target.id===peer.id))links.push({kind:to>from?'after':'before',target:peer,note:'Stored reciprocal wall annulus; the collecting lumen remains open. The distal bladder junction remains unfinished.'});
   graph.set(p.id,links);
  }
 }
 relationCaches.set(parts,graph);return graph;
}

/** A complete stored replacement retains the original semantic whole. A coat
 * alone must never inherit all of that organ's supplying vessels or nerves. */
export function renalWholeSource(selected:Part[],parts:Part[]):Part|undefined{
 if(!selected.length)return;
 const first=selected[0]as RenalPart,id=first.renalSourcePartId;
 if(!id){
  // The female source had no separate native pelvis mesh. Its new registered
  // whole is a semantic query proxy over two actual joined wall cells, never
  // an invented rendered source mesh or a reference in an unrelated frame.
  const parent=parts.find(p=>p.id===first.renalParentKidneyPartId),members=selected as RenalPart[];
  const layered=first.renalRepresentation==='renal-layered-v1',count=layered?3:1;
  if(members.length!==2*count||!parent?.suppressed||!/^LOCAL:female:[lr]_kidney:\d+$/.test(parent.id)||parent.replacementParts?.length!==24
   ||!first.renalWholeConceptId||['pelvis','upj'].some(region=>members.filter(p=>p.renalRegion===region).length!==count)
   ||layered&&['pelvis','upj'].some(region=>!same(members.filter(p=>p.renalRegion===region).map(p=>p.tissueWallLayer??''),coats))
   ||members.some(p=>p.suppressed||p.system!=='urinary'||p.renalSourcePartId||!detailed(p)||p.renalRepresentation!==first.renalRepresentation||p.renalParentKidneyPartId!==parent.id||p.renalWholeConceptId!==first.renalWholeConceptId
    ||!p.pairedTissue||!p.sectionSolid||p.sourceOffset||p.normalSigns===undefined||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount||p.tissueWallAssembly!==first.tissueWallAssembly||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin))
   ||members.some(p=>!p.tissueDrawGroups?.some(g=>members.some(q=>q.id===g.neighbor&&q.renalRegion!==p.renalRegion&&q.tissueWallLayer===p.tissueWallLayer))))return;
  const side=parent.id.startsWith('LOCAL:female:l_')?'left':'right';
  return {...parent,id:`LOCAL:female:registered_renal_pelvis:${side}`,conceptId:first.renalWholeConceptId,name:`Renal pelvis (${side})`,replacementParts:members.map(p=>p.id)};
 }
 if(selected.some(p=>(p as RenalPart).renalSourcePartId!==id))return;
 const source=parts.find(p=>p.id===id);
 if(!source?.suppressed||source.system!=='urinary'||!same(source.replacementParts??[],selected.map(p=>p.id)))return;
 if(!/^(?:ZA:(?:Kidney|Renal pelvis|Ureter)\.[lr]|LOCAL:female:[lr]_(?:kidney|ureter):\d+)$/.test(id))return;
 if(selected.some(p=>p.suppressed||!p.pairedTissue||!p.sectionSolid||p.system!=='urinary'||p.sourceOffset
  ||p.normalSigns===undefined||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals
  ||p.vertexCount!==first.vertexCount||p.tissueWallAssembly!==first.tissueWallAssembly
  ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return;
 const roles=selected.map(p=>p.tissueWallLayer??'');
 if(detailed(first)){
  if(selected.some(p=>(p as RenalPart).renalRepresentation!==first.renalRepresentation))return;
  if(/kidney/i.test(source.name)&&!(roles.length===24&&roles.filter(r=>r==='capsule').length===1&&roles.filter(r=>r==='cortex').length===1&&roles.filter(r=>r==='pyramid').length===11&&roles.filter(r=>r==='papilla').length===11))return;
  const layered=first.renalRepresentation==='renal-layered-v1';
  if(/renal pelvis/i.test(source.name)&&!(layered?roles.length===6&&coats.every(r=>roles.filter(q=>q===r).length===2):roles.length===2&&roles.every(r=>r==='wall')))return;
  if(/ureter/i.test(source.name)&&!(layered?same(roles,coats):roles.length===1&&roles[0]==='wall'))return;
  return source;
 }
 if(/kidney/i.test(source.name)&&!same(roles,['capsule','renal-tissue']))return;
 if(/renal pelvis/i.test(source.name)&&!same(roles,coats))return;
 if(/ureter/i.test(source.name)&&(selected.length!==6||coats.some(role=>roles.filter(r=>r===role).length!==2)))return;
 return source;
}

export function renalWallRegistry(atlas:Atlas):Registry|null{
 if(registries.has(atlas))return registries.get(atlas)!;
 const reject=()=>{registries.set(atlas,null);return null;};
 const marker=(atlas as Atlas&{finiteRenalIntegration?:Marker}).finiteRenalIntegration;
 if(!marker||!['renal-local-tissues-v281','renal-local-detail-v420','renal-local-layers-v1','renal-local-layers-v2'].includes(marker.schema)||marker.sex!==atlas.sex||!['male','female'].includes(marker.sex)
  ||!Array.isArray(marker.assemblies)||!same(marker.assemblies.map(g=>g.side),['left','right'])
  ||!Array.isArray(marker.partIds)||!same(marker.partIds,marker.assemblies.flatMap(g=>g.partIds??[])))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),assemblies=[];
 for(const g of marker.assemblies){
  const detailed=marker.schema!=='renal-local-tissues-v281';
  const separateUreters=marker.schema==='renal-local-layers-v2';
  if(!Array.isArray(g.partIds)||g.partIds.length!==(separateUreters?30:marker.schema==='renal-local-layers-v1'?33:detailed?27:11)||!Array.isArray(g.sources)||g.sources.length!==(separateUreters?2:3))return reject();
  const members=g.partIds.map(id=>byId.get(id)),first=members[0];
  if(!first||members.some(p=>!p||p.suppressed||p.system!=='urinary'||!p.pairedTissue||!p.sectionSolid||!p.tissueAppearance
   ||p.sourceOffset||p.normalSigns===undefined||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals
   ||p.vertexCount!==first.vertexCount||p.tissueWallAssembly!==`renal:${marker.sex}:${g.side}`
   ||!regions.includes((p as RenalPart).renalRegion??'')||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return reject();
  const parts=members as Part[];
  for(const s of g.sources){
   const p=s.partIds.map(id=>byId.get(id));
   if(p.some(q=>!q))return reject();
   if(s.sourcePartId){if(renalWholeSource(p as Part[],atlas.parts)?.id!==s.sourcePartId)return reject();}
   else if(!detailed||marker.sex!=='female'||renalWholeSource(p as Part[],atlas.parts)?.conceptId!==s.wholeConceptId)return reject();
   const whole=atlas.concepts.find(c=>c.id===s.wholeConceptId);
   if(!whole||!same(whole.elements,s.partIds))return reject();
  }
  if(!same(g.sources.flatMap(s=>s.partIds),g.partIds))return reject();
  for(const p of parts){
   let end=0;
   for(const group of p.tissueDrawGroups??[]){
    if(group.start!==end||!Number.isInteger(group.count)||group.count<=0||group.count%3)return reject();
    end+=group.count;
    if(group.neighbor){const peer=parts.find(q=>q.id===group.neighbor);
     const adjoining=peer??urinaryExternalPeer(atlas,p,group.neighbor);
     if(!adjoining||adjoining.id===p.id||adjoining.tissueDrawGroups?.filter(q=>q.neighbor===p.id).reduce((n,q)=>n+q.count,0)!==group.count)return reject();}
   }
   if(end!==p.indexCount)return reject();
  }
  assemblies.push({definition:g,parts});
 }
 const value={marker,assemblies};registries.set(atlas,value);return value;
}

export function renalTissuePeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const registry=renalWallRegistry(atlas);if(!registry||choice.id==='selection-set'||!choice.elements.length)return null;
 const ids=new Set(choice.elements),parts=registry.assemblies.flatMap(g=>g.parts).filter(p=>ids.has(p.id));
 if(parts.length!==ids.size)return null;
 if(registry.marker.schema!=='renal-local-tissues-v281'){
  const wholes=registry.assemblies.flatMap(g=>g.definition.sources).filter(s=>s.partIds.length===24&&s.partIds.some(id=>ids.has(id)));
  if(wholes.length){
   if(wholes.some(s=>!s.partIds.every(id=>ids.has(id)))||wholes.reduce((n,s)=>n+s.partIds.length,0)!==parts.length)return null;
   return [{name:'Fibrous capsule',ids:parts.filter(p=>p.tissueWallLayer==='capsule').map(p=>p.id)},
   {name:'Peripheral cortex',ids:parts.filter(p=>p.tissueWallLayer==='cortex').map(p=>p.id)},
   {name:'Medullary pyramids and papillae',ids:parts.filter(p=>['pyramid','papilla'].includes(p.tissueWallLayer??'')).map(p=>p.id)}];
  }
  if(!['renal-local-layers-v1','renal-local-layers-v2'].includes(registry.marker.schema))return null;
 }
 const keys=new Set(parts.map(p=>(p as RenalPart).renalRegion!));
 const kidney=[...keys].every(k=>k==='capsule'||k==='renal-tissue');
 if(kidney){
  const sides=registry.assemblies.filter(g=>g.parts.filter(p=>['capsule','renal-tissue'].includes((p as RenalPart).renalRegion!)).every(p=>ids.has(p.id)));
  if(sides.length*2!==parts.length)return null;
  return ['capsule','renal-tissue'].map(role=>({name:role==='capsule'?'Fibrous capsule':'Renal tissue (unpartitioned)',ids:parts.filter(p=>p.tissueWallLayer===role).map(p=>p.id)}));
 }
 const wallRegions=registry.marker.schema==='renal-local-layers-v2'?['pelvis','upj']:registry.marker.schema==='renal-local-layers-v1'?['pelvis','upj','ureter']:['pelvis','ureter-abdominal','ureter-pelvic'];
 if([...keys].some(k=>!wallRegions.includes(k)))return null;
 for(const g of registry.assemblies)for(const region of wallRegions){
  const members=g.parts.filter(p=>(p as RenalPart).renalRegion===region),count=members.filter(p=>ids.has(p.id)).length;
  if(count&&count!==3)return null;
 }
 return coats.map(role=>({name:role==='adventitia'?'Adventitia':role==='muscularis'?'Muscular layer':'Mucosa',ids:parts.filter(p=>p.tissueWallLayer===role).map(p=>p.id)}));
}
