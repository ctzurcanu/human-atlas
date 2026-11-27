import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';

const nativeIdentity=/^LOCAL:(female|male):(CAECOWALL:caecum|BOWELWALL:ileum):(serosa|subserosa|longitudinal-muscle|circular-muscle|submucosa|muscularis-mucosae|lamina-propria|epithelium)$/;
function nativeWall(p:Part){
 const match=p.id.match(nativeIdentity);
 const role=match?.[2]==='CAECOWALL:caecum'?'caecum':'ileum';
 return match&&!p.suppressed&&p.system==='digestive'&&p.pairedTissue&&p.sectionSolid&&!p.sourceOffset&&
  p.ileocaecalJunctionAssembly===`ileocaecal:${match[1]}:v138`&&p.ileocaecalJunctionRole===role&&p.ileocaecalJunctionLayer===match[3]&&
  p.geometryOrigin?.length===3&&p.geometryOrigin.every(Number.isFinite)&&p.tissueDrawGroups?.length?match:null;
}
/** A stored tissue annulus around a real opening, never proximity or a plug. */
export function ileocaecalWallInterface(source:Part,target:Part):boolean{
 const a=nativeWall(source),b=nativeWall(target);
 return !!a&&!!b&&a[1]===b[1]&&a[2]!==b[2]&&a[3]===b[3]&&source.chunk===target.chunk&&source.positions===target.positions&&source.normals===target.normals&&source.vertexCount===target.vertexCount&&
  JSON.stringify(source.geometryOrigin)===JSON.stringify(target.geometryOrigin)&&source.tissueDrawGroups!.some(g=>g.neighbor===target.id&&g.count>0&&target.tissueDrawGroups!.some(h=>h.neighbor===source.id&&h.count===g.count));
}
const nativeCache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function ileocaecalWallRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=nativeCache.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),byId=new Map(parts.map(p=>[p.id,p]));
 for(const source of parts){
  if(!nativeWall(source))continue;
  const linked=new Set<string>();
  for(const g of source.tissueDrawGroups!){
   const target=g.neighbor?byId.get(g.neighbor):undefined;
   if(!target||linked.has(target.id)||!ileocaecalWallInterface(source,target))continue;
   linked.add(target.id);const relations=graph.get(source.id)??[];
   relations.push({kind:'continuous',target,note:'Matching ileal and caecal wall coats meet around the open ileal orifice. The assembled annular interface stays hidden; separate coats retain their tissue rims.'});graph.set(source.id,relations);
  }
 }
 nativeCache.set(parts,graph);return graph;
}

/** Original united-female reference identities; these are not native-body
 * attachments or an accepted finite valve assembly. */
const identities:Record<string,[string,string]>={
 'HRA:VH_F_ileocecal_valve':['VH_F_ileocecal_valve','Ileocecal valve'],
 'HRA:VH_F_caecum':['VH_F_caecum','Caecum'],
 'HRA:VH_F_ileum_terminal':['VH_F_ileum_terminal','Distal part of ileum'],
 'HRA:VH_F_vermiform_appendix':['VH_F_vermiform_appendix','Vermiform appendix'],
};
function reference(p:Part):boolean {
 const identity=identities[p.id];
 return !!identity&&!p.suppressed&&p.system==='digestive'&&p.sourceId===identity[0]&&p.name===identity[1]&&
  (p as Part&{reviewSource?:string}).reviewSource==='female-HRA'&&!p.sourceOffset&&!p.geometryOrigin;
}
export function ileocaecalReferenceLimitation(p:Part):string|undefined {
 if(p.id!=='HRA:VH_F_ileocecal_valve'||!reference(p))return;
 return 'Original HRA reference reviewed as an open annular sheet. It does not establish finite valve tissue, separate ileocolic/ileocaecal lips, frenula or wall layers. Adjoining terminal ileal source overlaps this sheet; registration to the native body remains unverified. Broader caecal placement is source correspondence only.';
}
export function ileocaecalReferenceDescription(p:Part):{text:string;url:string}|undefined {
 if(!ileocaecalReferenceLimitation(p))return;
 return {text:'The ileocaecal junction leads from the terminal ileum into the caecum. Projecting lips and their joining frenula surround the ileal orifice. Their muscular and mucosal tissues help limit backward flow into the ileum.',url:'https://pubmed.ncbi.nlm.nih.gov/7278461/'};
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Source-bound anatomical associations. Primary human morphology evidence:
 * https://pubmed.ncbi.nlm.nih.gov/7278461/ and /38978915/.
 * No contains/partOf/continuous relation is inferred from proximity. */
export function ileocaecalReferenceRelations(parts:Part[]):Map<string,ResolvedRelation[]> {
 const old=cache.get(parts);if(old)return old;
 const graph=new Map<string,ResolvedRelation[]>(),valve=parts.find(p=>p.id==='HRA:VH_F_ileocecal_valve'&&reference(p));
 if(valve)for(const [id,note,via] of [
  ['HRA:VH_F_caecum','Ileocaecal junction source association. The annular valve reference shares source controls with the caecal sheet; complete tissue layers, wall junction and native-body registration remain unverified.',undefined],
  ['HRA:VH_F_ileum_terminal','Terminal ileum–ileocaecal junction source association. The separate original sheets overlap; this is not a continuous, accepted tissue connection. Named lips, frenula and all-layer invagination remain unrepresented.',undefined],
  ['HRA:VH_F_vermiform_appendix','Caecal neighborhood association through the caecum. The appendiceal opening is separate from the ileal orifice; this does not make the appendix a valve component or establish tissue contact.',['Caecum']],
 ] as const){
  const target=parts.find(p=>p.id===id&&reference(p));if(!target)continue;
  for(const [from,to] of [[valve,target],[target,valve]]){
   const links=graph.get(from.id)??[];links.push({kind:'adjacent',target:to,note,via:via?[...via]:undefined,viaModeled:via?parts.some(p=>p.id==='HRA:VH_F_caecum'&&reference(p)):undefined});graph.set(from.id,links);
  }
 }
 cache.set(parts,graph);return graph;
}
