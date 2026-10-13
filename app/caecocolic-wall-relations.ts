import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
const identity=/^LOCAL:(female|male):CAECOWALL:(caecum|appendix|ascending|transverse|descending|sigmoid|rectum):(serosa|subserosa|longitudinal-muscle|circular-muscle|submucosa|muscularis-mucosae|lamina-propria|epithelium)(?::(interband|free-taenia|mesocolic-taenia|omental-taenia))?$/;
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
function qualified(p:Part){
 const match=p.id.match(identity);
 return match&&!p.suppressed&&p.system==='digestive'&&p.pairedTissue&&p.sectionSolid&&!p.sourceOffset&&
 (!p.caecocolicLayer||p.caecocolicLayer===match[3])&&(!match[4]||match[3]==='longitudinal-muscle'&&p.caecocolicBand===match[4])&&p.caecocolicWallAssembly===`caecocolic:${match[1]}`&&p.caecocolicRegion===match[2]&&p.bowelWallSource&&
 p.geometryOrigin?.length===3&&p.geometryOrigin.every(Number.isFinite)&&p.tissueDrawGroups?.length?match:null;
}
export function caecocolicWallInterface(source:Part,target:Part):boolean{
 const a=qualified(source),b=qualified(target);
 return !!a&&!!b&&a[1]===b[1]&&(a[2]!==b[2]||a[4]!==b[4])&&a[3]===b[3]&&source.chunk===target.chunk&&source.positions===target.positions&&source.normals===target.normals&&source.vertexCount===target.vertexCount&&JSON.stringify(source.geometryOrigin)===JSON.stringify(target.geometryOrigin)&&source.tissueDrawGroups!.some(g=>g.neighbor===target.id&&g.count>0&&target.tissueDrawGroups!.some(h=>h.neighbor===source.id&&h.count===g.count));
}
/** Reciprocal, stored same-coat interfaces, not inferred organ proximity.
 * Clinical tissue allocation and external native attachments remain unverified. */
export function caecocolicWallRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const previous=cache.get(parts);if(previous)return previous;
 const graph=new Map<string,ResolvedRelation[]>(),byId=new Map(parts.map(p=>[p.id,p]));
 for(const source of parts){
  const a=qualified(source);if(!a)continue;
  const linked=new Set<string>();
  for(const group of source.tissueDrawGroups!){
   if(!group.neighbor||group.count<=0||linked.has(group.neighbor))continue;
   const target=byId.get(group.neighbor);
   if(!target||!caecocolicWallInterface(source,target))continue;
   linked.add(target.id);
   const relations=graph.get(source.id)??[];
   relations.push({kind:'continuous',target,note:'Shared regional interface in the stored large-intestinal tissue bank. The matching coat uses common surface controls and normals; clinical regional boundaries, peritoneal distribution and native adjoining-organ attachments remain unverified.'});graph.set(source.id,relations);
  }
 }
 cache.set(parts,graph);return graph;
}
