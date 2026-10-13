import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';

const names:Record<string,string>={
 'ZA:Free taenia':'Free taenia',
 'ZA:Mesocolic taenia':'Mesocolic taenia',
 'ZA:Omental taenia':'Omental taenia',
};
/** Named native reference sheets, not accepted finite muscle coats. */
export function nativeColonicTaenia(p:Part):boolean {
 return !p.suppressed&&p.system==='digestive'&&names[p.id]===p.name&&p.sourceId===p.name;
}
export function colonicTaeniaSourceLimitation(p:Part):string|undefined {
 if(nativeColonicTaenia(p)&&'colonicSurfaceReference' in p)return 'Conforming display footprint of a native open smooth-muscle band reference. It follows the existing exterior to avoid depth competition; this is not finite muscle tissue. Complete course, thickness and colonic attachment remain unverified.';
 return nativeColonicTaenia(p)?'Named native open smooth-muscle band reference. Finite muscle thickness, regional segmentation, complete course and source colon attachment remain unverified. The reference does not establish a complete longitudinal wall coat.':undefined;
}
/** Smooth muscle finish, isolated from the donor's shared ligament material.
 * No copied bitmap or additional mesh is required. */
export function colonicTaeniaFinish(p:Part) {
 return nativeColonicTaenia(p)?{color:[174,91,111] as [number,number,number],map:undefined,normalMap:undefined,roughness:.82,metalness:0,polygonOffset:true,polygonOffsetFactor:-1,polygonOffsetUnits:-1}:undefined;
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** OpenStax A&P 23.5: three bands of longitudinal smooth muscle.
 * Only same-source anatomical associations are shown. A complete source band
 * spans several regions; these links do not identify a regional tissue patch.
 */
export function colonicTaeniaRelations(parts:Part[]) {
 const old=cache.get(parts);if(old)return old;
 const graph=new Map<string,ResolvedRelation[]>();
 const pair=(a:Part,b:Part,note:string,via?:string[])=>{
  for(const [from,to] of [[a,b],[b,a]]){const rows=graph.get(from.id)??[];rows.push({kind:'adjacent',target:to,note,via,viaModeled:via?false:undefined});graph.set(from.id,rows);}
 };
 const native=(p:Part,name:string)=>!p.suppressed&&p.system==='digestive'&&p.id==='ZA:'+name&&p.name===name&&p.sourceId===name;
 for(const band of parts.filter(nativeColonicTaenia)) {
  for(const colon of parts)if(['Ascending colon','Transverse colon','Descending colon','Sigmoid colon'].some(name=>native(colon,name)))
   pair(band,colon,'Colon-level longitudinal muscle association. This whole named source band has not been partitioned into regional portions; its thickness, complete course and source tissue attachment remain unverified.',['Regional portion of the longitudinal muscle layer']);
  for(const support of parts){
   if(band.id==='ZA:Mesocolic taenia'&&native(support,'Mesocolon'))pair(band,support,'Mesocolic attachment landmark association. Source sheets, regional correspondence and a continuous tissue junction remain unverified.');
   if(band.id==='ZA:Omental taenia'&&native(support,'Greater omentum'))pair(band,support,'Omental attachment landmark association at the transverse colon. Source attachment coordinates and a continuous tissue junction remain unverified.',['Transverse colonic attachment']);
  }
 }
 cache.set(parts,graph);return graph;
}
