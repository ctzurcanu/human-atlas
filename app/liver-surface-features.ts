import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {liverControlRegistry} from './liver-controls';

type Point=[number,number,number];
export type LiverSurfaceFeature={term:string;conceptId:string;partIds:string[];representation:'surface-region';surface:{partId:string;faces:number[];sourceBankSha256:string};anchor:{partId:string;face:number;barycentric:Point;point:Point};boundaryProvisional:boolean};
type Marker={schema:string;hostPartId:string;bankSha256:string;sourcePlanSha256:string;termGroups:LiverSurfaceFeature[];provisional:boolean};
const cache=new WeakMap<Atlas,Marker|null>();
const terms=[2,3,4,5,6,7,8,9,10,12,13,14,18,19,20,21,22,23].map(n=>`A05.8.01.${String(n).padStart(3,'0')}`);
const sha=(s:unknown)=>typeof s==='string'&&/^[a-f0-9]{64}$/.test(s);

/** Surface annotations bind actual carrier facets. They add neither tissue
 * solids nor microscopic representations, and retain provisional boundaries. */
export function liverFeatureRegistry(atlas:Atlas):Marker|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;},liver=liverControlRegistry(atlas);
 const m=(atlas as Atlas&{liverFeatureIntegration?:Marker}).liverFeatureIntegration;
 const host=liver?.byId.get('TA98PART:ZA:Liver:capsule');
 if(!m||!liver||!host||m.schema!=='liver-surface-features/v1'||m.hostPartId!==host.id
  ||m.bankSha256!==liver.marker.bankSha256||!sha(m.sourcePlanSha256)||m.provisional!==true
  ||!Array.isArray(m.termGroups)||![terms.length,terms.length+1].includes(m.termGroups.length)||new Set(m.termGroups.map(g=>g.term)).size!==m.termGroups.length)return reject();
 const concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const outerCount=(host as Part&{liverReferenceTriangleCount?:number}).liverReferenceTriangleCount??((host.tissueDrawGroups?.find(g=>g.start===0&&!g.neighbor)?.count??0)/3);
 if(!Number.isSafeInteger(outerCount)||outerCount<1||outerCount>host.indexCount/3)return reject();
 for(const g of m.termGroups){const c=concepts.get(g.conceptId),surface=g.surface,a=g.anchor;
  if(![...terms,'A05.8.01.017'].includes(g.term)||g.representation!=='surface-region'||g.boundaryProvisional!==true
   ||!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||c.elements.length!==1||c.elements[0]!==host.id
   ||g.partIds.length!==1||g.partIds[0]!==host.id||surface?.partId!==host.id||surface.sourceBankSha256!==m.bankSha256
   ||!Array.isArray(surface.faces)||!surface.faces.length||new Set(surface.faces).size!==surface.faces.length
   ||surface.faces.some(f=>!Number.isSafeInteger(f)||f<0||f>=outerCount)
   ||a?.partId!==host.id||!surface.faces.includes(a.face)||!Array.isArray(a.point)||a.point.length!==3||!a.point.every(Number.isFinite)
   ||!Array.isArray(a.barycentric)||a.barycentric.length!==3||a.barycentric.some(w=>!Number.isFinite(w)||w<0||w>1)||Math.abs(a.barycentric.reduce((s,w)=>s+w,0)-1)>1e-12)return reject();
 }
 cache.set(atlas,m);return m;
}

export function liverSurfaceTermParts(atlas:Atlas):Map<string,Part[]>{
 const m=liverFeatureRegistry(atlas),host=m&&atlas.parts.find(p=>p.id===m.hostPartId);
 return new Map(m&&host?m.termGroups.map(g=>[g.term,[host]]):[]);
}
export function selectedLiverSurfaceFeature(atlas:Atlas,choice:HierarchyChoice|undefined){
 const m=liverFeatureRegistry(atlas);if(!m||!choice||choice.elements.length!==1||choice.elements[0]!==m.hostPartId)return;
 return m.termGroups.find(g=>g.conceptId===choice.id||g.term===(choice.ta98Term??choice.terminology?.ta98));
}
export function liverSurfaceLabelReference(atlas:Atlas,choice:HierarchyChoice){
 const g=selectedLiverSurfaceFeature(atlas,choice);return g?{partId:g.anchor.partId,point:g.anchor.point}:undefined;
}
