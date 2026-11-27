import type {Atlas} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {liverControlRegistry} from './liver-controls';
type Point=[number,number,number];
type Group={term:string;conceptId:string;partIds:string[];representation:string;completeOpeningRepresented:boolean;anchor:{partId:string;face:number;barycentric:Point;point:Point}};
type Marker={schema:string;hostPartId:string;bankSha256:string;termGroups:Group[]};

export function liverLandmarkRegistry(atlas:Atlas){
 const liver=liverControlRegistry(atlas),m=(atlas as Atlas&{liverLandmarkIntegration?:Marker}).liverLandmarkIntegration;
 const host=liver?.byId.get('TA98PART:ZA:Liver:capsule');
 if(!m||!host||m.schema!=='liver-hilar-landmark/v1'||m.hostPartId!==host.id||m.bankSha256!==liver!.marker.bankSha256||m.termGroups?.length!==1)return;
 const g=m.termGroups[0],a=g.anchor,c=atlas.concepts.find(c=>c.id===g.conceptId);
 if(g.term!=='A05.8.01.016'||g.representation!=='point-landmark'||g.completeOpeningRepresented!==false
  ||g.partIds.length!==1||g.partIds[0]!==host.id||!c||c.ta98Term!==g.term||c.elements.length!==1||c.elements[0]!==host.id
  ||a?.partId!==host.id||!Number.isSafeInteger(a.face)||a.face<0||a.face>=host.indexCount/3
  ||a.point?.length!==3||!a.point.every(Number.isFinite)||a.barycentric?.length!==3
  ||a.barycentric.some(w=>!Number.isFinite(w)||w<0||w>1)||Math.abs(a.barycentric.reduce((s,w)=>s+w,0)-1)>1e-12)return;
 return {marker:m,host};
}
export function liverLandmarkTermParts(atlas:Atlas){
 const r=liverLandmarkRegistry(atlas);return new Map(r?r.marker.termGroups.map(g=>[g.term,[r.host]]):[]);
}
export function liverLandmarkLabelReference(atlas:Atlas,choice:HierarchyChoice){
 const r=liverLandmarkRegistry(atlas);if(!r)return;
 const g=r.marker.termGroups[0];if((choice.id!==g.conceptId&&(choice.ta98Term??choice.terminology?.ta98)!==g.term)||choice.elements.length!==1||choice.elements[0]!==r.host.id)return;
 return {partId:g.anchor.partId,point:g.anchor.point};
}
