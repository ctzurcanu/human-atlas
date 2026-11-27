import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {LiverSurfaceFeature} from './liver-surface-features';
import {cardiacIntegration} from './heart-wall';

type Marker={schema:'cardiac-surface-features/v1';hostPartId:string;bankSha256:string;termGroups:LiverSurfaceFeature[]};
const cache=new WeakMap<Atlas,Marker|null>();
const allowed=new Set(['A12.1.00.002','A12.1.00.003','A12.1.00.004','A12.1.00.005','A12.1.00.007','A12.1.00.031']);
export function cardiacFeatureRegistry(atlas:Atlas):Marker|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;},heart=cardiacIntegration(atlas);
 const m=(atlas as Atlas&{cardiacFeatureIntegration?:Marker}).cardiacFeatureIntegration;
 const host=m&&atlas.parts.find(p=>p.id===m.hostPartId),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 if(!m||!heart||!host||host.suppressed||m.schema!=='cardiac-surface-features/v1'||!heart.wallPartIds.includes(host.id)||host.tissueWallLayer!=='epicardium'||m.bankSha256!==heart.coupledBankSha256||!/^[a-f0-9]{64}$/.test(m.bankSha256)||m.termGroups.length!==allowed.size||new Set(m.termGroups.map(g=>g.term)).size!==allowed.size)return reject();
 for(const g of m.termGroups){const c=concepts.get(g.conceptId),s=g.surface,a=g.anchor;
  if(!allowed.has(g.term)||g.representation!=='surface-region'||!g.boundaryProvisional||!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||c.elements.length!==1||c.elements[0]!==host.id||g.partIds.length!==1||g.partIds[0]!==host.id||s.partId!==host.id||s.sourceBankSha256!==m.bankSha256||!s.faces.length||new Set(s.faces).size!==s.faces.length||s.faces.some(f=>!Number.isSafeInteger(f)||f<0||f>=host.indexCount/3)||a.partId!==host.id||!s.faces.includes(a.face)||a.barycentric.length!==3||a.barycentric.some(w=>!Number.isFinite(w)||w<0||w>1)||Math.abs(a.barycentric.reduce((sum,w)=>sum+w,0)-1)>1e-12||a.point.length!==3||!a.point.every(Number.isFinite))return reject();
  for(const f of s.faces){const offset=f*3;if(!host.tissueDrawGroups?.some(d=>!d.neighbor&&offset>=d.start&&offset<d.start+d.count))return reject();}
 }
 cache.set(atlas,m);return m;
}
export function selectedCardiacSurfaceFeature(atlas:Atlas,choice:HierarchyChoice|undefined){
 const m=cardiacFeatureRegistry(atlas);if(!m||!choice||choice.elements.length!==1||choice.elements[0]!==m.hostPartId)return;
 return m.termGroups.find(g=>g.conceptId===choice.id||g.term===(choice.ta98Term??choice.terminology?.ta98));
}
export function cardiacSurfaceLabelReference(atlas:Atlas,choice:HierarchyChoice){
 const g=selectedCardiacSurfaceFeature(atlas,choice);return g?{partId:g.anchor.partId,point:g.anchor.point}:undefined;
}
