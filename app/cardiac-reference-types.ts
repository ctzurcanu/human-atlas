import type {Part} from './anatomy';
export type ReferencePoint=[number,number,number];
export type CardiacReferenceGeometry=
 | {type:'boundary-curve'|'opening-rim';paths:{points:ReferencePoint[];closed:boolean}[];filled?:false}
 | {type:'point-set';points:ReferencePoint[]}
 | {type:'surface-region';surfaces:{partId:string;faces:number[];areaMetres2:number;sourceBankSha256:string}[];state?:string}
 | {type:'finite-mesh';partIds:string[]}
 | {type:'finite-mesh'|'space-region';vertices:ReferencePoint[];faces:ReferencePoint[];filled?:false}
 | {type:'collection';terms:string[]};
const point=(p:unknown):p is ReferencePoint=>Array.isArray(p)&&p.length===3&&p.every(Number.isFinite);
export function validCardiacReferenceGeometry(g:CardiacReferenceGeometry,parts:Map<string,Part>,bank:string):boolean{
 if(!g||typeof g!=='object')return false;
 if(g.type==='collection')return Array.isArray(g.terms)&&g.terms.length>0&&g.terms.length<32&&new Set(g.terms).size===g.terms.length&&g.terms.every(t=>/^A12\.1\.\d{2}\.\d{3}$/.test(t));
 if(g.type==='point-set')return Array.isArray(g.points)&&g.points.length>0&&g.points.length<=32&&g.points.every(point);
 if(g.type==='boundary-curve'||g.type==='opening-rim')return (g.type!=='opening-rim'||g.filled===false)&&Array.isArray(g.paths)&&g.paths.length>0&&g.paths.length<=128&&g.paths.every(p=>typeof p.closed==='boolean'&&Array.isArray(p.points)&&p.points.length>=(p.closed?3:2)&&p.points.length<=16384&&p.points.every(point));
 if(g.type==='surface-region')return Array.isArray(g.surfaces)&&g.surfaces.length>0&&g.surfaces.length<=16&&g.surfaces.every(s=>{
  if(!s||typeof s!=='object')return false;
  const p=parts.get(s.partId);return p&&!p.suppressed&&s.sourceBankSha256===bank&&s.areaMetres2>0&&Number.isFinite(s.areaMetres2)&&Array.isArray(s.faces)&&s.faces.length>0&&new Set(s.faces).size===s.faces.length&&s.faces.every(f=>Number.isSafeInteger(f)&&f>=0&&f<p.indexCount/3);
 });
 if(g.type==='finite-mesh'&&'partIds'in g)return Array.isArray(g.partIds)&&g.partIds.length>0&&g.partIds.every(id=>parts.has(id)&&!parts.get(id)!.suppressed);
 if(g.type==='finite-mesh'||g.type==='space-region')return (g.type!=='space-region'||g.filled===false)&&Array.isArray(g.vertices)&&g.vertices.length>=4&&g.vertices.length<250000&&g.vertices.every(point)&&Array.isArray(g.faces)&&g.faces.length>=4&&g.faces.length<500000&&g.faces.every(f=>point(f)&&f.every(i=>Number.isSafeInteger(i)&&i>=0&&i<g.vertices.length)&&new Set(f).size===3);
 return false;
}
