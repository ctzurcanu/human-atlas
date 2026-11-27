import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';
import type {ResolvedRelation} from './anatomical-relations';

type Group={term:string;conceptId:string;partIds:string[]};
export type CordSegment={name:string;conceptId:string;partIds:string[];coveringPartIds:string[];parentTerm:string;yRange:[number,number]};
type CordCellRegion={partId:string;segment:string;nativeRegion:string};
export type SpinalBoundaryReference={kind:'lumen-boundary';ownerPartId:string;sourceBankSha256:string;sourceIndexStart:number;sourceIndexCount:number;physicalOwnerRetained:true};
type BoundaryPart=Part&{spinalBoundaryReference?:SpinalBoundaryReference};
type Marker={schema:'spinal-cord-native-core/v1'|'spinal-cord-native-core/v2'|'spinal-cord-native-core/v3'|'spinal-cord-native-core/v4';sex:'male'|'female';sourcePartId:string;wholeConceptId:string;partIds:string[];coveringPartIds?:string[];referencePartIds?:string[];bankSha256:string;sourcePlanSha256:string;termGroups:Group[];segments?:CordSegment[];cellRegions?:CordCellRegion[]};
const segmentNames=[...Array.from({length:8},(_,i)=>`C${i+1}`),...Array.from({length:12},(_,i)=>`T${i+1}`),...Array.from({length:5},(_,i)=>`L${i+1}`),...Array.from({length:5},(_,i)=>`S${i+1}`),'Co1'];
const cache=new WeakMap<Atlas,{marker:Marker;parts:Part[];referenceParts:Part[]}|null>();
const same=(a:readonly string[],b:readonly string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));

/** Only an explicit, reciprocal native-frame tissue assembly gets controls. */
export function spinalCordRegistry(atlas:Atlas){
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const m=(atlas as Atlas&{finiteSpinalIntegration?:Marker}).finiteSpinalIntegration;
 if(!m||!['spinal-cord-native-core/v1','spinal-cord-native-core/v2','spinal-cord-native-core/v3','spinal-cord-native-core/v4'].includes(m.schema)||m.sex!==atlas.sex
  ||(m.schema==='spinal-cord-native-core/v4'?!(m.partIds?.length>=372&&m.partIds.length<=4096):m.partIds?.length!==({'spinal-cord-native-core/v1':12,'spinal-cord-native-core/v2':96,'spinal-cord-native-core/v3':372}[m.schema]))
  ||!['bankSha256','sourcePlanSha256'].every(key=>/^[a-f0-9]{64}$/.test(m[key as 'bankSha256']))||!Array.isArray(m.termGroups))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const ids=[...m.partIds,...(m.coveringPartIds??[])];
 if(new Set(ids).size!==ids.length)return reject();
 const parts=ids.map(id=>byId.get(id)),first=parts[0];
 if(!first||parts.some(p=>!p||p.suppressed||p.system!=='nervous'||!p.pairedTissue||!p.sectionSolid||p.normalSigns===undefined
  ||p.tissueWallAssembly!==`spinal-cord:${m.sex}`||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals
  ||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return reject();
 const actual=parts as Part[];
 const actualById=new Map(actual.map(p=>[p.id,p]));
 const referenceIds=m.referencePartIds??[],referenceParts:Part[]=[];
 if(referenceIds.length>512||new Set([...ids,...referenceIds]).size!==ids.length+referenceIds.length)return reject();
 for(const id of referenceIds){
  const p=byId.get(id)as BoundaryPart|undefined,g=p?.spinalBoundaryReference,owner=g&&actualById.get(g.ownerPartId);
  if(!p||p.suppressed||p.pairedTissue||p.sectionSolid!==false||!g||g.kind!=='lumen-boundary'||g.physicalOwnerRetained!==true||g.sourceBankSha256!==m.bankSha256||!owner
   ||!Number.isSafeInteger(g.sourceIndexStart)||g.sourceIndexStart<0||g.sourceIndexStart%3||g.sourceIndexCount!==p.indexCount||p.indexCount<=0||p.indexCount%3||g.sourceIndexStart+p.indexCount>owner.indexCount
   ||p.indices!==owner.indices+g.sourceIndexStart*4||['chunk','positions','normals','normalSigns','vertexCount'].some(k=>p[k as keyof Part]!==owner[k as keyof Part])
   ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(owner.geometryOrigin))return reject();
  referenceParts.push(p);
 }
 const selectable=[...ids,...referenceIds];
 for(const p of actual){
  let end=0;
  for(const g of p.tissueDrawGroups??[]){
   if(g.start!==end||g.count<=0||g.count%3)return reject();end+=g.count;
   if(g.neighbor&&actualById.get(g.neighbor)?.tissueDrawGroups?.filter(r=>r.neighbor===p.id).reduce((n,r)=>n+r.count,0)!==g.count)return reject();
  }
  if(end!==p.indexCount)return reject();
 }
 for(const g of m.termGroups){
  const c=concepts.get(g.conceptId);
  if(!/^A14\.1\.0[12]\.\d{3}$/.test(g.term)||!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'
   ||!same(c.elements,g.partIds)||g.partIds.some(id=>!selectable.includes(id)))return reject();
 }
 const whole=concepts.get(m.wholeConceptId);
 if(!whole||!same(whole.elements,m.partIds))return reject();
 if(m.schema==='spinal-cord-native-core/v3'||m.schema==='spinal-cord-native-core/v4'){
  if(m.segments?.length!==31||m.coveringPartIds?.length!==31||!same(m.segments.map(s=>s.name),segmentNames))return reject();
  for(const [i,s]of m.segments.entries()){
   const c=concepts.get(s.conceptId),parent=m.termGroups.find(g=>g.term===s.parentTerm);
   if(!c||!parent||c.partOf!==parent.conceptId||c.ta98Term||!same(c.elements,s.partIds)||(m.schema.endsWith('/v4')?s.partIds.length<12:s.partIds.length!==12)||s.coveringPartIds?.length!==1
    ||s.partIds.some(id=>!m.partIds.includes(id)||!id.endsWith('@'+s.name))
    ||s.coveringPartIds.some(id=>!m.coveringPartIds!.includes(id)||(m.schema.endsWith('/v4')?!(id.startsWith(`TA98CORD:${m.sex}:`)&&id.endsWith(`:pia@${s.name}`)):id!==`TA98CORD:${m.sex}:pia@${s.name}`))
    ||s.yRange?.length!==2||!s.yRange.every(Number.isFinite)||s.yRange[0]>=s.yRange[1]
    ||s.name!==segmentNames[i]||(i>0&&Math.abs(m.segments[i-1].yRange[0]-s.yRange[1])>1e-12))return reject();
  }
  if(!same(m.segments.flatMap(s=>s.partIds),m.partIds)||!same(m.segments.flatMap(s=>s.coveringPartIds),m.coveringPartIds))return reject();
 }
 if(m.schema==='spinal-cord-native-core/v4'){
  if(!m.cellRegions||!same(m.cellRegions.map(r=>r.partId),ids))return reject();
  const regions=new Map(m.cellRegions.map(r=>[r.partId,r]));
  for(const s of m.segments!){
   const native=new Set<string>();
   for(const id of [...s.partIds,...s.coveringPartIds]){
    const row=regions.get(id),part=actualById.get(id)!;
    if(!row||row.segment!==s.name||part.tissueWallLayer!==row.nativeRegion.split('-')[0])return reject();
    native.add(row.nativeRegion);
   }
   const expected=['pia',...['left','right'].flatMap(side=>['grey-anterior','grey-posterior','grey-intermediate','white-anterior','white-lateral','white-posterior'].map(region=>`${region}-${side}`))];
   if(!same([...native],expected))return reject();
  }
 }
 const result={marker:m,parts:actual,referenceParts};cache.set(atlas,result);return result;
}

export function spinalCordTermParts(atlas:Atlas):Map<string,Part[]>{
 const r=spinalCordRegistry(atlas),byId=new Map(r?[...r.parts,...r.referenceParts].map(p=>[p.id,p]):[]);
 return new Map(r?.marker.termGroups.map(g=>[g.term,g.partIds.map(id=>byId.get(id)!)])??[]);
}

export function spinalCordSegments(atlas:Atlas){return spinalCordRegistry(atlas)?.marker.segments??[];}

export function spinalCordPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const r=spinalCordRegistry(atlas);
 if(!r||choice.id==='selection-set')return null;
 const physical=choice.elements.filter(id=>!r.marker.referencePartIds?.includes(id));
 const covered=!!r.marker.coveringPartIds?.length&&same(physical,[...r.marker.partIds,...r.marker.coveringPartIds]);
 if(!covered&&!same(physical,r.marker.partIds))return null;
 return [...(covered?['pia']:[]),'white','grey'].map(key=>({name:key==='white'?'White matter':key==='grey'?'Grey matter':'Pia mater',ids:r.parts.filter(p=>p.tissueWallLayer===key).map(p=>p.id)}));
}

const coveringCache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Explicit stored pial ownership, without proximity-based root claims. */
export function spinalCordCoveringRelations(parts:Part[]){
 let result=coveringCache.get(parts);if(result)return result;
 result=new Map();coveringCache.set(parts,result);
 for(const sex of ['female','male']){
  const pial=parts.filter(p=>p.id.startsWith(`TA98CORD:${sex}:`)&&/(?:^|:)pia(?:@[^:]+)?$/.test(p.id)&&!p.suppressed&&p.pairedTissue&&p.sectionSolid&&p.tissueWallLayer==='pia'&&p.tissueWallAssembly===`spinal-cord:${sex}`);
  for(const pia of pial){
   const segment=pia.id.split('@')[1];
   const cord=parts.filter(p=>p.id.startsWith(`TA98CORD:${sex}:`)&&!p.suppressed&&p.pairedTissue&&p.sectionSolid&&['white','grey'].includes(p.tissueWallLayer??'')
    &&(!segment||p.id.endsWith('@'+segment))&&p.tissueWallAssembly===pia.tissueWallAssembly&&p.chunk===pia.chunk&&p.positions===pia.positions&&p.normals===pia.normals&&JSON.stringify(p.geometryOrigin)===JSON.stringify(pia.geometryOrigin));
   for(const part of cord){result.set(part.id,[{kind:'coveredBy',target:pia}]);}
   result.set(pia.id,cord.map(target=>({kind:'covers' as const,target})));
  }
 }
 return result;
}
