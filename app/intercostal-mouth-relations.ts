import type {Atlas,Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import signatures from './data/intercostal-mouth-signatures.json';

type ConnectedPart=Part&{verifiedIntercostalMouths?:string[]};
type Expected={id:string;bank:string;positions:number;indices:number;vertexCount:number;indexCount:number;origin:number[]};
const pairs=signatures as [Expected,Expected][];

/** Expose only the source-mouth correspondences verified on these exact banks. */
export function withVerifiedIntercostalMouths(atlas:Atlas):Atlas{
 const data=atlas as Atlas&{chunks:{sha256?:string}[]};
 const byId=new Map(atlas.parts.map(part=>[part.id,part]));
 const valid=(e:Expected)=>{const p=byId.get(e.id);return !!p&&!p.suppressed&&p.system==='arterial'&&p.positions===e.positions&&p.indices===e.indices&&p.vertexCount===e.vertexCount&&p.indexCount===e.indexCount&&p.geometryOrigin?.length===3&&p.geometryOrigin.every((v,i)=>v===e.origin[i])&&!p.geometryScale&&!p.geometryInstance&&!p.sourceOffset&&data.chunks[p.chunk]?.sha256===e.bank;};
 const links=new Map<string,string[]>();
 for(const [a,b] of pairs)if(valid(a)&&valid(b))for(const [from,to] of [[a.id,b.id],[b.id,a.id]])links.set(from,[...(links.get(from)??[]),to]);
 if(!links.size)return atlas;
 return {...atlas,parts:atlas.parts.map(part=>links.has(part.id)?{...part,verifiedIntercostalMouths:links.get(part.id)}:part)};
}

/** Matched openings are geometric correspondences, not a whole-lumen patency claim. */
export function intercostalMouthRelations(source:Part,parts:Part[]):ResolvedRelation[]{
 const sourcePart=source as ConnectedPart;
 if(source.suppressed||source.system!=='arterial')return [];
 return (sourcePart.verifiedIntercostalMouths??[]).flatMap(id=>{
  const target=parts.find(part=>part.id===id) as ConnectedPart|undefined;
  if(!target||target.suppressed||target.system!=='arterial'||!target.verifiedIntercostalMouths?.includes(source.id))return [];
  return [{kind:'continuous' as const,target,note:'These source vessel openings correspond within 0.001 mm; complete lumen continuity has not been verified.'}];
 });
}
