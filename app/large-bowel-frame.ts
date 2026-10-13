import type {Part} from './anatomy';
import type {RelationKind} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const names=new Set(['caecum','cecum','large intestine cecum','vermiform appendix','appendix','ascending colon','transverse colon','descending colon','large intestine descending colon','sigmoid colon','rectum','large intestine rectum']);
const nativeLargeBowel=(p:Part)=>!!nativePartFrame(p)&&p.system==='digestive'&&(!!p.caecocolicRegion||names.has(p.name.toLowerCase()));
const kinds=new Set<RelationKind>(['before','after','arterial','venous','innervation','supplies','drains','innervates','lymphaticDrainage','lymphaticTributaries']);
const replacementCache=new WeakMap<Part[],{sources:Map<string,Part>;outer:Map<string,Part>}>();
function replacements(parts:Part[]){
 const prior=replacementCache.get(parts);if(prior)return prior;
 const byId=new Map(parts.map(p=>[p.id,p])),sources=new Map<string,Part>(),outer=new Map<string,Part>();
 for(const source of parts){
  if(!source.suppressed||source.system!=='digestive'||(!source.caecocolicRegion&&!names.has(source.name.toLowerCase()))||!source.replacementParts?.length)continue;
  // Explicit retired references can name their registered local wall. This
  // does not assign the old FJ/ZA geometry the replacement's body frame.
  const parentFrame=nativePartFrame(source),firstReplacement=byId.get(source.replacementParts[0]);
  const replacementFrame=firstReplacement&&nativePartFrame(firstReplacement);
  const registeredMaleWall=replacementFrame==='LOCAL:male:'&&(parentFrame==='ZA:'||/^FJ\d+$/.test(source.id));
  const frame=registeredMaleWall?replacementFrame:parentFrame;
  if(!frame)continue;
  let valid=true;
  const expand=(id:string,path:Set<string>):Part[]=>{
   const p=byId.get(id);
   if(path.has(id)||!p||p.system!=='digestive'||!p.pairedTissue||p.sourceOffset||p.bowelWallSource!==source.id||nativePartFrame(p)!==frame||!p.caecocolicRegion||!p.caecocolicWallAssembly||
    registeredMaleWall&&(!/^LOCAL:male:CAECOWALL:/.test(p.id)||p.caecocolicWallAssembly!=='caecocolic:male')){valid=false;return [];}
   if(!p.suppressed)return [p];
   if(!p.replacementParts?.length){valid=false;return [];}
   const next=new Set(path);next.add(id);return p.replacementParts.flatMap(child=>expand(child,next));
  };
  const active=source.replacementParts.flatMap(id=>expand(id,new Set([source.id])));if(!valid||!active.length)continue;
  const first=active[0];
  if(active.some(p=>p.caecocolicRegion!==first.caecocolicRegion||p.caecocolicWallAssembly!==first.caecocolicWallAssembly||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))continue;
  const covering=active.filter(p=>p.caecocolicLayer==='serosa');if(covering.length!==1)continue;
  outer.set(source.id,covering[0]);for(const member of active)sources.set(member.id,source);
 }
 const result={sources,outer};replacementCache.set(parts,result);return result;
}
/** Retired source surfaces resolve to their registered active covering for an
 * organ-level route. This is not a claim of local coat-to-coat continuity.
 */
export function largeBowelRouteTargets(p:Part,parts:Part[]):Part[]{
 const target=replacements(parts).outer.get(p.id);return target?[target]:[p];
}
export function largeBowelWallSource(p:Part,parts:Part[]):Part|undefined{return replacements(parts).sources.get(p.id);}

/** Functional links belong to the complete registered organ assembly, not an
 * arbitrary coat, partial selection, or independently placed donor. */
export function largeBowelWholeSource(selected:Part[],parts:Part[]):Part|undefined{
 if(!selected.length||selected.some(p=>p.suppressed))return;
 const registry=replacements(parts),source=registry.sources.get(selected[0].id);
 if(!source||selected.some(p=>registry.sources.get(p.id)!==source))return;
 const ids=new Set(selected.map(p=>p.id));
 const members=[...registry.sources].filter(([,parent])=>parent===source).map(([id])=>id);
 return ids.size===members.length&&members.every(id=>ids.has(id))?source:undefined;
}

/** A matching organ name cannot join a native bowel to an independently
 * registered donor, or supply it from the other sex/source frame. This guard
 * does not certify registration of legacy unframed catalogue pairs.
 */
export function largeBowelFrameCompatible(kind:RelationKind,a:Part,b:Part):boolean{
 if(!kinds.has(kind)||!nativeLargeBowel(a)&&!nativeLargeBowel(b))return true;
 const frame=nativePartFrame(a);
 return !a.suppressed&&!b.suppressed&&!!frame&&frame===nativePartFrame(b);
}
