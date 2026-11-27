import type {Atlas,Concept,Part} from './anatomy';
/** A surface (a part with no volume) is shown selected together with the smallest whole unit it belongs to,
 * drawn transparent around it. Returns the ids of those whole units' parts (never the selected surfaces themselves). */
export const isSurfaceOnly=(part:Part|undefined)=>!!part&&(part.surfaceOnly===true||part.representationClass==='surface');
/** The smallest whole unit of each selected surface, with the parts it is drawn from. */
export function surfaceWholes(atlas:Pick<Atlas,'parts'|'concepts'>,selected:readonly string[]):{concept:Concept;ids:string[]}[]{
 const parts=new Map(atlas.parts.map(p=>[p.id,p])),chosen=new Set(selected),found=new Map<string,{concept:Concept;ids:string[]}>();
 for(const surface of selected.map(id=>parts.get(id)).filter(isSurfaceOnly) as Part[]){
  let best:{concept:Concept;ids:string[]}|undefined;
  for(const c of atlas.concepts){
   if(!c.elements.includes(surface.id))continue;
   const ids=c.elements.filter(id=>{const p=parts.get(id);return p&&!p.suppressed&&!isSurfaceOnly(p)&&!/^DMT[AM]:/.test(id)&&!chosen.has(id)});
   if(ids.length&&(!best||ids.length<best.ids.length))best={concept:c,ids};
  }
  if(best)found.set(best.concept.id,best);
 }
 return [...found.values()];
}
/** The whole's own subdivisions (other than the selected surfaces) stay out of the picture: the transparent whole is drawn from the bone itself. */
export const surfaceHidden=(atlas:Pick<Atlas,'parts'|'concepts'>,selected:readonly string[]):Set<string>=>{
 const chosen=new Set(selected),parts=new Map(atlas.parts.map(p=>[p.id,p]));
 return new Set(surfaceWholes(atlas,selected).flatMap(w=>w.concept.elements.filter(id=>!chosen.has(id)&&parts.get(id)&&(parts.get(id)!.ta98Term!==undefined&&/^DMT[AM]:/.test(id)))));
};
export const surfaceCompanions=(atlas:Pick<Atlas,'parts'|'concepts'>,selected:readonly string[]):Set<string>=>new Set(surfaceWholes(atlas,selected).flatMap(w=>w.ids));
