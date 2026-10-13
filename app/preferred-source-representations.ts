import type {Atlas,Part} from './anatomy';
import preferences from './data/preferred-source-representations.json';

type SourceDecision={id:string;name:string;conceptId:string;basis:string;term:string;parts:string[]};

/** Retire only alternatives for which this catalogue already records a complete
 * source selection. Unresolved alternatives, bilateral parts and wall cells stay
 * available. Geometry buffers and the original part order remain unchanged. */
export function withPreferredSourceRepresentations(atlas:Atlas):Atlas{
 const source=atlas as unknown as Record<string,unknown>;
 if(source.reviewOnly!==true||!atlas.sex)return atlas;
 const decisions=preferences[atlas.sex] as SourceDecision[],byId=new Map(atlas.parts.map(p=>[p.id,p])),hidden=new Set(atlas.defaultHidden??[]);
 const retired=new Map<string,string[]>();
 for(const d of decisions){
  const old=byId.get(d.id),selection=source[d.basis] as {rules?:Record<string,string[]>}|undefined;
  if(!old||old.suppressed||old.pairedTissue||old.name!==d.name||old.conceptId!==d.conceptId||!hidden.has(d.id))continue;
  const rule=selection?.rules?.[d.term];
  if(!rule||!d.parts.length||d.parts.some(id=>!rule.includes(id)||!byId.has(id)||byId.get(id)!.suppressed||hidden.has(id)))continue;
  retired.set(d.id,d.parts);
 }
 if(!retired.size)return atlas;
 const parts=atlas.parts.map(p=>retired.has(p.id)?{...p,suppressed:true,replacementParts:retired.get(p.id)!}:p);
 const concepts=atlas.concepts.map(c=>{
  if(!c.elements.some(id=>retired.has(id)))return c;
  return {...c,elements:[...new Set(c.elements.flatMap(id=>retired.get(id)??[id]))]};
 });
 // Hidden old references must never hide their active replacement after URL
 // restoration. Explicit legacy selections still resolve to the preferred set.
 return {...atlas,parts:parts as Part[],concepts,defaultHidden:atlas.defaultHidden?.filter(id=>!retired.has(id))};
}

/** A combined TA98 catalogue keeps unfinished donor meshes for local repair.
 * They are not members of the displayed body. A visibility reset or parent
 * selection must not turn them into a second skeleton or organ set. */
export function withHiddenReviewAlternatives(atlas:Atlas,inspectAlternatives=false):Atlas{
 if(inspectAlternatives)return atlas;
 const hidden=new Set(atlas.defaultHidden??[]);
 if(!atlas.parts.some(p=>hidden.has(p.id)&&!p.suppressed))return atlas;
 const parts=atlas.parts.map(p=>hidden.has(p.id)&&!p.suppressed?{...p,suppressed:true}:p);
 const unavailable=new Set(parts.filter(p=>p.suppressed&&!p.replacementParts?.length).map(p=>p.id));
 const concepts=atlas.concepts.flatMap(c=>{
  if(!c.elements.some(id=>unavailable.has(id)))return [c];
  const elements=c.elements.filter(id=>!unavailable.has(id));
  return elements.length?[{...c,elements}]:[];
 });
 return {...atlas,parts,concepts};
}
