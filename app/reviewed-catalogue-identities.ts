import type {Atlas} from './anatomy';
import metadata from './data/ta98-metadata.json';

// Only reviews explicitly bound to source concept IDs can supersede a stale
// embedded catalogue tag. Generic name rules cannot overwrite authored tags.
const scoped=metadata.reviewedCatalogueIdentities as Record<string,{term:string;kind:'exact'|'parent'}>;
export function withReviewedCatalogueIdentities(atlas:Atlas):Atlas{
 let changed=false;
 const concepts=atlas.concepts.map(concept=>{
  const rule=scoped[concept.id],match=(metadata.byConceptMatch as Record<string,{term:string;kind:string;method:string}>)[concept.id];
  if(!rule||match?.method!=='reviewed'||match.term!==rule.term||match.kind!==rule.kind)return concept;
  if(concept.ta98Term===rule.term&&concept.ta98Kind===rule.kind)return concept;
  changed=true;return {...concept,ta98Term:rule.term,ta98Kind:rule.kind as 'exact'|'parent'};
 });
 return changed?{...atlas,concepts}:atlas;
}
