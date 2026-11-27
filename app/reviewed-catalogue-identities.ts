import type {Atlas} from './anatomy';
import metadata from './data/ta98-metadata.json';
import {withLaryngealSystems} from './laryngeal-relations';
import {cardiacIntegration,cardiacTermParts} from './heart-wall';
import {bloodReferenceConcepts} from './blood-reference-identities';
import {withImportedHandVeinIdentities} from './imported-hand-identities';

// Only reviews explicitly bound to source concept IDs can supersede a stale
// embedded catalogue tag. Generic name rules cannot overwrite authored tags.
const scoped=metadata.reviewedCatalogueIdentities as Record<string,{term:string;kind:'exact'|'parent'}>;
export function withReviewedCatalogueIdentities(atlas:Atlas):Atlas{
 atlas=withImportedHandVeinIdentities(atlas);
 let changed=false;
 const installedCardiacTerms=cardiacTermParts(atlas);
 const installedCardiacConcepts=new Set(cardiacIntegration(atlas)?.termGroups?.filter(g=>g.kind==='exact'&&installedCardiacTerms.has(g.term)).map(g=>g.conceptId));
 const concepts=atlas.concepts.map(concept=>{
  // A later installed, explicitly owned cardiac partition supersedes the
  // older donor-wide identity correction, without changing other sources.
  if(installedCardiacConcepts.has(concept.id))return concept;
  const rule=scoped[concept.id],match=(metadata.byConceptMatch as Record<string,{term:string;kind:string;method:string}>)[concept.id];
  if(!rule||match?.method!=='reviewed'||match.term!==rule.term||match.kind!==rule.kind)return concept;
  if(concept.ta98Term===rule.term&&concept.ta98Kind===rule.kind)return concept;
  changed=true;return {...concept,ta98Term:rule.term,ta98Kind:rule.kind as 'exact'|'parent'};
 });
 const bloodConcepts=bloodReferenceConcepts(atlas,concepts);
 return withLaryngealSystems(changed||bloodConcepts!==concepts?{...atlas,concepts:bloodConcepts}:atlas);
}
