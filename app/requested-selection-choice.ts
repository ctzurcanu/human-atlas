import type {Atlas,Concept} from './anatomy';

/** Preserve an explicitly requested concept only for its complete mesh set.
 * A partial region, mixed selection or ambiguous name must keep its own identity.
 */
export function requestedSelectionChoice(atlas:Atlas,terms:string[],selected:string[]):Concept|null{
 const requested=[...new Set(terms)];
 if(requested.length!==1||!selected.length)return null;
 const term=requested[0];
 const exact=atlas.concepts.filter(concept=>concept.id===term);
 const matches=exact.length?exact:atlas.concepts.filter(concept=>concept.name.toLowerCase()===term.toLowerCase());
 if(matches.length!==1)return null;
 const concept=matches[0],pieces=new Set(selected),elements=new Set(concept.elements);
 return elements.size===pieces.size&&[...elements].every(id=>pieces.has(id))?concept:null;
}
