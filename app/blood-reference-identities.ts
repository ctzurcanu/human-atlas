import type {Atlas,Concept} from './anatomy';
import metadata from './data/ta98-metadata.json';

const matches=metadata.byConceptMatch as Record<string,{term:string;kind:string}>;
const paths=metadata.byConceptPath as Record<string,string[]>;

/** Remove contradicted legacy blood-vessel aggregates without reassigning tissue. */
export function bloodReferenceConcepts(atlas:Atlas,concepts:Concept[]):Concept[]{
 const parts=new Map(atlas.parts.map(part=>[part.id,part]));
 let changed=false;
 const result=concepts.map(concept=>{
  // Authored identities and bounded feature carriers have their own contracts.
  // This guard applies only to older, automatically matched vessel concepts.
  const inferred=matches[concept.id];
  if(concept.ta98Term||inferred?.kind!=='exact'||!/^A12\.[23]\./.test(inferred.term))return concept;
  const elements=concept.elements.filter(id=>{
   const part=parts.get(id),term=part?.ta98Term;
   if(!part||!term||term===inferred.term)return true;
   const identity=matches[part.conceptId],path=paths[part.conceptId];
   // Missing or disagreeing identity evidence is not grounds for deleting a
   // reference. A verified descendant remains a legitimate whole assembly.
   if(identity?.kind!=='exact'||identity.term!==term||!path)return true;
   return path.includes(inferred.term);
  });
  if(elements.length===concept.elements.length)return concept;
  changed=true;
  const allowed=new Set(elements);
  return {...concept,elements,...(concept.componentElements?{componentElements:concept.componentElements.filter(id=>allowed.has(id))}:{})};
 });
 return changed?result:concepts;
}
