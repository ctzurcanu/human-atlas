import index from './data/ta98-parents.json';
import {taExactEntityForConcept} from './anatomical-terminology';
import type {Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
const parents=index.parents as Record<string,string>;
const conceptIndex=new WeakMap<Concept[],Map<string,Concept[]>>();
function conceptsByTerm(concepts:Concept[]):Map<string,Concept[]>{
 const cached=conceptIndex.get(concepts);if(cached)return cached;
 const byTerm=new Map<string,Concept[]>();
 for(const concept of concepts){
  const code=concept.ta98Term??taExactEntityForConcept(concept.id);if(!code)continue;
  const group=byTerm.get(code)??[];group.push(concept);byTerm.set(code,group);
 }
 conceptIndex.set(concepts,byTerm);return byTerm;
}

/** Resolve a direct TA98 parent to its stored whole-organ assembly outside
 * the TA98 browser. Do not invent a parent from a name suffix or spatial bounds,
 * skip an absent parent, or mistake a placement mapping for exact identity.
 */
export function ta98ModeledParent(choice:HierarchyChoice,concepts:Concept[],parts:Map<string,Part>):HierarchyChoice|undefined{
 if(choice.id==='selection-set'||choice.children)return;
 const code=choice.ta98Term??taExactEntityForConcept(choice.id),parent=code&&parents[code];
 if(!parent)return;
 const members=choice.elements.map(id=>parts.get(id)).filter((part):part is Part=>!!part&&!part.suppressed);
 if(!members.length)return;
 const systems=new Set(members.map(part=>part.system));
 const candidates=(conceptsByTerm(concepts).get(parent)??[]).filter(concept=>concept.id!==choice.id).flatMap(concept=>{
  const elements=concept.elements.filter(id=>{const part=parts.get(id);return part&&!part.suppressed;});
  if(!elements.length||!elements.some(id=>systems.has(parts.get(id)!.system)))return [];
  return [{...concept,elements}];
 });
 // Prefer the actual stored assembly containing the selected region. Another
 // source organ with the same term must not displace that coordinated assembly.
 candidates.sort((a,b)=>Number(members.every(part=>b.elements.includes(part.id)))-Number(members.every(part=>a.elements.includes(part.id))));
 return candidates[0];
}
