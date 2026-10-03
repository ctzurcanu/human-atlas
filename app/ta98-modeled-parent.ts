import {nonComponentRelationship} from './anatomical-containment';
import index from './data/ta98-parents.json';
import {taExactEntityForConcept} from './anatomical-terminology';
import type {Concept,Part} from './anatomy';
import {displayLaterality} from './laterality';
import type {HierarchyChoice} from './hierarchy-choice';
const parents=index.parents as Record<string,string>;
const exactTerm=(concept:Concept)=>concept.ta98Kind==='parent'?undefined:concept.ta98Term??taExactEntityForConcept(concept.id);
const sides=(concept:Concept,parts:Map<string,Part>)=>{const named=displayLaterality(concept.name).side;return named?new Set([named]):new Set(concept.elements.flatMap(id=>{const part=parts.get(id),side=part&&displayLaterality(part.name).side;return side?[side]:[];}));};
const identities=new WeakMap<Concept[],Map<string,Concept>>();
const conceptIndex=new WeakMap<Concept[],Map<string,Concept[]>>();
function conceptsByTerm(concepts:Concept[]):Map<string,Concept[]>{
 const cached=conceptIndex.get(concepts);if(cached)return cached;
 const byTerm=new Map<string,Concept[]>();
 for(const concept of concepts){
  const code=exactTerm(concept);if(!code)continue;
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
 let byId=identities.get(concepts);if(!byId){byId=new Map(concepts.map(concept=>[concept.id,concept]));identities.set(concepts,byId);}
 const identity=choice.ta98Term||choice.ta98Kind?choice:byId.get(choice.id)??choice;
 const code=exactTerm(identity),parent=code&&parents[code];
 if(!parent)return;
 const members=choice.elements.map(id=>parts.get(id)).filter((part):part is Part=>!!part&&!part.suppressed);
 if(!members.length)return;
 const systems=new Set(members.map(part=>part.system)),sourceSides=sides(choice,parts);
 const candidates=(conceptsByTerm(concepts).get(parent)??[]).filter(concept=>concept.id!==choice.id).flatMap(concept=>{
  if(nonComponentRelationship(identity.name,concept.name))return [];
  const elements=concept.elements.filter(id=>{const part=parts.get(id);return part&&!part.suppressed;});
  if(!elements.length||!elements.some(id=>systems.has(parts.get(id)!.system)))return [];
  // A shared TA98 parent is not evidence that separately sourced organs fit.
  // Stored assemblies containing the child remain valid across source IDs.
  if(!members.every(part=>elements.includes(part.id))){
   const frame=nativePartFrame(members[0]);
   if(!frame||!members.every(part=>nativePartFrame(part)===frame)||!elements.every(id=>nativePartFrame(parts.get(id)!)===frame))return [];
  }
  const parentSides=sides({...concept,elements},parts);
  if(parentSides.size&&[...sourceSides].some(side=>!parentSides.has(side)))return [];
  return [{...concept,elements}];
 });
 // Prefer the actual stored assembly containing the selected region. Another
 // source organ with the same term must not displace that coordinated assembly.
 candidates.sort((a,b)=>Number(members.every(part=>b.elements.includes(part.id)))-Number(members.every(part=>a.elements.includes(part.id))));
 return candidates[0];
}

/** Known native catalogue frames only. HRA donors are independently registered
 * reference organs, so a shared HRA prefix does not establish a common frame. */
export function nativePartFrame(part:Part):string|undefined{
 return /^LOCAL:(?:male|female):/.exec(part.id)?.[0]??(part.id.startsWith('ZA:')?'ZA:':undefined);
}
const componentIndex=new WeakMap<Concept[],WeakMap<Map<string,Part>,Map<string,string[]>>>();
/** Invert the exact direct-parent resolver. Shared terminology alone cannot
 * bring an independently registered donor into a native organ assembly. */
export function ta98ModeledComponentIds(choice:HierarchyChoice,concepts:Concept[],parts:Map<string,Part>):string[]{
 if(choice.id==='selection-set'||choice.children)return [];
 let byParts=componentIndex.get(concepts);if(!byParts){byParts=new WeakMap();componentIndex.set(concepts,byParts);}
 let graph=byParts.get(parts);if(!graph){
  graph=new Map();
  for(const child of concepts){
   const parent=ta98ModeledParent(child,concepts,parts);if(!parent)continue;
   const parentIds=new Set(parent.elements),frames=new Set(parent.elements.flatMap(id=>{const part=parts.get(id),frame=part&&nativePartFrame(part);return frame?[frame]:[];}));
   const ids=child.elements.filter(id=>{const part=parts.get(id);if(!part||part.suppressed)return false;if(parentIds.has(id))return true;const frame=nativePartFrame(part);return !!frame&&frames.has(frame);});
   if(!ids.length)continue;const prior=graph.get(parent.id)??[];graph.set(parent.id,[...new Set([...prior,...ids])]);
  }
  byParts.set(parts,graph);
 }
 return graph.get(choice.id)??[];
}
