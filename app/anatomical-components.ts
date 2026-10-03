import {nonComponentRelationship} from './anatomical-containment';
import {nativePartFrame,ta98ModeledComponentIds} from './ta98-modeled-parent';
import type {Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';

/** Explicit catalogue assemblies must actually contain every selected piece. */
export function storedAssemblyParent(choice:HierarchyChoice,concepts:Concept[],parts:Map<string,Part>):HierarchyChoice|undefined{
 if(choice.children||choice.id==='selection-set')return;
 const identity=concepts.find(concept=>concept.id===choice.id);
 const parent=identity?.partOf&&concepts.find(concept=>concept.id===identity.partOf);
 if(!parent||parent.id===choice.id)return;
 const elements=parent.elements.filter(id=>parts.has(id)&&!parts.get(id)!.suppressed);
 if(!choice.elements.length||!choice.elements.every(id=>elements.includes(id)))return;
 return {...parent,elements};
}

/** Parent navigation is separate from the related-structure expansion action. */
export function hasExpandableAnatomicalRelations(relations:ResolvedRelation[],parts:Map<string,Part>):boolean{
 return relations.some(relation=>relation.kind!=='partOf'&&!!parts.get(relation.target.id)&&!parts.get(relation.target.id)!.suppressed);
}

/** Use component membership, never adjacency, supply or spatial overlap.
 * TA98 nodes also carry a direct source-organ leaf; that leaf is the whole,
 * not another component, and must be excluded from the component view.
 */
export function anatomicalComponentIds(choice:HierarchyChoice,relations:ResolvedRelation[],parts:Map<string,Part>,concepts:Concept[]=[]):string[]{
 const ta98=choice.id.startsWith('hierarchy:guest:ta98:');
 const children=(choice.children??[]).filter(child=>!ta98||/^hierarchy:guest:ta98:A\d{2}\.\d\.\d{2}\.\d{3}$/.test(child.id));
 const members=choice.elements.map(id=>parts.get(id)).filter((part):part is Part=>!!part&&!part.suppressed);
 // A stored, named multi-part assembly can be selected directly, outside its
 // hierarchy node. Distinct member concepts establish subdivisions; several
 // chunks of a single concept do not. Mixed ad-hoc selections have no parent.
 const assembly=!choice.children&&choice.id!=='selection-set'&&new Set(members.map(part=>part.conceptId)).size>1
  ?members.filter(part=>part.conceptId!==choice.id).map(part=>part.id):[];
 const frames=new Set(members.flatMap(part=>{const frame=nativePartFrame(part);return frame?[frame]:[];}));
 const contains=!choice.children?relations.filter(relation=>relation.kind==='contains'&&(!frames.size||members.some(part=>part.id===relation.target.id)||frames.has(nativePartFrame(relation.target)??''))).map(relation=>relation.target.id):[];
 const ids=[...ta98ModeledComponentIds(choice,concepts,parts),...children.flatMap(child=>child.elements),...contains,...assembly];
 return [...new Set(ids)].filter(id=>{const part=parts.get(id);return part&&!part.suppressed&&!nonComponentRelationship(part.name,choice.name);});
}

/** The inspector and relationship expansion share the same component evidence.
 * Reconcile lexical Contains links against native-frame membership, then add
 * exact TA98 children and explicit stored subdivisions without duplicating them. */
export function anatomicalComponentRelations(choice:HierarchyChoice,relations:ResolvedRelation[],parts:Map<string,Part>,concepts:Concept[]):ResolvedRelation[]{
 if(choice.children||choice.id==='selection-set')return relations;
 const ids=anatomicalComponentIds(choice,relations,parts,concepts),allowed=new Set(ids);
 const result=relations.filter(relation=>relation.kind!=='contains'||allowed.has(relation.target.id));
 const present=new Set(result.filter(relation=>relation.kind==='contains').map(relation=>relation.target.id));
 for(const id of ids){if(present.has(id))continue;const target=parts.get(id);if(target){result.push({kind:'contains',target});present.add(id);}}
 return result;
}
