import type {Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {ResolvedRelation} from './anatomical-relations';

/** Parent navigation is separate from the related-structure expansion action. */
export function hasExpandableAnatomicalRelations(relations:ResolvedRelation[],parts:Map<string,Part>):boolean{
 return relations.some(relation=>relation.kind!=='partOf'&&!!parts.get(relation.target.id)&&!parts.get(relation.target.id)!.suppressed);
}

/** Use component membership, never adjacency, supply or spatial overlap.
 * TA98 nodes also carry a direct source-organ leaf; that leaf is the whole,
 * not another component, and must be excluded from the component view.
 */
export function anatomicalComponentIds(choice:HierarchyChoice,relations:ResolvedRelation[],parts:Map<string,Part>):string[]{
 const ta98=choice.id.startsWith('hierarchy:guest:ta98:');
 const children=(choice.children??[]).filter(child=>!ta98||/^hierarchy:guest:ta98:A\d{2}\.\d\.\d{2}\.\d{3}$/.test(child.id));
 const members=choice.elements.map(id=>parts.get(id)).filter((part):part is Part=>!!part&&!part.suppressed);
 // A stored, named multi-part assembly can be selected directly, outside its
 // hierarchy node. Distinct member concepts establish subdivisions; several
 // chunks of a single concept do not. Mixed ad-hoc selections have no parent.
 const assembly=!choice.children&&choice.id!=='selection-set'&&new Set(members.map(part=>part.conceptId)).size>1
  ?members.filter(part=>part.conceptId!==choice.id).map(part=>part.id):[];
 const ids=[...children.flatMap(child=>child.elements),...(!choice.children?relations.filter(relation=>relation.kind==='contains').map(relation=>relation.target.id):[]),...assembly];
 return [...new Set(ids)].filter(id=>{const part=parts.get(id);return part&&!part.suppressed;});
}
