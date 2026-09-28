import {structureName,type Part} from './anatomy';
import {terminologyForConcept,terminologyForGroup} from './anatomical-terminology';
import {hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';
import type {ResolvedGuestNode} from './guest-hierarchy';

const cache=new WeakMap<ResolvedGuestNode,Map<string,HierarchyChoice>>();
/** A substance and its single identically named drug record need only one row. */
export function guestVisibleChildren(hierarchyId:string,node:ResolvedGuestNode):ResolvedGuestNode[]{
 const drug=node.children[0];
 return hierarchyId==='drugs'&&/^ATC:[A-Z]\d{2}[A-Z]{2}\d{2}$/.test(node.id)&&node.children.length===1&&/^CHEMBL\d+$/.test(drug.id)&&!drug.directParts.length&&node.name.toLowerCase()===drug.name.toLowerCase()?drug.children:node.children;
}
/** Reuse DAG choices and omit unmodeled descendants from selection navigation. */
export function guestNodeChoice(hierarchyId:string,node:ResolvedGuestNode):HierarchyChoice{
 let choices=cache.get(node);if(!choices){choices=new Map();cache.set(node,choices);}
 const existing=choices.get(hierarchyId);if(existing)return existing;
 const children=guestVisibleChildren(hierarchyId,node),nested=children.filter(child=>child.parts.length).map(child=>guestNodeChoice(hierarchyId,child));
 const inherited=new Set(children.flatMap(child=>child.parts.map(part=>part.id))),direct=new Map<string,Part[]>();
 for(const part of node.directParts)if(!inherited.has(part.id)){const parts=direct.get(part.conceptId)??[];parts.push(part);direct.set(part.conceptId,parts);}
 const leaves=[...direct].map(([id,parts])=>({id,name:structureName(parts[0].name),elements:parts.map(part=>part.id),terminology:terminologyForConcept(id)}));
 const result=hierarchyChoice(`guest:${hierarchyId}:${node.id}`,node.name,node.parts,[...nested,...leaves],terminologyForGroup(node.name));
 choices.set(hierarchyId,result);return result;
}
