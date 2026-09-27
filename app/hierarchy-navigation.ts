import {SYSTEMS,structureName,type Atlas,type Part} from './anatomy';
import {MAJOR_SYSTEMS,REGION_ORDER,buildAnatomyNodes,depthPathFor,hierarchyEntries,majorSystemFor,regionPathFor,systemPathFor} from './anatomy-hierarchy';
import {terminologyForConcept,terminologyForGroup} from './anatomical-terminology';
import {DEPTH_LAYERS,depthLayerFor} from './depth-layers';
import {createDepthOrder} from './depth-sort';
import {resolveGuestHierarchy,type GuestHierarchy,type ResolvedGuestNode} from './guest-hierarchy';
import {anatomyNodeChoice,hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';

export type NavigationMode='systems'|'regions'|'depth'|`guest:${string}`;

/** Build the same modeled parentage used by the Layers browser for details navigation. */
export function hierarchyNavigation(atlas:Atlas,mode:NavigationMode,guest?:GuestHierarchy):HierarchyChoice{
 const available=atlas.parts.filter(part=>!part.suppressed);
 if(mode.startsWith('guest:')&&guest){
  const guestChoice=(node:ResolvedGuestNode):HierarchyChoice=>{
   const children=node.children.map(guestChoice);
   const inherited=new Set(node.children.flatMap(child=>child.parts.map(part=>part.id)));
   const direct=new Map<string,Part[]>();
   for(const part of node.directParts)if(!inherited.has(part.id)){const list=direct.get(part.conceptId)??[];list.push(part);direct.set(part.conceptId,list);}
   const leaves=[...direct].map(([id,parts])=>({id,name:structureName(parts[0].name),elements:parts.map(part=>part.id),terminology:terminologyForConcept(id)}));
   return hierarchyChoice(`guest:${guest.id}:${node.id}`,node.name,node.parts,[...children,...leaves],terminologyForGroup(node.name));
  };
  const nodes=resolveGuestHierarchy(atlas,guest);
  return hierarchyChoice(`guest:${guest.id}:all`,'All',available,nodes.map(guestChoice),terminologyForGroup('All'));
 }
 const entries=hierarchyEntries(atlas);
 if(mode==='depth'){
  const compare=createDepthOrder(atlas);
  const layers=DEPTH_LAYERS.flatMap(layer=>{
   const members=entries.flatMap(entry=>{const parts=entry.parts.filter(part=>depthLayerFor(part)===layer.id);return parts.length?[{...entry,parts}]:[];});
   if(!members.length)return [];
   const nodes=buildAnatomyNodes(members,entry=>depthPathFor(entry,layer.name),compare);
   return [hierarchyChoice(`depth:${layer.id}`,layer.name,members.flatMap(entry=>entry.parts),nodes.map((node,index)=>anatomyNodeChoice(node,`depth:${layer.id}:${index}`)),terminologyForGroup(layer.name))];
  });
  return hierarchyChoice('depth:all','All',available,layers,terminologyForGroup('All'));
 }
 const categories=mode==='systems'
  ?atlas.scope==='cell'?SYSTEMS:MAJOR_SYSTEMS.map(system=>({id:system.id,name:system.name}))
  :(atlas.scope==='cell'?['Cell boundary','Nucleus','Cytoplasm']:REGION_ORDER).map(name=>({id:name,name}));
 const branches=categories.flatMap(category=>{
  const members=entries.filter(entry=>mode==='systems'?(atlas.scope==='cell'?entry.system===category.id:majorSystemFor(entry)===category.id):entry.region===category.id);
  if(!members.length)return [];
  const path=(entry:typeof entries[number])=>mode==='systems'?systemPathFor(entry,atlas.scope):regionPathFor(entry);
  const nodes=buildAnatomyNodes(members,path);
  return [hierarchyChoice(`${mode}:${category.id}`,category.name,members.flatMap(entry=>entry.parts),nodes.map((node,index)=>anatomyNodeChoice(node,`${mode}:${category.id}:${index}`)),terminologyForGroup(category.name))];
 });
 return hierarchyChoice(`${mode}:all`,'All',available,branches,terminologyForGroup('All'));
}

/** The returned choices are selectable ancestors, from All to the immediate parent. */
export function hierarchyAncestors(root:HierarchyChoice,choice:HierarchyChoice):HierarchyChoice[]{
 if(choice.id==='selection-set')return [];
 const selected=new Set(choice.elements);
 const exact=(node:HierarchyChoice)=>node.id===choice.id&&(!!node.children?node.elements.length===choice.elements.length&&node.elements.every(id=>selected.has(id)):node.elements.some(id=>selected.has(id)));
 const equivalent=(node:HierarchyChoice)=>!!choice.children&&!!node.children&&node.name===choice.name&&node.elements.length===choice.elements.length&&node.elements.every(id=>selected.has(id));
 const visit=(node:HierarchyChoice,path:HierarchyChoice[]):HierarchyChoice[]|null=>{
  if(exact(node)||equivalent(node))return path;
  if(!node.children)return null;
  if(selected.size&&!choice.children&&!node.elements.some(id=>selected.has(id)))return null;
  for(const child of node.children){const result=visit(child,[...path,node]);if(result)return result;}
  return null;
 };
 return visit(root,[])??[];
}
