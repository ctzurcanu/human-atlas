import {SYSTEMS,structureName,type Atlas,type Part} from './anatomy';
import {MAJOR_SYSTEMS,REGION_ORDER,buildAnatomyNodes,depthPathFor,hierarchyEntries,majorSystemFor,regionPathFor,systemPathFor} from './anatomy-hierarchy';
import {terminologyForConcept,terminologyForGroup} from './anatomical-terminology';
import {DEPTH_LAYERS,depthLayerFor} from './depth-layers';
import {createDepthOrder} from './depth-sort';
import {guestNodeChoice} from './guest-choice';
import {resolveGuestHierarchy,type GuestHierarchy,type ResolvedGuestNode} from './guest-hierarchy';
import {anatomyNodeChoice,hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';

export type NavigationMode='systems'|'regions'|'depth'|`guest:${string}`;

/** Build the same modeled parentage used by the Layers browser for details navigation. */
export function hierarchyNavigation(atlas:Atlas,mode:NavigationMode,guest?:GuestHierarchy):HierarchyChoice{
 const available=atlas.parts.filter(part=>!part.suppressed);
 if(mode.startsWith('guest:')&&guest){
  const nodes=resolveGuestHierarchy(atlas,guest);
  return hierarchyChoice(`guest:${guest.id}:all`,'All',[...new Map(nodes.flatMap(node=>node.parts).map(part=>[part.id,part])).values()],nodes.filter(node=>node.parts.length).map(node=>guestNodeChoice(guest.id,node)),terminologyForGroup('All'));
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
 const visited=new Set<HierarchyChoice>();
 const visit=(node:HierarchyChoice,path:HierarchyChoice[],match:(node:HierarchyChoice)=>boolean):HierarchyChoice[]|null=>{
  if(visited.has(node))return null;visited.add(node);
  if(match(node))return path;
  if(!node.children)return null;
  if(selected.size&&!choice.children&&!node.elements.some(id=>selected.has(id)))return null;
  for(const child of node.children){const result=visit(child,[...path,node],match);if(result)return result;}
  return null;
 };
 const exactPath=visit(root,[],exact);if(exactPath)return exactPath;
 visited.clear();return visit(root,[],equivalent)??[];
}
