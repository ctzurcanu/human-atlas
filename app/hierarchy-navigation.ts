import {SYSTEMS,structureName,type Atlas,type Part} from './anatomy';
import {MAJOR_SYSTEMS,REGION_ORDER,buildAnatomyNodes,depthPathFor,hierarchyEntries,majorSystemFor,regionPathFor,systemPathFor} from './anatomy-hierarchy';
import {terminologyForConcept,terminologyForGroup} from './anatomical-terminology';
import {DEPTH_LAYERS,depthLayerFor} from './depth-layers';
import {createDepthOrder} from './depth-sort';
import {guestNodeChoice} from './guest-choice';
import {resolveGuestHierarchy,type GuestHierarchy,type ResolvedGuestNode} from './guest-hierarchy';
import {anatomyNodeChoice,hierarchyChoice,type HierarchyChoice} from './hierarchy-choice';
import {withPulmonarySystem} from './pulmonary-system';
import {withHeartSystem} from './heart-system';
import {withLiverSystem} from './liver-system';
import {withTongueSystem} from './tongue-system';
import {withSpinalCordSystem} from './spinal-cord-system';

export type NavigationMode='systems'|'regions'|'depth'|`guest:${string}`;

/** Direct partonomic parent in the reviewed TA98 tree. Chapter and view
 * containers are navigation aids and must never become anatomical parents.
 * Keep the parent choice intact so selecting it selects its whole assembly,
 * including when the inspected term has no mesh in the active sex model.
 */
export function ta98PartOfChoice(choice:HierarchyChoice,ancestors:HierarchyChoice[]):HierarchyChoice|undefined{
 const term=/^hierarchy:guest:ta98:A\d{2}\.\d\.\d{2}\.\d{3}$/;
 if(!term.test(choice.id))return;
 const parent=ancestors.at(-1);
 return parent&&term.test(parent.id)&&parent.children?.some(child=>child.id===choice.id)?parent:undefined;
}

/** Build browser parentage, retaining unmodeled TA98 terms for inspection. */
export function hierarchyNavigation(atlas:Atlas,mode:NavigationMode,guest?:GuestHierarchy):HierarchyChoice{
 const available=atlas.parts.filter(part=>!part.suppressed);
 if(mode.startsWith('guest:')&&guest){
  const nodes=resolveGuestHierarchy(atlas,guest);
  return hierarchyChoice(`guest:${guest.id}:all`,'All',[...new Map(nodes.flatMap(node=>node.parts).map(part=>[part.id,part])).values()],nodes.filter(node=>guest.id==='ta98'||node.parts.length).map(node=>guestNodeChoice(guest.id,node)),terminologyForGroup('All'));
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
  const built=buildAnatomyNodes(members,path);
  const nodes=mode==='systems'&&category.id==='respiratory'?withPulmonarySystem(atlas,built):mode==='systems'&&category.id==='digestive'?withTongueSystem(atlas,withLiverSystem(atlas,built)):mode==='systems'&&category.id==='circulatory'?withHeartSystem(atlas,built):mode==='systems'&&category.id==='nervous'?withSpinalCordSystem(atlas,built):built;
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
 visited.clear();const byName=visit(root,[],equivalent);if(byName)return byName;
 // Reverse lookup for hierarchies that name the same anatomy under their own ids: the node holding exactly this
 // anatomy under any name. A node that merely contains the anatomy among others is not a parent of it.
 const sameAnatomy=(node:HierarchyChoice)=>selected.size>0&&node.id!=='selection-set'&&node.elements.length===selected.size&&node.elements.every(id=>selected.has(id));
 visited.clear();return visit(root,[],sameAnatomy)??[];
}

/** Where an item sits in a hierarchy that does not list it as a node of its own (genes, drugs, exercises…):
 * the most specific nodes whose anatomy includes all of the item's anatomy, each as a path from All to that node.
 * `more` counts the further nodes left out. */
export function hierarchyCovering(root:HierarchyChoice,choice:HierarchyChoice,limit=3):{paths:HierarchyChoice[][];more:number}{
 const selected=choice.elements;
 if(!selected.length||choice.id==='selection-set')return {paths:[],more:0};
 const sets=new WeakMap<HierarchyChoice,Set<string>>();
 const covers=(node:HierarchyChoice)=>{let set=sets.get(node);if(!set){set=new Set(node.elements);sets.set(node,set);}return selected.every(id=>set.has(id));};
 const found:{path:HierarchyChoice[];size:number}[]=[];
 const visited=new Set<HierarchyChoice>();
 const visit=(node:HierarchyChoice,path:HierarchyChoice[])=>{
  if(visited.has(node))return;visited.add(node);
  const kids=(node.children??[]).filter(covers);
  if(!kids.length){if(path.length)found.push({path:[...path,node],size:node.elements.length});return;}
  for(const kid of kids)visit(kid,[...path,node]);
 };
 if(covers(root))visit(root,[]);
 if(found.length){
  found.sort((a,b)=>a.size-b.size||b.path.length-a.path.length);
  return {paths:found.slice(0,limit).map(item=>item.path),more:Math.max(0,found.length-limit)};
 }
 // A larger item (the heart, a limb) is never wholly inside one node: list the nodes that lie within it
 // (a gene, drug or exercise tied to part of it), the ones covering most of it first.
 const inside=new Set(selected),within:{path:HierarchyChoice[];size:number}[]=[];
 const seen=new Set<HierarchyChoice>();
 const walk=(node:HierarchyChoice,path:HierarchyChoice[])=>{
  if(seen.has(node))return;seen.add(node);
  if(!node.elements.some(id=>inside.has(id)))return;
  if(path.length&&node.elements.every(id=>inside.has(id))){within.push({path:[...path,node],size:node.elements.length});return;}
  for(const kid of node.children??[])walk(kid,[...path,node]);
 };
 walk(root,[]);
 within.sort((a,b)=>b.size-a.size||a.path.length-b.path.length);
 return {paths:within.slice(0,limit).map(item=>item.path),more:Math.max(0,within.length-limit)};
}
