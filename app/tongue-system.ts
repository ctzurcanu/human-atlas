import type {Atlas,Part} from './anatomy';
import type {AnatomyEntry,AnatomyNode} from './anatomy-hierarchy';
import {terminologyForConcept} from './anatomical-terminology';
import {tongueRegistry,type TongueTermGroup} from './tongue-controls';
const unique=(parts:Part[])=>[...new Map(parts.map(p=>[p.id,p])).values()];

/** The Systems organ assembly is independent of canonical TA98 chapter parents. */
export function withTongueSystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const registry=tongueRegistry(atlas);if(!registry)return nodes;
 const groups=registry.marker.termGroups,root=groups.find(g=>g.conceptId===registry.marker.wholeConceptId)!;
 const build=(g:TongueTermGroup):AnatomyNode=>{
  const children=groups.filter(child=>child.parent===g.key||g===root&&!child.parent&&child!==root).map(build);
  const parts=unique([...g.partIds.map(id=>registry.byId.get(id)!),...children.flatMap(n=>n.parts)]),terminology=terminologyForConcept(g.conceptId,g.term,'exact');
  const entry:AnatomyEntry={id:g.conceptId,name:g.name,system:'digestive',region:'Head & neck',location:[],ancestry:[],parts,ta98Term:g.term,ta98Kind:'exact',terminology};
  return children.length?{kind:'group',id:'tongue:'+g.key,name:g.name,nodes:children,parts,terminology}:{kind:'entry',id:'tongue:'+g.key,name:g.name,parts,entry};
 };
 let replaced=false;
 const visit=(node:AnatomyNode):AnatomyNode=>{
  if(!replaced&&node.name.toLowerCase()==='tongue'&&registry.parts.every(p=>node.parts.some(q=>q.id===p.id))){replaced=true;return build(root);}
  return node.kind==='group'?{...node,nodes:node.nodes.map(visit)}:node;
 };
 const result=nodes.map(visit);
 return replaced?result:[...result,build(root)];
}
