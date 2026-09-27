import type {Concept,Part} from './anatomy';
import type {AnatomyNode} from './anatomy-hierarchy';
import type {Terminology} from './anatomical-terminology';

export interface HierarchyChoice extends Concept {
 terminology?:Terminology;
 children?:HierarchyChoice[];
}

export function hierarchyChoice(id:string,name:string,parts:Part[],children:HierarchyChoice[],terminology?:Terminology):HierarchyChoice {
 return {id:`hierarchy:${id}`,name,elements:[...new Set(parts.map(part=>part.id))],children,terminology};
}

export function anatomyNodeChoice(node:AnatomyNode,key:string):HierarchyChoice {
 if(node.kind==='entry')return {id:node.entry.id,name:node.name,elements:node.parts.map(part=>part.id),terminology:node.entry.terminology};
 const children=node.kind==='group'
  ?node.nodes.map((child,index)=>anatomyNodeChoice(child,`${key}:${index}:${child.id}`))
  :node.entries.map(entry=>({id:entry.id,name:entry.name,elements:entry.parts.map(part=>part.id),terminology:entry.terminology}));
 return hierarchyChoice(key,node.name,node.parts,children,node.terminology);
}
