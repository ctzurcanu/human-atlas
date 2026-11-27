import type {Atlas,Part} from './anatomy';
import type {AnatomyEntry,AnatomyNode} from './anatomy-hierarchy';
import {terminologyForConcept} from './anatomical-terminology';
import {cardiacIntegration} from './heart-wall';
import {cardiacLabelRegistry} from './cardiac-label-sites';
import {hasSourceCardiacIntegration,sourceCardiacIntegration} from './biodigital-heart';
const unique=(parts:Part[])=>[...new Map(parts.map(p=>[p.id,p])).values()];

function withoutOwnedParts(nodes:AnatomyNode[],owned:Set<string>):AnatomyNode[]{
 const prune=(node:AnatomyNode):AnatomyNode[]=>{
  if(node.kind==='entry'){const parts=node.parts.filter(p=>!owned.has(p.id));return parts.length?[{...node,parts,entry:{...node.entry,parts}}]:[];}
  if(node.kind==='bilateral'){const entries=node.entries.flatMap(e=>{const parts=e.parts.filter(p=>!owned.has(p.id));return parts.length?[{...e,parts}]:[];});return entries.length?[{...node,entries,parts:unique(entries.flatMap(e=>e.parts))}]:[];}
  const children=node.nodes.flatMap(prune);return children.length?[{...node,nodes:children,parts:unique(children.flatMap(child=>child.parts))}]:[];
 };
 return nodes.flatMap(prune).map(node=>node.kind==='group'&&node.name==='Heart'?{...node,name:'Other cardiac structures'}:node);
}

/** Systems groups the installed organ independently of canonical TA98 parents. */
export function withHeartSystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 if(hasSourceCardiacIntegration(atlas)){
  const source=sourceCardiacIntegration(atlas);if(!source)return nodes;
  const whole=source.termGroups.find(g=>g.term==='A12.1.00.001');if(!whole)return nodes;
  const available=new Map(atlas.parts.map(p=>[p.id,p])),owned=new Set(source.sourceParts.map(p=>p.partId));
  // These are explicit selectors over original material, not additional meshes.
  // The opening selects its active carrier; its hidden source patch stays out.
  const controls:AnatomyNode[]=source.termGroups.filter(g=>g!==whole).map(g=>{
   const parts=g.partIds.map(id=>available.get(id)!),terminology=terminologyForConcept(g.conceptId,g.term,g.kind);
   const entry:AnatomyEntry={id:g.conceptId,name:g.name,system:'cardiac',region:'Torso & pelvis',location:[],ancestry:[],parts,ta98Term:g.term,ta98Kind:g.kind,terminology};
   return {kind:'entry',id:'cardiac:source:'+g.key,name:g.name,parts,entry};
  });
  const heart:AnatomyNode={kind:'group',id:'cardiac:source:heart',name:'Heart',nodes:controls,parts:unique(whole.partIds.map(id=>available.get(id)!)),terminology:terminologyForConcept(whole.conceptId,whole.term,whole.kind),itemCount:source.termGroups.length};
  return [...withoutOwnedParts(nodes,owned),heart].sort((a,b)=>a.name.localeCompare(b.name));
 }
 const marker=cardiacIntegration(atlas);if(!marker?.termGroups)return nodes;
 const termGroups=[...marker.termGroups,...cardiacLabelRegistry(atlas)?.termGroups??[]];
 const available=new Map(atlas.parts.filter(p=>!p.suppressed).map(p=>[p.id,p]));
 const owned=new Set([...marker.wallPartIds,...marker.valvePartIds,...marker.papillaryPartIds,...marker.chordalPartIds??[],...marker.rootPartIds??[]]);
 if([...owned].some(id=>!available.has(id)))return nodes;
 const nodeFor=(key:string):AnatomyNode=>{
  const g=termGroups.find(row=>row.key===key)!,children=termGroups.filter(row=>row.parent===key).map(row=>nodeFor(row.key));
  const parts=unique([...g.partIds.map(id=>available.get(id)!).filter(Boolean),...children.flatMap(child=>child.parts)]),terminology=terminologyForConcept(g.conceptId,g.term,g.kind);
  const entry:AnatomyEntry={id:g.conceptId,name:g.name,system:'cardiac',region:'Torso & pelvis',location:[],ancestry:[],parts,ta98Term:g.term,ta98Kind:g.kind,terminology};
  const terms=new Set<string>(),pending=[g];
  while(pending.length){const row=pending.pop()!;if(row.kind==='exact')terms.add(row.term);pending.push(...termGroups.filter(child=>child.parent===row.key));}
  return children.length?{kind:'group',id:'cardiac:'+key,name:g.name,nodes:children,parts,terminology,itemCount:terms.size}:{kind:'entry',id:'cardiac:'+key,name:g.name,parts,entry};
 };
 if(!marker.termGroups.some(g=>g.key==='heart'))return nodes;
 const heart=nodeFor('heart');
 if(heart.kind==='group'&&marker.rootPartIds?.length){
  const roots=marker.rootPartIds.map(id=>available.get(id)!);
  heart.nodes.push({kind:'group',id:'cardiac:vessel-roots',name:'Vessel roots',parts:roots,terminology:terminologyForConcept(''),nodes:roots.map(p=>({kind:'entry' as const,id:'cardiac:root:'+p.id,name:p.name,parts:[p],entry:{id:p.conceptId,name:p.name,system:p.system,region:'Torso & pelvis',location:[],ancestry:[],parts:[p],terminology:terminologyForConcept(p.conceptId,p.ta98Term)}}))});
  heart.parts=unique([...heart.parts,...roots]);
 }
 const rest=withoutOwnedParts(nodes,owned);
 return [...rest,heart].sort((a,b)=>a.name.localeCompare(b.name));
}
