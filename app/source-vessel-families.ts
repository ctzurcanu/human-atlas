import type {Atlas,Part} from './anatomy';
import {buildAnatomyNodes,type AnatomyEntry,type AnatomyNode} from './anatomy-hierarchy';
import {catalogueTa98TermParts} from './ta98-term-bindings';
import {terminologyForConcept} from './anatomical-terminology';

type Family={conceptId:string;collectionId:string;term:string;displayName:string;partIds:string[]};
type FamilyMarker={version:1;groups:Family[]};
const unique=(parts:Part[])=>[...new Map(parts.map(p=>[p.id,p])).values()];

/** Named source families are selectors, not extra geometry or vessel walls.
 * Replace their existing selector in situ and keep each physical leaf once.
 * TA98's independent parent graph is never rewritten by this Systems view. */
export function withSourceVesselFamilies(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const marker=(atlas as Atlas&{sourceVesselFamilies?:FamilyMarker}).sourceVesselFamilies;
 if(marker?.version!==1||!Array.isArray(marker.groups)||!marker.groups.length)return nodes;
 if(marker.groups.some(g=>!g||typeof g.conceptId!=='string'||typeof g.collectionId!=='string'||typeof g.term!=='string'||typeof g.displayName!=='string'||!Array.isArray(g.partIds)||!g.partIds.length||g.partIds.some(id=>typeof id!=='string')||new Set(g.partIds).size!==g.partIds.length))return nodes;
 const concepts=new Map(atlas.concepts.map(c=>[c.id,c])),parts=new Map(atlas.parts.map(p=>[p.id,p]));
 const bindings=marker.groups.map(g=>atlas.ta98TermBindings?.find(b=>b.conceptId===g.conceptId&&b.term===g.term&&b.sourceCollectionId===g.collectionId&&b.representation==='source-collection-identity'));
 if(bindings.some(b=>!b))return nodes;
 const validated=catalogueTa98TermParts({...atlas,ta98TermBindings:bindings.filter((b):b is NonNullable<typeof b>=>!!b)});
 const allIds=marker.groups.flatMap(g=>g.partIds);
 if(new Set(allIds).size!==allIds.length||new Set(marker.groups.map(g=>g.conceptId)).size!==marker.groups.length)return nodes;
 const physical=new Map<string,AnatomyEntry>(),selectors=new Map<string,AnatomyEntry>();
 const collect=(ns:AnatomyNode[])=>{for(const n of ns){
  if(n.kind==='group')collect(n.nodes);
  else for(const e of n.kind==='entry'?[n.entry]:n.entries)(e.selectionOnly?selectors:physical).set(e.id,e);
 }};
 collect(nodes);
 const groups=new Map<string,AnatomyNode>();
 for(const g of marker.groups){
  const c=concepts.get(g.conceptId),selector=selectors.get(g.conceptId),members=g.partIds.map(id=>parts.get(id));
  const valid=validated.get(g.term)??[];
  if(!c||!selector||!g.partIds.length||!g.displayName||members.some(p=>!p||p.suppressed)||valid.length!==g.partIds.length||!valid.every(p=>g.partIds.includes(p.id))||c.elements.length!==g.partIds.length||!c.elements.every(id=>g.partIds.includes(id)))return nodes;
  const entries=[...new Set(members.map(p=>p!.conceptId))].map(id=>physical.get(id));
  if(entries.some(e=>!e)||entries.some(e=>e!.parts.some(p=>!g.partIds.includes(p.id))))return nodes;
  groups.set(g.conceptId,{kind:'group',id:'source-vessel-family:'+g.conceptId,name:g.displayName,nodes:buildAnatomyNodes(entries as AnatomyEntry[],()=>[]),parts:members as Part[],terminology:terminologyForConcept(c.id,c.ta98Term,c.ta98Kind),itemCount:entries.length+1});
 }
 const owned=new Set(allIds);
 const visit=(n:AnatomyNode):AnatomyNode[]=>{
  if(n.kind==='entry'){
   const group=groups.get(n.entry.id);if(group)return [group];
   if(n.entry.selectionOnly)return [n];
   const kept=n.parts.filter(p=>!owned.has(p.id));return kept.length?[{...n,parts:kept,entry:{...n.entry,parts:kept}}]:[];
  }
  if(n.kind==='bilateral'){
   const entries=n.entries.flatMap(e=>{const kept=e.parts.filter(p=>!owned.has(p.id));return kept.length?[{...e,parts:kept}]:[];});
   return entries.length?[{...n,entries,parts:unique(entries.flatMap(e=>e.parts))}]:[];
  }
  const children=n.nodes.flatMap(visit);return children.length?[{...n,nodes:children,parts:unique(children.flatMap(c=>c.parts))}]:[];
 };
 return nodes.flatMap(visit);
}
