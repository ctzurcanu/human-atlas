import {liverCompletionTermParts,liverWholePresentationParts} from './liver-completion-controls';
import type {Atlas,Part} from './anatomy';
import type {AnatomyEntry,AnatomyNode} from './anatomy-hierarchy';
import {terminologyForConcept} from './anatomical-terminology';
import {liverControlRegistry,liverTermParts} from './liver-controls';
import {hepaticConnectionRegistry,hepaticConnectionTermParts} from './hepatic-connection-controls';
import {liverFeatureRegistry,liverSurfaceTermParts} from './liver-surface-features';
import {liverLandmarkTermParts} from './liver-landmark';
import inventory from './liver-system-inventory';

const unique=(parts:Part[])=>[...new Map(parts.map(p=>[p.id,p])).values()];
/** Systems exposes the whole 69-term organ tree. Missing terms retain empty
 * controls; surface references remain references to the explicit host facets. */
export function withLiverSystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const liver=liverControlRegistry(atlas),features=liverFeatureRegistry(atlas);if(!liver||!features)return nodes;
 const hepatic=hepaticConnectionRegistry(atlas);
 const direct=new Map([...liverTermParts(atlas),...hepaticConnectionTermParts(atlas),...liverSurfaceTermParts(atlas),...liverLandmarkTermParts(atlas),...liverCompletionTermParts(atlas)]);
 const identities=new Map(atlas.concepts.filter(c=>c.ta98Kind==='exact'&&c.ta98Term).map(c=>[c.ta98Term!,c.id]));
 const nodeFor=(term:typeof inventory.terms[number]):AnatomyNode=>{
  const children=inventory.terms.filter(t=>t.parent===term.id).map(nodeFor);
  const parts=liverWholePresentationParts(atlas,term.id,unique([...(direct.get(term.id)??[]),...children.flatMap(c=>c.parts)])),name=term.name.replace(/^./,c=>c.toUpperCase());
  const id=identities.get(term.id)??`systems:liver:${term.id}`;
  const entry:AnatomyEntry={id,name,system:'digestive',region:'Torso & pelvis',location:[],ancestry:[],parts,
   ta98Term:term.id,ta98Kind:parts.length?'exact':undefined,terminology:terminologyForConcept(id,term.id,'exact')};
  return children.length?{kind:'group',id:`liver:${term.id}`,name,nodes:children,parts,terminology:entry.terminology}
   :{kind:'entry',id:`liver:${term.id}`,name,parts,entry};
 };
 let replaced=false;
 const visit=(node:AnatomyNode):AnatomyNode=>{
  if(!replaced&&node.kind==='group'&&node.name==='Liver'&&liver.marker.partIds.every(id=>node.parts.some(p=>p.id===id))){
   replaced=true;const organ=nodeFor(inventory.terms[0])as Extract<AnatomyNode,{kind:'group'}>;
   const retained=node.nodes.filter(n=>['Hepatic arteries','Portal vein','Hepatic veins and cava'].includes(n.name));
   const conduitParts=hepatic?[...hepatic.marker.wholePartIds,...hepatic.marker.regionPartIds].map(id=>hepatic.byId.get(id)!):[];
   organ.nodes.push(...retained);organ.parts=unique([...organ.parts,...retained.flatMap(n=>n.parts),...conduitParts]);organ.itemCount=inventory.terms.length;return organ;
  }
  return node.kind==='group'?{...node,nodes:node.nodes.map(visit)}:node;
 };
 return nodes.map(visit);
}
