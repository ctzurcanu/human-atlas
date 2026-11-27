import type {Atlas,Part} from './anatomy';
import type {AnatomyEntry,AnatomyNode} from './anatomy-hierarchy';
import {taExactEntityForConcept,terminologyForConcept} from './anatomical-terminology';
import inventory from './pulmonary-system-inventory';

export const pulmonarySystemInventory=inventory;

const whole='A06.5.01.001';
const terms=new Set(inventory.terms.map(term=>term.id));
const unique=(parts:Part[])=>[...new Map(parts.map(part=>[part.id,part])).values()];
type PulmonaryAtlas=Atlas&{
 pulmonarySegmentIntegration?:{schema?:string;partIds?:string[]};
 pulmonaryBranchIntegration?:{schema?:string;partIds?:string[]};
};

/** Complete organ navigation in Systems, independent of the canonical TA98 tree.
 * Empty terms stay visible without borrowing a whole lung as their geometry.
 * Repeated collective controls reference existing part IDs, never extra meshes.
 */
export function withPulmonarySystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const marker=(atlas as PulmonaryAtlas).pulmonarySegmentIntegration;
 if(marker?.schema!=='owned-bronchopulmonary-territories/v1'||!Array.isArray(marker.partIds))return nodes;
 const hidden=new Set(atlas.defaultHidden??[]);
 const available=new Map(atlas.parts.filter(part=>!part.suppressed&&!part.defaultHidden&&!hidden.has(part.id)).map(part=>[part.id,part]));
 if(marker.partIds.length!==19||marker.partIds.some(id=>!available.has(id)))return nodes;
 const direct=new Map<string,Part[]>(),identities=new Map<string,string>();
 for(const concept of atlas.concepts){
  if(concept.ta98Kind==='parent')continue;
  const term=concept.ta98Term??taExactEntityForConcept(concept.id);
  if(!term||!terms.has(term))continue;
  const parts=concept.elements.flatMap(id=>available.has(id)?[available.get(id)!]:[]);
  if(!parts.length)continue;
  direct.set(term,unique([...(direct.get(term)??[]),...parts]));
  if(!identities.has(term))identities.set(term,concept.id);
 }
 // Collective gross terms refer to their actual named lobes/lingular regions.
 for(const [term,children]of [
  ['A06.5.01.017',['A06.5.02.002','A06.5.02.015']],
  ['A06.5.01.018',['A06.5.02.018','A06.5.02.019']],
  ['A06.5.01.019',['A06.5.02.006']],
  ['A06.5.01.020',['A06.5.02.009','A06.5.02.020']],
 ] as const)direct.set(term,unique(children.flatMap(id=>direct.get(id)??[])));
 direct.set('A06.5.02.001',marker.partIds.map(id=>available.get(id)!));
 const owned=new Set([...direct.values()].flatMap(parts=>parts.map(part=>part.id)));
 const branchMarker=(atlas as PulmonaryAtlas).pulmonaryBranchIntegration;
 const branches=branchMarker?.schema==='continuous-pulmonary-source/v1'&&Array.isArray(branchMarker.partIds)
  ?branchMarker.partIds.flatMap(id=>available.has(id)?[available.get(id)!]:[]):[];
 branches.forEach(part=>owned.add(part.id));
 const entry=(term:string,name:string,parts:Part[],id=identities.get(term)??`systems:pulmonary:${term}`):AnatomyEntry=>({
  id,name,system:'respiratory',region:'Torso & pelvis',location:[],ancestry:[],parts,
  ta98Term:term,ta98Kind:parts.length?'exact':undefined,terminology:terminologyForConcept(id,term,'exact'),
 });
 const nodeFor=(term:typeof inventory.terms[number]):AnatomyNode=>{
  const children=inventory.terms.filter(child=>child.parent===term.id).map(nodeFor);
  const parts=unique([...(direct.get(term.id)??[]),...children.flatMap(child=>child.parts)]);
  const name=term.name.replace(/^./,letter=>letter.toUpperCase());
  return children.length?{kind:'group',id:`pulmonary:${term.id}`,name,nodes:children,parts,terminology:terminologyForConcept('',term.id,'exact')}
   :{kind:'entry',id:`pulmonary:${term.id}`,name,parts,entry:entry(term.id,name,parts)};
 };
 const lung=nodeFor(inventory.terms.find(term=>term.id===whole)!) as Extract<AnatomyNode,{kind:'group'}>;
 lung.itemCount=inventory.terms.length;
 if(branches.length){
  const bronchial=branches.map(part=>({kind:'entry' as const,id:`pulmonary:bronchus:${part.id}`,name:part.name,parts:[part],
   entry:entry(part.ta98Term??'A06.4.01.002',part.name,[part],part.conceptId)}));
  lung.nodes.push({kind:'group',id:'pulmonary:bronchial-tree',name:'Bronchial tree',nodes:bronchial,parts:branches,terminology:terminologyForConcept('','A06.4.01.002','exact')});
 }
 // Remove the old scattered occurrences from Airways, source folders and
 // region paths while retaining every unrelated structure in those branches.
 const prune=(node:AnatomyNode):AnatomyNode[]=>{
  if(node.kind==='entry'){
   const parts=node.parts.filter(part=>!owned.has(part.id));
   return parts.length?[{...node,parts,entry:{...node.entry,parts}}]:[];
  }
  if(node.kind==='bilateral'){
   const entries=node.entries.flatMap(member=>{const parts=member.parts.filter(part=>!owned.has(part.id));return parts.length?[{...member,parts}]:[];});
   return entries.length?[{...node,entries,parts:unique(entries.flatMap(member=>member.parts))}]:[];
  }
  const children=node.nodes.flatMap(prune);
  return children.length?[{...node,nodes:children,parts:unique(children.flatMap(child=>child.parts))}]:[];
 };
 const rest=nodes.flatMap(prune),previous=rest.find(node=>node.kind==='group'&&node.name==='Lungs');
 if(previous?.kind==='group')lung.nodes.push(...previous.nodes);
 lung.parts=unique(lung.nodes.flatMap(node=>node.parts).concat(direct.get(whole)??[]));
 return [...rest.filter(node=>node!==previous),lung].sort((a,b)=>a.name.localeCompare(b.name));
}
