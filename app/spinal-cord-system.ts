import type {Atlas,Part} from './anatomy';
import type {AnatomyEntry,AnatomyNode} from './anatomy-hierarchy';
import {terminologyForConcept} from './anatomical-terminology';
import {spinalCordRegistry,spinalCordTermParts,spinalCordSegments} from './spinal-cord-controls';
import inventory from './spinal-cord-system-inventory';

const unique=(parts:Part[])=>[...new Map(parts.map(p=>[p.id,p])).values()];
/** Systems exposes named regional controls as well as the physical tissues.
 * Unconstructed terms remain disabled; a parent selection is not coverage.
 */
export function withSpinalCordSystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const cord=spinalCordRegistry(atlas);if(!cord)return nodes;
 const direct=spinalCordTermParts(atlas),identities=new Map(cord.marker.termGroups.map(g=>[g.term,g.conceptId]));
 const byId=new Map(cord.parts.map(p=>[p.id,p])),segments=spinalCordSegments(atlas);
 const cordTerms=inventory.terms.filter(t=>t.id.startsWith('A14.1.02.'));
 const nodeFor=(row:typeof inventory.terms[number],rows:readonly typeof row[]):AnatomyNode=>{
  const children=rows.filter(t=>t.parent===row.id).map(t=>nodeFor(t,rows));
  for(const segment of segments.filter(s=>s.parentTerm===row.id)){
   const parts=segment.partIds.map(id=>byId.get(id)!);
   const entry:AnatomyEntry={id:segment.conceptId,name:`${segment.name} spinal cord segment`,system:'nervous',region:'Torso & pelvis',location:[],ancestry:[],parts,terminology:terminologyForConcept(segment.conceptId)};
   children.push({kind:'entry',id:segment.conceptId,name:entry.name,entry,parts});
  }
  const parts=unique([...(direct.get(row.id)??[]),...children.flatMap(c=>c.parts)]),name=row.name.replace(/^./,c=>c.toUpperCase()),id=identities.get(row.id)??`systems:spinal:${row.id}`;
  const entry:AnatomyEntry={id,name,system:'nervous',region:'Torso & pelvis',location:[],ancestry:[],parts,ta98Term:row.id,
   ta98Kind:direct.has(row.id)?'exact':undefined,terminology:terminologyForConcept(id,row.id,direct.has(row.id)?'exact':undefined)};
  return children.length?{kind:'group',id:'spinal:'+row.id,name,nodes:children,parts,terminology:entry.terminology}
   :{kind:'entry',id:'spinal:'+row.id,name,entry,parts};
 };
 const whole=nodeFor(cordTerms[0],cordTerms)as Extract<AnatomyNode,{kind:'group'}>;
 const coverings=inventory.terms.filter(t=>t.id.startsWith('A14.1.01.'));
 const roots=inventory.terms.filter(t=>t.id.startsWith('A14.2.00.'));
 for(const [name,rows]of [['Coverings',coverings],['Nerve roots',roots]]as const){
  const nested=rows.filter(row=>!rows.some(parent=>parent.id===row.parent)).map(row=>nodeFor(row,rows));
  if(nested.length)whole.nodes.push({kind:'group',id:'spinal:'+name,name,nodes:nested,parts:unique(nested.flatMap(n=>n.parts)),terminology:terminologyForConcept('')});
 }
 whole.parts=unique([...whole.parts,...(cord.marker.coveringPartIds??[]).map(id=>cord.parts.find(p=>p.id===id)!)]);
 whole.itemCount=cordTerms.length;
 const visit=(node:AnatomyNode):AnatomyNode=>{
  if(node.kind==='group'&&node.name==='Spinal cord'&&cord.marker.partIds.every(id=>node.parts.some(p=>p.id===id)))return whole;
  return node.kind==='group'?{...node,nodes:node.nodes.map(visit)}:node;
 };
 return nodes.map(visit);
}
