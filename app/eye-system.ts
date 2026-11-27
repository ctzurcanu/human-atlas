import type {Atlas,Part} from './anatomy';
import type {AnatomyEntry,AnatomyNode} from './anatomy-hierarchy';
import {terminologyForConcept,terminologyForGroup} from './anatomical-terminology';
import {eyePresentationParts} from './eye-membership';

type EyeGroup={side:string;eyeConceptId:string;canonicalEyeConceptId?:string;activeComponentIds:string[];ocularVesselConceptIds?:string[];ocularVesselComponentIds?:string[]};
type EyeAtlas=Atlas&{eyePhysicalComponentMembership?:{wholeEyeConceptId:string;wholeOcularConceptId?:string;groups:EyeGroup[]}};
const unique=(parts:Part[])=>[...new Map(parts.map(part=>[part.id,part])).values()];

/** Eye is an organ assembly as well as a set of blood-system structures. */
export function withEyeSystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const registry=(atlas as EyeAtlas).eyePhysicalComponentMembership;
 if(!registry?.groups.every(group=>group.ocularVesselConceptIds?.length))return nodes;
 const parts=new Map(atlas.parts.filter(part=>!part.suppressed).map(part=>[part.id,part]));
 const concepts=new Map(atlas.concepts.map(concept=>[concept.id,concept]));
 const owned=new Set(registry.groups.flatMap(group=>group.activeComponentIds));
 const entryFor=(id:string,allowed:Set<string>):AnatomyNode=>{
  const concept=concepts.get(id)!;
  const members=unique(concept.elements.filter(pid=>allowed.has(pid)).map(pid=>parts.get(pid)).filter((p):p is Part=>!!p));
  const entry:AnatomyEntry={id,name:concept.name,system:members[0].system,region:'Head & neck',location:[],ancestry:[],parts:members,
   ta98Term:concept.ta98Term,ta98Kind:concept.ta98Kind,terminology:terminologyForConcept(id,concept.ta98Term,concept.ta98Kind),
   ...(members.every(part=>part.system==='arterial'||part.system==='venous')?{selectionOnly:true as const}:{})};
  return {kind:'entry',id:'eye:component:'+id,name:concept.name,parts:members,entry};
 };
 const eyes:AnatomyNode[]=registry.groups.map(group=>{
  const eyeId=group.canonicalEyeConceptId??group.eyeConceptId,concept=concepts.get(eyeId)!;
  const active=new Set(group.activeComponentIds);
  const blood=new Set(group.ocularVesselComponentIds);
  const nativeIds=[...new Set(group.activeComponentIds.filter(pid=>!blood.has(pid)).map(pid=>parts.get(pid)!.conceptId))];
  const circulation:AnatomyNode={kind:'group',id:'eye:circulation:'+eyeId,name:'Blood vessels',
   nodes:group.ocularVesselConceptIds!.map(id=>entryFor(id,active)),parts:group.ocularVesselComponentIds!.map(pid=>parts.get(pid)!),terminology:terminologyForGroup('Blood vessels')};
  return {kind:'group',id:'eye:whole:'+eyeId,name:concept.name,nodes:[...nativeIds.map(id=>entryFor(id,active)),circulation],
   parts:group.activeComponentIds.map(pid=>parts.get(pid)!),terminology:terminologyForConcept(concept.id,concept.ta98Term,concept.ta98Kind)};
 });
 const hidden=new Set(atlas.defaultHidden??[]);
 const prune=(node:AnatomyNode,insideEye=false):AnatomyNode[]=>{
  insideEye=insideEye||node.name==='Eye';
  const keepPart=(part:Part)=>!owned.has(part.id)&&(!insideEye||!hidden.has(part.id));
  if(node.kind==='entry'){const keep=node.parts.filter(keepPart);return keep.length?[{...node,parts:keep,entry:{...node.entry,parts:keep}}]:[];}
  if(node.kind==='bilateral'){const entries=node.entries.flatMap(entry=>{const keep=entry.parts.filter(keepPart);return keep.length?[{...entry,parts:keep}]:[];});return entries.length?[{...node,entries,parts:unique(entries.flatMap(entry=>entry.parts))}]:[];}
  const children=node.nodes.flatMap(child=>prune(child,insideEye));return children.length?[{...node,nodes:children,parts:unique(children.flatMap(child=>child.parts))}]:[];
 };
 const ocular:AnatomyNode={kind:'group',id:'eye:both',name:'Eyes',nodes:eyes,parts:unique(eyes.flatMap(eye=>eye.parts)),
  terminology:registry.wholeOcularConceptId?terminologyForGroup('Eyes'):terminologyForConcept(registry.wholeEyeConceptId,'A15.2.01.001','parent')};
 const rest=nodes.flatMap(node=>prune(node));
 const extendLegacy=(rows:AnatomyNode[]):boolean=>{
  for(const node of rows)if(node.kind==='group'){
   if(node.name==='Eye'){
    node.nodes.push(ocular);node.parts=eyePresentationParts(atlas,unique([...node.parts,...ocular.parts]),{id:node.id,name:node.name,term:node.terminology.ta98??undefined});return true;
   }
   if(extendLegacy(node.nodes)){node.parts=unique(node.nodes.flatMap(child=>child.parts));return true;}
  }
  return false;
 };
 const insert=(rows:AnatomyNode[]):boolean=>{
  for(const node of rows)if(node.kind==='group'){
   if(node.name==='Sensory organs'){node.nodes.push(ocular);node.parts=unique([...node.parts,...ocular.parts]);return true;}
   if(insert(node.nodes)){node.parts=unique(node.nodes.flatMap(child=>child.parts));return true;}
  }
  return false;
 };
 if(!extendLegacy(rest)&&!insert(rest))rest.push({kind:'group',id:'eye:sensory',name:'Sensory organs',nodes:[ocular],parts:ocular.parts,terminology:terminologyForGroup('Sensory organs')});
 return rest;
}
