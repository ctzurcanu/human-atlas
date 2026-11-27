import type {Atlas,Concept,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {hierarchyEntries,majorSystemFor,systemPathFor} from './anatomy-hierarchy';

type EyeGroup={side:string;eyeConceptId:string;canonicalEyeConceptId?:string;activeComponentIds:string[];ocularVesselComponentIds?:string[]};
type EyeAtlas=Atlas&{eyePhysicalComponentMembership?:{groups:EyeGroup[]}};
type Membership={native:Set<string>;globe:Set<string>;blood:Part[];active:Set<string>;activeParts:Part[]};
const cache=new WeakMap<Atlas,Membership|null>();
const unique=(parts:Part[])=>[...new Map(parts.map(part=>[part.id,part])).values()];
const sameSet=(ids:string[],set:Set<string>)=>ids.length===set.size&&ids.every(id=>set.has(id));

/** The existing physical registry owns the ocular vessels. TA98 ancestry only
 * identifies the already-present native Eye/adnexa presentation around it. */
function membership(atlas:Atlas):Membership|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const groups=(atlas as EyeAtlas).eyePhysicalComponentMembership?.groups;
 const parts=new Map(atlas.parts.filter(part=>!part.suppressed).map(part=>[part.id,part]));
 if(!groups?.length||groups.some(group=>!group.ocularVesselComponentIds?.length||group.activeComponentIds.some(id=>!parts.has(id))||group.ocularVesselComponentIds.some(id=>!group.activeComponentIds.includes(id)||!parts.has(id)))){cache.set(atlas,null);return null;}
 const hidden=new Set(atlas.defaultHidden??[]),native=new Set<string>(),globe=new Set<string>();
 for(const entry of hierarchyEntries(atlas)){
  if(majorSystemFor(entry)!=='nervous')continue;
  const path=systemPathFor(entry),eye=path.includes('Eye'),ball=path.includes('Eyeball');
  if(!eye)continue;
  for(const part of entry.parts){if(hidden.has(part.id)||part.referenceOnly||['arterial','venous'].includes(part.system))continue;native.add(part.id);if(ball)globe.add(part.id);}
 }
 const blood=unique(groups.flatMap(group=>group.ocularVesselComponentIds!.map(id=>parts.get(id)!))),active=new Set(groups.flatMap(group=>group.activeComponentIds));
 const result={native,globe,blood,active,activeParts:[...active].map(id=>parts.get(id)!)};cache.set(atlas,result);return result;
}

/** Extend a legacy whole-Eye node, retaining native adnexa and source regions.
 * This membership never reclassifies an extraocular vessel as intraocular. */
export function eyePresentationParts(atlas:Atlas,existing:Part[],identity:{id?:string;name?:string;term?:string}):Part[]{
 const registry=membership(atlas);if(!registry)return existing;
 const term=identity.term??identity.id?.replace(/^(?:TA98GROUP:|TA98:|hierarchy:guest:ta98:)/,'');
 const whole=term==='A15.2.00.001'||identity.name?.toLowerCase()==='eye';
 const ball=term==='A15.2.01.001'||identity.name?.toLowerCase()==='eyeball';
 if(!whole&&!ball||!existing.some(part=>registry.native.has(part.id)||registry.active.has(part.id)))return existing;
 const hidden=new Set(atlas.defaultHidden??[]);
 const native=existing.filter(part=>!hidden.has(part.id)&&!part.suppressed&&!part.referenceOnly);
 const sides=new Set(native.flatMap(part=>/\.l$|\bleft\b/i.test(part.id+' '+part.name)?['left']:/\.r$|\bright\b/i.test(part.id+' '+part.name)?['right']:[]));
 const modeled=registry.activeParts.filter(part=>sides.size!==1||(/\bleft\b|left_/i.test(part.name+' '+part.id)?sides.has('left'):/\bright\b|right_/i.test(part.name+' '+part.id)?sides.has('right'):false));
 return unique([...native,...modeled]);
}

/** Old packed whole-Eye selections contain the exact native set but no choice
 * ID. Upgrade only those complete known sets; partial/custom selections stay. */
export function restoreEyeAssemblySelection(atlas:Atlas,selected:string[]):string[]{
 const registry=membership(atlas);if(!registry)return selected;
 const ids=[...new Set(selected)];
 if(!sameSet(ids,registry.native)&&!sameSet(ids,registry.globe))return selected;
 return [...new Set([...selected,...registry.blood.map(part=>part.id)])];
}

/** A literal generic Eye query includes the modeled native accessory organs.
 * The named FMA Eye/eyeball aliases retain their qualified 18-part subsets. */
export function eyeSelectionForTerm(atlas:Atlas,term:string):string[]|undefined{
 const registry=membership(atlas);if(!registry)return;
 const code=term.replace(/^(?:TA98GROUP:|TA98:|hierarchy:guest:ta98:)/,'');
 if(code==='A15.2.00.001'||term.toLowerCase()==='eye')return [...registry.native,...registry.blood.map(part=>part.id)];
}

/** A restored legacy whole-Eye view needs its organ label, not "58 structures".
 * This is a presentation assembly, not an additional TA98 identity or mesh. */
export function restoredEyeAssemblyChoice(atlas:Atlas,selected:string[]):Concept|undefined{
 const registry=membership(atlas);if(!registry)return;
 const expected=new Set([...registry.native,...registry.blood.map(part=>part.id)]);
 if(!sameSet([...new Set(selected)],expected))return;
 return {id:'LOCAL:modeled-eye-and-adnexa',name:'Eye',elements:[...selected],componentElements:[...selected],coverageLimitation:'Currently modeled native Eye and accessory structures, including the explicitly registered ocular circulation. Intraocular beds and extraocular vessel courses retain their separate source roles. This selection does not establish complete eye anatomy, accepted vessel placement, finite walls or lumen patency.'};
}

/** Complete named Eye choices already carry their physical member IDs. The
 * registered side concepts provide independent ownership for their vessels. */
export function eyeChoiceComponentIds(choice:HierarchyChoice,parts:Map<string,Part>,concepts:Concept[]):string[]|undefined{
 const term=choice.ta98Term??choice.terminology?.ta98;
 if(choice.id!=='LOCAL:modeled-eye-and-adnexa'&&term!=='A15.2.00.001'&&term!=='A15.2.01.001'&&choice.name.toLowerCase()!=='eye'&&choice.name.toLowerCase()!=='eyeball')return;
 const declared=new Set(concepts.filter(concept=>['FMA54449','FMA54450'].includes(concept.id)).flatMap(concept=>concept.componentElements??concept.elements));
 if(!choice.elements.some(id=>declared.has(id)&&['arterial','venous'].includes(parts.get(id)?.system??'')))return;
 return [...new Set(choice.elements)].filter(id=>parts.has(id)&&!parts.get(id)!.suppressed);
}
