import type {StencilCapCandidate} from './stencil-caps';

type ClosureMatrix=Pick<StencilCapCandidate['matrix'],'elements'>;
export type SectionClosureCandidate=Pick<StencilCapCandidate,'part'|'cutIndex'>&{matrix?:ClosureMatrix};

/**
 * Eligibility only: a declaration does not prove a surface or solid is closed.
 * Apply to all visible candidates before cut-range culling, assembly grouping,
 * or GPU/CPU jobs. Required peers must repeat the same complete declaration.
 */
export function filterSectionClosureCandidates<T extends SectionClosureCandidate>(
 visible:readonly T[],matrixFor:(candidate:T)=>ClosureMatrix|undefined=candidate=>candidate.matrix,
):T[]{
 if(!visible.some(candidate=>candidate.part.sectionClosureParts!==undefined))return [...visible];
 const declaredAssemblies=new Set(visible.filter(candidate=>candidate.part.sectionClosureParts!==undefined&&candidate.part.sectionAssembly)
  .map(candidate=>JSON.stringify([candidate.cutIndex,candidate.part.sectionAssembly])));
 const byCut=new Map<number,Map<string,T[]>>();
 for(const candidate of visible){
  let parts=byCut.get(candidate.cutIndex);if(!parts){parts=new Map();byCut.set(candidate.cutIndex,parts);}
  const peers=parts.get(candidate.part.id)??[];peers.push(candidate);parts.set(candidate.part.id,peers);
 }
 const declarations=new Map<T,{ids:string[];key:string}|null>();
 const declaration=(candidate:T)=>{
  if(declarations.has(candidate))return declarations.get(candidate)!;
  const ids=candidate.part.sectionClosureParts;
  const valid=Array.isArray(ids)&&ids.length>0&&ids.every(id=>typeof id==='string'&&id.trim().length>0)
   &&new Set(ids).size===ids.length&&ids.includes(candidate.part.id);
  const value=valid?{ids:[...ids].sort(),key:JSON.stringify([...ids].sort())}:null;
  declarations.set(candidate,value);return value;
 };
 const matrices=new Map<T,readonly number[]|null>();
 const matrix=(candidate:T)=>{
  if(matrices.has(candidate))return matrices.get(candidate)!;
  const elements=matrixFor(candidate)?.elements;
  const value=elements?.length===16&&elements.every(Number.isFinite)?Array.from(elements):null;
  matrices.set(candidate,value);return value;
 };
 return visible.filter(candidate=>{
  if(candidate.part.sectionClosureParts===undefined)return !declaredAssemblies.has(JSON.stringify([candidate.cutIndex,candidate.part.sectionAssembly]));
  const required=declaration(candidate),assembly=candidate.part.sectionAssembly;
  if(!required||!assembly?.trim())return false;
  const reference=matrix(candidate);if(!reference)return false;
  return required.ids.every(id=>{
   const peers=byCut.get(candidate.cutIndex)?.get(id);
   if(peers?.length!==1)return false;
   const peer=peers[0];
   if(peer.part.sectionAssembly!==assembly||declaration(peer)?.key!==required.key)return false;
   const transform=matrix(peer);
   return !!transform&&reference.every((value,index)=>value===transform[index]);
  });
 });
}
