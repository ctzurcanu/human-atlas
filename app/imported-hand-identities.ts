import type {Atlas} from './anatomy';

const handVeins=new Set(['left','right'].map(side=>`BDM:male:human_20_male_cardiovascular_system-${side}_dorsal_digital_veins_ID`));
const foot='A12.3.11.014',handParent='A12.3.08.023';

/** The authority's digitales dorsales pedis term cannot identify hand vessels. */
export function withImportedHandVeinIdentities(atlas:Atlas):Atlas{
 let changed=false;
 const concepts=atlas.concepts.map(c=>{
  if(!handVeins.has(c.id)||c.ta98Term!==foot)return c;
  changed=true;return {...c,ta98Term:handParent,ta98Kind:'parent' as const};
 });
 const parts=atlas.parts.map(p=>{
  if(!handVeins.has(p.conceptId)||p.ta98Term!==foot)return p;
  changed=true;return {...p,ta98Term:undefined};
 });
 const ta98TermBindings=atlas.ta98TermBindings?.map(b=>{
  if(!handVeins.has(b.conceptId)||b.term!==foot)return b;
  changed=true;return {...b,term:handParent,kind:'parent' as const};
 });
 return changed?{...atlas,concepts,parts,ta98TermBindings}:atlas;
}
