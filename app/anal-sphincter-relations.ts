import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const name=(p:Part)=>p.name.toLowerCase().replace(/\s*\((?:left|right)\)$/,'');
const side=(p:Part)=>p.name.match(/\((left|right)\)$/i)?.[1]?.toLowerCase();
const inferiorNames=new Set(['inferior anal nerves','inferior anal nerve','inferior rectal nerves','inferior rectal nerve']);
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const wholeExternal=(p:Part)=>!p.suppressed&&p.system==='muscular'&&!!nativePartFrame(p)&&['external anal sphincter','m. sphincter ani externus','sphincter ani externus'].includes(name(p));

/** Parent nerve route for the innervation hierarchy; no spinal level is inferred. */
export function analSphincterInnervationNames(p:Part):string[]{
 return wholeExternal(p)?['Pudendal nerve']:[];
}

/** TA98 A04.5.04.012 receives an inferior anal (rectal) nerve route.
 * PMID 3097851 reports additional perineal/posterior branches: this route
 * does not establish exclusive supply, reconstructed fascicles or attachment.
 * Resolve existing peers only; never borrow independently registered donors.
 */
export function analSphincterRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const available=parts.filter(p=>!p.suppressed&&nativePartFrame(p));
 const nerves=available.filter(p=>p.system==='nervous');
 const add=(p:Part,r:ResolvedRelation)=>graph.set(p.id,[...(graph.get(p.id)??[]),r]);
 for(const muscle of available.filter(wholeExternal)){
  const peers=nerves.filter(p=>nativePartFrame(p)===nativePartFrame(muscle)&&(!side(muscle)||!side(p)||side(p)===side(muscle)));
  const inferior=peers.filter(p=>inferiorNames.has(name(p)));
  const targets=inferior.length?inferior:peers.filter(p=>name(p)==='pudendal nerve');
  for(const nerve of targets){
   const route=inferior.length?{}:{via:['Inferior anal nerves'],viaModeled:false};
   const note='Inferior anal (inferior rectal) nerve route to the external anal sphincter. Additional perineal or direct sacral branches and individual variants are not fully represented. This functional relationship does not certify nerve–muscle attachment or source placement.';
   add(muscle,{kind:'innervation',target:nerve,...route,note});
   add(nerve,{kind:'innervates',target:muscle,...route,note});
  }
 }
 cache.set(parts,graph);return graph;
}
