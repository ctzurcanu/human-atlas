import {catalogueSourceReference} from "./source-reference";
import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';

const prefix='BDM:male:human_20_male_cardiovascular_system-';
/** Only the reviewed source collection and its five splenic recipient sites.
 * This identifies drainage; it does not certify finite walls or patent joins. */
export function reviewedPancreaticVenousPair(parts:readonly Part[]):[Part,Part]|undefined{
 const source=parts.find(p=>p.id===prefix+'pancreatic_veins_ID'&&p.ta98Term==='A12.3.12.029');
 const recipient=parts.find(p=>p.id===prefix+'splenic_vein_ID'&&p.ta98Term==='A12.3.12.028');
 if(!source||!recipient)return;
 if(![source,recipient].every(p=>p.system==='venous'&&!p.suppressed&&!p.defaultHidden
  &&p.sourceId&&p.id==='BDM:male:'+p.sourceId&&p.sourceGroups?.includes(prefix+'root')
  &&catalogueSourceReference(p.provenance?.url)==='/ATTRIBUTION.md#registered-collection-source'))return;
 return [source,recipient];
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function pancreaticVenousRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),pair=reviewedPancreaticVenousPair(parts);
 if(pair){const [source,recipient]=pair,note='Named pancreatic tributaries drain into the splenic vein at five source recipient sites. Finite vessel walls and patent lumen junctions remain unverified.';
  graph.set(source.id,[{kind:'after',target:recipient,note}]);
  graph.set(recipient.id,[{kind:'before',target:source,note}]);
 }
 cache.set(parts,graph);return graph;
}
