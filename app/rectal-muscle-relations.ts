import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';
import {largeBowelWallSource} from './large-bowel-frame';

const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const name=(part:Part)=>part.name.toLowerCase().replace(/^(?:left|right)\s+/,'').replace(/\s*\((?:left|right)\)$/,'');

/** The puborectal sling surrounds the anorectal junction. This topographic
 * relation is not part-of, a validated tissue attachment, or a movement group.
 * Source: Wu et al., PLOS ONE 2015, doi:10.1371/journal.pone.0132226.
 * Resolve only active native peers in the same established source frame.
 */
export function rectalMuscleRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const active=parts.filter(part=>!part.suppressed&&nativePartFrame(part));
 const rectums=active.filter(part=>part.system==='digestive'&&(
  ['rectum','large intestine rectum'].includes(name(part))||
  part.caecocolicRegion==='rectum'&&part.caecocolicLayer==='serosa'&&!!largeBowelWallSource(part,parts)
 ));
 const muscles=active.filter(part=>part.system==='muscular'&&name(part)==='puborectalis');
 const note='Puborectalis forms a sling around the anorectal junction and contributes to the anorectal bend.';
 const add=(source:Part,target:Part)=>graph.set(source.id,[...(graph.get(source.id)??[]),{kind:'adjacent',target,note}]);
 for(const rectum of rectums)for(const muscle of muscles){
  if(nativePartFrame(rectum)!==nativePartFrame(muscle))continue;
  add(rectum,muscle);add(muscle,rectum);
 }
 cache.set(parts,graph);return graph;
}
