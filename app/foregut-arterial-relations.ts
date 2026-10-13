import type {Part} from './anatomy';
import type {RelationKind,ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const names=new Set(['Coeliac trunk','Celiac trunk','Left gastric artery','Gastric artery (left)','Common hepatic artery','Proper hepatic artery','Gastroduodenal artery','Splenic artery']);
const nativeLeaf=(p:Part)=>!!p.sourceId&&(p.id===`ZA:${p.sourceId}`||!!nativePartFrame(p)?.startsWith('LOCAL:')&&p.id.startsWith(`${nativePartFrame(p)}${p.sourceId}:`));
/** A native arterial itinerary cannot jump into an independently placed donor.
 * This guard also applies to reciprocal generic routes in the final resolver. */
export function foregutArterialFrameCompatible(kind:RelationKind,a:Part,b:Part):boolean{
 if(!['before','after','adjacent'].includes(kind)||![a,b].some(p=>names.has(p.name)))return true;
 const frame=nativePartFrame(a),other=nativePartFrame(b);
 if(!frame&&!other)return true;
 return !a.suppressed&&!b.suppressed&&a.system==='arterial'&&b.system==='arterial'&&!!frame&&frame===other&&[a,b].every(p=>!names.has(p.name)||nativeLeaf(p));
}
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** NBK459241: usual foregut arterial branching, with anatomical variants.
 * Named sources and shared rims do not prove complete walls or lumen joins. */
export function foregutArterialRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),active=parts.filter(p=>!p.suppressed&&p.system==='arterial'&&nativePartFrame(p)&&nativeLeaf(p)&&names.has(p.name));
 for(const [parents,children] of [
  [['Coeliac trunk','Celiac trunk'],['Left gastric artery','Gastric artery (left)','Common hepatic artery','Splenic artery']],
  [['Common hepatic artery'],['Proper hepatic artery','Gastroduodenal artery']],
 ])for(const a of active.filter(p=>parents.includes(p.name)))for(const b of active.filter(p=>children.includes(p.name)&&nativePartFrame(p)===nativePartFrame(a))){
  const note='Usual named foregut arterial branching in this native model. Arterial variants occur; complete vessel walls, biological lumen continuity and source attachment anatomy remain unverified.';
  for(const [source,target,kind] of [[a,b,'after'],[b,a,'before']] as const){const rows=graph.get(source.id)??[];rows.push({kind,target,note});graph.set(source.id,rows);}
 }
 cache.set(parts,graph);return graph;
}
