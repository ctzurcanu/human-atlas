import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const key=(p:Part)=>p.name.toLowerCase().trim();
const branches:Record<string,string[]>={
 'ascending aorta':['left coronary artery','right coronary artery'],
 'left coronary artery':['anterior interventricular artery','anterior interventricular branch','circumflex artery of heart','circumflex branch of left coronary artery'],
 'anterior interventricular artery':['septal branches of anterior interventricular artery'],
 'right coronary artery':['right inferolateral branch of right coronary artery','posterior interventricular branch of right coronary artery','marginal branch of right coronary artery'],
 'circumflex artery of heart':['left marginal artery'],
 'circumflex branch of left coronary artery':['left marginal artery'],
};
const names=new Set(Object.values(branches).flat());
export const isCoronaryArtery=(p:Part)=>p.system==='arterial'&&names.has(key(p));
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
/** Conventional named branching. Native geometry alone does not establish
 * coronary dominance, anastomoses, ostia, walls or continuous lumens.
 * OpenStax Anatomy & Physiology 2e 19.1, Coronary Arteries. */
export function coronaryRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),arteries=parts.filter(p=>!p.suppressed&&p.system==='arterial'&&nativePartFrame(p));
 for(const source of arteries){
  const children=branches[key(source)];if(!children)continue;
  for(const target of arteries.filter(p=>nativePartFrame(p)===nativePartFrame(source)&&children.includes(key(p)))){
   const note='Named native coronary branch relationship. Coronary ostia, wall layers, lumen continuity and individual perfusion territories remain unverified. Posterior interventricular origin and coronary dominance are not inferred from an unqualified vessel.';
   for(const [from,to,kind] of [[source,target,'after'],[target,source,'before']] as const){const links=graph.get(from.id)??[];links.push({kind,target:to,note});graph.set(from.id,links);}
  }
 }
 cache.set(parts,graph);return graph;
}
export function coronarySourceLimitation(p:Part):string|undefined{
 if(p.id==='ZA:Right inferolateral branch of right coronary artery')return 'Named right-coronary source branch. The inferolateral source label has not been established as the TA98 right posterolateral branch; placement under the right coronary artery is provisional. Extent, wall, lumen and perfusion territory remain unverified.';
 if(isCoronaryArtery(p)&&nativePartFrame(p))return 'Named native coronary arterial exterior. Separate wall layers, continuous lumens, ostia, coronary dominance and complete perfusion territories remain unverified.';
}
