import {SYSTEMS,type Part} from './anatomy';

export interface AnatomyPathNode {label:string;search:string}

const REGION=/^(head|neck|trunk|torso|thorax|abdomen|pelvis|back|left (upper|lower) limb|right (upper|lower) limb)$/i;
const NERVOUS_ORDER:Record<string,number>={
 'central nervous system':1,'peripheral nervous system':1,
 'brain':2,'spinal cord':2,
 'cerebrum':3,'cerebellum':3,'brainstem':3,
 'telencephalon':4,'diencephalon':4,
 'left cerebral hemisphere':5,'right cerebral hemisphere':5,'right cerebal hemisphere':5,
 'neo-cortex':6,
};
const pathKey=(label:string)=>label.replace(/^\s*(?:left|right)\s+/i,'').replace(/\s*(?:\((?:left|right)\)|,\s*(?:left|right))\s*$/i,'').replace(/\s*·.*$/u,'').trim().toLocaleLowerCase();

export function anatomyPath(name:string,selectedParts:Part[],allParts:Iterable<Part>,sex:'male'|'female',scope?:string):AnatomyPathNode[]{
 const root=scope==='cell'?'Human Cell':scope==='embryo'?'Embryo':sex==='male'?'Human Male':'Human Female';
 const path:AnatomyPathNode[]=[{label:root,search:scope==='cell'?'Cell_(biology)':scope==='embryo'?'Embryo':sex==='male'?'Man':'Woman'}];
 const systemIds=[...new Set(selectedParts.map(part=>part.system))];
 const system=systemIds.length===1?SYSTEMS.find(item=>item.id===systemIds[0]):undefined;
 if(system)path.push({label:system.name.replace(/\b\w/g,letter=>letter.toUpperCase()),search:system.name});

 const first=selectedParts[0];
 const shared=first?.groups?.filter(group=>selectedParts.every(part=>part.groups?.includes(group)))??[];
 const groups=[...new Set(shared)].filter(group=>!/^\d+:\s/.test(group)&&group!=='Bonus collection');
 const counts=new Map<string,number>();
 for(const part of allParts)for(const group of new Set(part.groups??[]))counts.set(group,(counts.get(group)??0)+1);
 const tier=(group:string)=>REGION.test(group)?0:system?.id==='nervous'?(NERVOUS_ORDER[group.toLowerCase()]??7):1;
 groups.sort((a,b)=>tier(a)-tier(b)||(counts.get(b)??0)-(counts.get(a)??0)||a.localeCompare(b));
 const seen=new Set([pathKey(root),pathKey(name),...(system?[pathKey(system.name)]:[])]);
 for(const group of groups){
  const label=group.replace(/cerebal/i,'cerebral');
  const key=pathKey(label);if(!key||seen.has(key))continue;seen.add(key);
  path.push({label,search:label});
 }
 return path;
}
