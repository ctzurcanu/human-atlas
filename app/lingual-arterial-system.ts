import type {Atlas,Part} from './anatomy';
import type {AnatomyNode} from './anatomy-hierarchy';
import {isInstalledLingualArterialPart,registeredLingualFamily} from './lingual-arterial-relations';

/** Reconcile only the two reviewed lingual banks, retaining their trunk cells. */
export function withLingualArterialSystem(atlas:Atlas,nodes:AnatomyNode[]):AnatomyNode[]{
 const families=(['left','right'] as const).map(side=>registeredLingualFamily(atlas.parts,side));
 if(families.some(family=>!family))return nodes;
 const members=families.flatMap(family=>[...family!.values()]);
 if(members.length!==12||members.some(part=>!isInstalledLingualArterialPart(part)))return nodes;
 const owned=new Set(members.map(part=>part.id));
 const concepts=new Map(atlas.concepts.map(concept=>[concept.id,concept]));
 for(const side of ['left','right'] as const){
  const family=families[side==='left'?0:1]!,trunk=family.get('lingual')!,whole=concepts.get(trunk.conceptId);
  if(whole?.ta98Kind!=='exact'||whole.ta98Term!=='A12.2.05.015'||whole.partOf!==`ZA:External carotid artery.${side[0]}`
   ||whole.elements.length!==1||whole.elements[0]!==trunk.id)return nodes;
  for(const part of family.values())if(part!==trunk){
   const concept=concepts.get(part.conceptId);
   if(concept?.ta98Kind!=='exact'||concept.ta98Term!==part.ta98Term||concept.partOf!==whole.id||!concept.elements.includes(part.id))return nodes;
  }
 }
 const unique=(parts:Part[])=>[...new Map(parts.map(part=>[part.id,part])).values()];
 const valid=(node:AnatomyNode)=>node.kind!=='entry'&&node.terminology?.ta98==='A12.2.05.015'&&node.parts.length>0&&node.parts.every(part=>owned.has(part.id));
 const visit=(siblings:AnatomyNode[]):AnatomyNode[]=>{
  const branch=siblings.find(node=>node.kind==='group'&&valid(node)&&node.parts.length===10);
  const trunk=siblings.find(node=>node.kind==='bilateral'&&valid(node)&&node.parts.length===2);
  if(branch?.kind==='group'&&trunk?.kind==='bilateral'){
   const parts=unique([...branch.parts,...trunk.parts]);
   if(parts.length===12&&members.every(part=>parts.includes(part)))return siblings.flatMap(node=>node===trunk?[]:node===branch
    ?[{...branch,id:'lingual-arterial:male:whole',parts,itemCount:5,nodes:[{...trunk,id:'lingual-arterial:male:trunks',name:'Main trunk'},...branch.nodes]}]
    :[node]);
  }
  return siblings.map(node=>node.kind==='group'?{...node,nodes:visit(node.nodes)}:node);
 };
 return visit(nodes);
}
