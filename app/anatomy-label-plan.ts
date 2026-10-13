import {structureName,type Atlas,type Part,type Concept} from './anatomy';
import {anatomicalComponentIds} from './anatomical-components';
import {anatomicalRelations} from './anatomical-relations';
import type {HierarchyChoice} from './hierarchy-choice';

export interface AnatomyLabelItem {id:string;name:string;members:string[];components:string[];choice?:HierarchyChoice}
export interface LayerLabelTarget {axis:0|1|2;min:number;max:number;coordinate:number}
/** Select the item actually named by this label, not its surrounding set or
 * other regions absent from the callout. This also bypasses Add selection. */
export function anatomyLabelChoice(item:AnatomyLabelItem):HierarchyChoice{
 return {...item.choice,id:item.id,name:item.name,elements:[...item.members]};
}
const catalogueParts=new WeakMap<Atlas,Map<string,Part>>();
const catalogueComponents=new WeakMap<Atlas,Map<string,string[]>>();

/** Assemblies come from catalogue membership, never from overlapping bounds. */
export function selectedAnatomyLabels(atlas:Atlas,selected:string[],choice:HierarchyChoice|undefined,expanded:boolean):AnatomyLabelItem[]{
 let parts=catalogueParts.get(atlas);if(!parts){parts=new Map(atlas.parts.map(p=>[p.id,p]));catalogueParts.set(atlas,parts);}
 const ids=new Set(selected.filter(id=>parts.has(id)&&!parts.get(id)!.suppressed));
 const components=(group:HierarchyChoice,resolveRelations=false)=>{
  const stored=anatomicalComponentIds(group,[],parts,atlas.concepts);
  if(stored.length||!resolveRelations)return stored;
  let cache=catalogueComponents.get(atlas);if(!cache){cache=new Map();catalogueComponents.set(atlas,cache);}
  const key=`${group.id}:${group.elements.join('|')}`,prior=cache.get(key);if(prior)return prior;
  const resolved=anatomicalComponentIds(group,group.elements.flatMap(id=>{const p=parts.get(id);return p?anatomicalRelations(p,atlas.parts):[];}),parts,atlas.concepts);cache.set(key,resolved);return resolved;
 };
 const result:AnatomyLabelItem[]=[],used=new Set<string>();
 const addGroup=(group:HierarchyChoice,children:string[])=>{
  const members=group.elements.filter(id=>ids.has(id)&&!used.has(id));
  if(!members.length)return;
  result.push({id:group.id,name:group.name,members,components:children,choice:group});members.forEach(id=>used.add(id));
 };
 if(!expanded&&choice&&choice.id!=='selection-set'&&choice.elements.every(id=>ids.has(id))){
  const children=components(choice,true);
  if(children.length)addGroup(choice,children);
 }
 if(!expanded){
  // Largest fully selected stored assembly first. Partial selections must not
  // silently acquire their unselected parent or another organ's tissues.
  const candidates=atlas.concepts.filter(c=>c.elements.length>1&&c.elements.every(id=>ids.has(id)&&!used.has(id))).sort((a,b)=>b.elements.length-a.elements.length||a.id.localeCompare(b.id));
  for(const group of candidates){if(group.elements.some(id=>used.has(id)))continue;const children=components(group);if(children.length)addGroup(group,children);}
 }
 const byConcept=new Map<string,AnatomyLabelItem>(),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 for(const id of ids){if(used.has(id))continue;const part=parts.get(id)!;let item=byConcept.get(part.conceptId);
  if(!item){item={id:part.conceptId,name:concepts.get(part.conceptId)?.name??structureName(part.name),members:[],components:[]};byConcept.set(part.conceptId,item);result.push(item);}item.members.push(id);
 }
 return result;
}

/** Source/assembly identities distinguish neighbouring organs sharing a bank. */
export function layerLabelFamily(part:Part):string|undefined{
 if(part.suppressed)return;
 if(part.tissueWallAssembly?.startsWith('gallbladder:')){
  const region=part.gallbladderRegion??/(A05\.8\.02\.00[2-5])$/.exec(part.id)?.[1];
  if(region)return `${part.tissueWallAssembly}:${region}`;
 }
 if(part.caecocolicWallAssembly)return `${part.caecocolicWallAssembly}:${part.caecocolicRegion??part.bowelWallSource??''}`;
 if(part.tissueWallAssembly)return part.tissueWallAssembly;
 if(part.bowelWallSource)return part.bowelWallSource;
 if(/^GASTRIC:(male|female):/.test(part.id))return part.id.split(':').slice(0,2).join(':');
 if(!part.pairedTissue)return;
 return `bank:${part.chunk}:${part.positions}:${part.geometryOrigin?.join(',')??''}`;
}

/** Space semantic layers through 25–75% of the organ's longest measure. */
export function layerLabelTargets(items:AnatomyLabelItem[],parts:Part[],concepts:Concept[]=[]):Map<string,LayerLabelTarget>{
 const byId=new Map(parts.map(p=>[p.id,p])),byConcept=new Map(concepts.map(c=>[c.id,c])),families=new Map<string,AnatomyLabelItem[]>();
 const familyFor=(part:Part)=>{
  const stored=layerLabelFamily(part);if(stored||part.suppressed)return stored;
  const identity=byConcept.get(part.conceptId),parent=identity?.partOf&&byConcept.get(identity.partOf);
  // Native catalogue layers can predate finite-wall metadata. Their explicit
  // parent membership still supplies the assembly; no proximity inference.
  if(parent&&parent.elements.includes(part.id)&&/\b(?:tunica|tela|stratum|lamina|epithelium|serosa|subserosa|mucosa|submucosa|muscularis|capsula|capsule|coat|layer)\b/i.test(identity!.name))return `parent:${parent.id}`;
 };
 for(const item of items){const keys=new Set(item.members.map(id=>byId.get(id)).filter((p):p is Part=>!!p).map(familyFor));
  if(keys.size!==1||keys.has(undefined))continue;const key=[...keys][0]!;const family=families.get(key)??[];family.push(item);families.set(key,family);
 }
 const result=new Map<string,LayerLabelTarget>();
 for(const [family,members] of families){if(members.length<2)continue;
  const bounds=[ [Infinity,Infinity,Infinity],[-Infinity,-Infinity,-Infinity] ];
  const parent=family.startsWith('parent:')?byConcept.get(family.slice(7)):undefined;
  for(const part of parts){if(part.suppressed||!(parent?.elements.includes(part.id)||familyFor(part)===family))continue;for(let a=0;a<3;a++){bounds[0][a]=Math.min(bounds[0][a],part.bounds[0][a]);bounds[1][a]=Math.max(bounds[1][a],part.bounds[1][a]);}}
  const spans=bounds[0].map((v,a)=>bounds[1][a]-v),axis=spans.indexOf(Math.max(...spans)) as 0|1|2;
  const min=bounds[0][axis]+spans[axis]*.25,max=bounds[0][axis]+spans[axis]*.75;
  members.sort((a,b)=>a.id.localeCompare(b.id));
  members.forEach((item,i)=>result.set(item.id,{axis,min,max,coordinate:min+(max-min)*i/(members.length-1)}));
 }
 return result;
}

export type LabelSurfacePoint={partId:string;point:[number,number,number]};
/** Prefer a retained point within the requested interval. Sectioning can
 * remove that interval; the nearest remaining tissue is the fallback. */
export function chooseLayerLabelPoint(samples:LabelSurfacePoint[],target:LayerLabelTarget|undefined,center:[number,number,number]):LabelSurfacePoint|undefined{
 let best:LabelSurfacePoint|undefined,score=Infinity;
 const central=target&&samples.some(s=>s.point[target.axis]>=target.min&&s.point[target.axis]<=target.max);
 for(const sample of samples){const p=sample.point;if(central&&target&&(p[target.axis]<target.min||p[target.axis]>target.max))continue;
  const distance=p.reduce((sum,v,a)=>sum+(v-center[a])**2,0);
  const cost=target?100*(p[target.axis]-target.coordinate)**2+distance*.01:distance;
  if(cost<score){score=cost;best=sample;}
 }
 return best;
}
