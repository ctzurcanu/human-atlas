import {isSurfaceOnly,surfaceWholes} from './surface-companion';
import {liverCompletionLabelReference} from './liver-completion-controls';
import type {Atlas,Part,Concept} from './anatomy';
import {structureLabelName} from './laterality';
import {anatomicalComponentIds} from './anatomical-components';
import {anatomicalRelations} from './anatomical-relations';
import type {HierarchyChoice} from './hierarchy-choice';
import {catalogueTa98TermParts} from './ta98-term-bindings';
import {urinaryLabelReference} from './urinary-wall';
import {diaphragmLabelReference} from './diaphragm-controls';
import {liverSurfaceLabelReference} from './liver-surface-features';
import {carotidSurfaceLabelReference} from './carotid-surface-features';
import {cardiacSurfaceLabelReference} from './cardiac-surface-features';
import {cardiacSiteLabelReference,cardiacComponentChoices,cardiacLabelRegistry} from './cardiac-label-sites';
import {liverLandmarkLabelReference} from './liver-landmark';
import {spinalCordComponentChoices} from './spinal-cord-presentation';
import {tongueComponentChoices,tongueLabelReference} from './tongue-controls';

export interface AnatomyLabelItem {id:string;name:string;members:string[];components:string[];choice?:HierarchyChoice;anchor?:LabelSurfacePoint}
export interface LayerLabelTarget {axis:0|1|2;min:number;max:number;coordinate:number}
/** A saved presentation preference must not suppress the current selection's
 * identifying callout. General exploded-group labels remain optional. */
export {selectedLabelsVisible} from './label-visibility';
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
 // A reviewed singular course names the existing carrier without inventing component anatomy.
 const course=choice?atlas.ta98TermBindings?.find(b=>b.representation==='source-course-identity'&&(choice.id===b.conceptId||choice.id===`hierarchy:guest:ta98:${b.term}`)):undefined;
 const courseParts=course?catalogueTa98TermParts({...atlas,ta98TermBindings:[course]}).get(course.term):undefined;
 if(choice&&courseParts?.length===1&&choice.elements.length===1&&choice.elements[0]===courseParts[0].id&&ids.has(courseParts[0].id)){
  result.push({id:choice.id,name:choice.name,members:[courseParts[0].id],components:[],choice});used.add(courseParts[0].id);
 }
 const cardiacComponents=expanded&&choice?cardiacComponentChoices(atlas,choice):undefined;
 if(cardiacComponents){
  // Validate one registry for the whole expansion. Revalidating all reference
  // geometry for each of 90 callouts needlessly stalls the main thread.
  const sites=new Map(cardiacLabelRegistry(atlas)?.termGroups.map(g=>[g.conceptId,g]));
  return cardiacComponents.filter(c=>c.elements.every(id=>ids.has(id))).map(c=>{
   const g=sites.get(c.id),site=g&&g.partIds.length===c.elements.length&&g.partIds.every(id=>c.elements.includes(id))?g:undefined;
   return {id:c.id,name:c.name,members:[...c.elements],components:[],choice:c,anchor:cardiacSurfaceLabelReference(atlas,c)??(site?{partId:site.anchor.partId,point:site.anchor.point}:undefined)};
  });
 }
 const tongueComponents=expanded&&choice?tongueComponentChoices(atlas,choice):undefined;
 if(tongueComponents)return tongueComponents.filter(c=>c.elements.every(id=>ids.has(id))).map(c=>({id:c.id,name:c.name,members:[...c.elements],components:[],choice:c,anchor:tongueLabelReference(atlas,c)}));
 const spinalComponents=expanded&&choice?spinalCordComponentChoices(atlas,choice):undefined;
 if(spinalComponents)return spinalComponents.length?spinalComponents.filter(c=>c.elements.every(id=>ids.has(id))).map(c=>({id:c.id,name:c.name,members:[...c.elements],components:[],choice:c})):[{id:choice!.id,name:choice!.name,members:[...choice!.elements],components:[],choice}];
 const reference=choice&&(carotidSurfaceLabelReference(atlas,choice)??tongueLabelReference(atlas,choice)??urinaryLabelReference(atlas,choice)??diaphragmLabelReference(atlas,choice)??cardiacSurfaceLabelReference(atlas,choice)??cardiacSiteLabelReference(atlas,choice)??liverSurfaceLabelReference(atlas,choice)??liverLandmarkLabelReference(atlas,choice)??liverCompletionLabelReference(atlas,choice));
 if(reference&&choice&&choice.elements.every(id=>ids.has(id))){
  result.push({id:choice.id,name:choice.name,members:[...choice.elements],components:[],choice,anchor:reference});
  choice.elements.forEach(id=>used.add(id));
 }
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
 for(const id of ids){if(used.has(id))continue;const part=parts.get(id)!;if(part.conceptId.startsWith('TA98COURSE:'))continue;
  // A surface (no volume) is labelled as itself; its smallest whole is labelled on its own geometry below.
  if(isSurfaceOnly(part)){result.push({id:part.id,name:structureLabelName(part.name),members:[id],components:[]});continue;}
  // Several selected declared subdivisions of one whole are separate named items, not one merged label.
  const decl=concepts.get(part.conceptId)?.componentElements;
  if(decl?.includes(id)&&decl.filter(c=>ids.has(c)).length>1){result.push({id:part.id,name:structureLabelName(part.name),members:[id],components:[]});continue;}
  const whole0=concepts.get(part.conceptId),decl0=whole0?.componentElements;
  if(decl0&&!decl0.includes(id)&&decl0.filter(c=>ids.has(c)).length>1)continue;   // the whole itself is shown, labelled through its subdivisions
  let item=byConcept.get(part.conceptId);
  if(!item){const whole=concepts.get(part.conceptId),declared=whole?.componentElements;
   // A selected declared subdivision keeps its whole's label and dropdown.
   item={id:part.conceptId,name:whole?.name??structureLabelName(part.name),members:[],components:declared?.includes(id)&&declared.every(c=>parts.has(c))?[...declared]:[],...(declared?.includes(id)&&whole?{choice:{...whole,elements:[...whole.elements]}}:{})};byConcept.set(part.conceptId,item);result.push(item);}item.members.push(id);
 }
 // The whole unit drawn transparent around a selected surface carries the whole's own label and dropdown, anchored on its geometry.
 for(const {concept,ids:wholeIds} of surfaceWholes(atlas,[...ids])){
  const declared=concept.componentElements;
  result.push({id:concept.id,name:structureLabelName(concept.name),members:wholeIds,components:declared&&declared.every(c=>parts.has(c))?[...declared]:[],choice:{...concept,elements:[...concept.elements]}});
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
 return `bank:${part.chunk}:${part.positions}:${part.geometryOrigin?.join(',')??''}:${part.geometryScale??1}:${JSON.stringify(part.geometryInstance??null)}`;
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
 for(const item of items){if(item.anchor)continue;const keys=new Set(item.members.map(id=>byId.get(id)).filter((p):p is Part=>!!p).map(familyFor));
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
