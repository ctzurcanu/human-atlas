import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import {gallbladderWallPeelSteps} from './gallbladder-wall';
import {vascularWallPeelSteps} from './vascular-wall';
export type WallPeelStep={name:string;ids:string[]};
const coatNames:Record<string,string>={serosa:'Tunica serosa',subserosa:'Tela subserosa','longitudinal-muscle':'Stratum longitudinale','circular-muscle':'Stratum circulare',submucosa:'Tela submucosa','muscularis-mucosae':'Lamina muscularis mucosae','lamina-propria':'Lamina propria mucosae',epithelium:'Epithelium mucosae'};
const layers=['serosa','subserosa','longitudinal-muscle','circular-muscle','submucosa','muscularis-mucosae','lamina-propria','epithelium'];
function coatName(key:string,parts:Part[]):string{
 if(parts.some(p=>p.caecocolicRegion==='rectum')){
  if(key==='serosa')return 'External covering (rectal distribution unresolved)';
  if(key==='subserosa')return 'External connective tissue (rectal allocation unresolved)';
 }
 return coatNames[key];
}
/** Only complete, explicitly stored wall/feature assemblies get tissue peeling. */
export function bowelWallPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const vascular=vascularWallPeelSteps(atlas,choice);if(vascular)return vascular;
 const gallbladder=gallbladderWallPeelSteps(atlas,choice);if(gallbladder)return gallbladder;
 if(!atlas.bowelWallPrototype||choice.id==='selection-set'||!choice.elements.length)return null;
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),members=choice.elements.map(id=>byId.get(id));
 const allowed=new Set(atlas.bowelWallPrototype.partIds);
 if(members.some(p=>!p||p.suppressed||!allowed.has(p.id)||!p.bowelWallSource||!p.pairedTissue))return null;
 const parts=members as Part[],sourceId=parts[0].bowelWallSource!,assemblyId=parts[0].caecocolicWallAssembly;
 const tissueAssembly=parts[0].tissueWallAssembly;
 if(tissueAssembly?.startsWith('oesophageal:')){
  if(parts.some(p=>p.tissueWallAssembly!==tissueAssembly||p.chunk!==parts[0].chunk||p.positions!==parts[0].positions||p.normals!==parts[0].normals||JSON.stringify(p.geometryOrigin)!==JSON.stringify(parts[0].geometryOrigin)))return null;
  const selectedSources=new Set(parts.map(p=>p.bowelWallSource));
  const cells=atlas.parts.filter(p=>allowed.has(p.id)&&!p.suppressed&&p.tissueWallAssembly===tissueAssembly&&selectedSources.has(p.bowelWallSource));
  const selected=new Set(parts.map(p=>p.id)),same=(expected:Part[])=>expected.length===selected.size&&expected.every(p=>selected.has(p.id));
  const order=['external-connective-coat','longitudinal-muscle','circular-muscle','submucosa','muscularis-mucosae','lamina-propria','epithelium'];
  const keys=same(cells)?order:same(cells.filter(p=>['longitudinal-muscle','circular-muscle'].includes(p.tissueWallLayer??'')))?order.slice(1,3):same(cells.filter(p=>['muscularis-mucosae','lamina-propria','epithelium'].includes(p.tissueWallLayer??'')))?order.slice(4):[];
  const steps=keys.map(key=>({name:key==='external-connective-coat'?'External connective coat':coatNames[key],ids:parts.filter(p=>p.tissueWallLayer===key).map(p=>p.id)}));
  return steps.length>0&&steps.every(step=>step.ids.length)&&steps.flatMap(step=>step.ids).length===parts.length?steps:null;
 }
 const group=choice.id.replace(/^hierarchy:guest:ta98:/,'').match(/^ALIMENTARY:(small-intestine|large-intestine):wall-assembly$/);
 if(group){
  const bank=group[1]==='small-intestine'?':BOWELWALL:':':CAECOWALL:';
  const expected=atlas.parts.filter(p=>allowed.has(p.id)&&!p.suppressed&&p.id.includes(bank));
  const selected=new Set(choice.elements);
  if(expected.length!==selected.size||!expected.every(p=>selected.has(p.id)))return null;
  const steps=layers.map(key=>({name:coatName(key,parts),ids:parts.filter(p=>p.caecocolicLayer===key||p.id.endsWith(':'+key)||key==='epithelium'&&p.id.endsWith(':circular-folds')).map(p=>p.id)}));
  return steps.every(s=>s.ids.length)&&steps.flatMap(s=>s.ids).length===parts.length?steps:null;
 }
 const sharedAssembly=!!assemblyId&&parts.every(p=>p.caecocolicWallAssembly===assemblyId&&p.chunk===parts[0].chunk&&p.positions===parts[0].positions&&p.normals===parts[0].normals&&JSON.stringify(p.geometryOrigin)===JSON.stringify(parts[0].geometryOrigin));
 if(!sharedAssembly&&parts.some(p=>p.bowelWallSource!==sourceId))return null;
 const selectedSources=new Set(parts.map(p=>p.bowelWallSource));
 const selected=new Set(choice.elements),same=(ids:string[])=>ids.length===selected.size&&ids.every(id=>selected.has(id));
 const identity=choice.id.replace(/^hierarchy:guest:ta98:/,'');
 const assembly=atlas.concepts.find(c=>c.id===identity);
 const fold=assembly&&/^SIW:(duodenum|jejunum|ileum):circular-folds$/.test(identity)?assembly:null;
 const sourceCells=atlas.parts.filter(p=>allowed.has(p.id)&&(sharedAssembly?p.caecocolicWallAssembly===assemblyId&&selectedSources.has(p.bowelWallSource):p.bowelWallSource===sourceId)&&!p.suppressed);
 const whole=same(sourceCells.map(p=>p.id))&&atlas.concepts.some(c=>same(c.elements));
 const aggregate=assembly&&/^(?:SIW:(?:duodenum|jejunum|ileum)|CW:colon|CAECOWALL:(?:female|male):(?:caecum|appendix|ascending|transverse|descending|sigmoid|rectum)):(?:mucosa|muscularis)$/.test(identity)?assembly:null;
 // The hierarchy's mucosa feature subtree also reaches the folds' submucosal
 // core. Include precisely that stored core, not the whole submucosal coat.
 const foldCore=identity.endsWith(':mucosa')?sourceCells.filter(p=>p.id.endsWith(':circular-folds:submucosa')).map(p=>p.id):[];
 const aggregateMembers=aggregate&&(same(aggregate.elements)||same([...aggregate.elements,...foldCore]));
 if(!whole&&(!fold||!same(fold.elements))&&!aggregateMembers)return null;
 const order=whole?layers:fold?[...layers.slice(4)].reverse():identity.endsWith(':muscularis')?layers.slice(2,4):layers.slice(parts.some(p=>p.id.endsWith(':circular-folds:submucosa'))?4:5);
 const steps=order.map(key=>{
  const cells=parts.filter(p=>p.caecocolicLayer===key||p.id.endsWith(':'+key)||key==='epithelium'&&p.id.endsWith(':circular-folds'));
  const tissue=cells.find(p=>!p.id.includes(':circular-folds'))??cells[0];
  const name=sharedAssembly?coatName(key,parts)+' of '+(selectedSources.size>1?'large intestine':tissue?.name.split(' of ').at(-1)??parts[0].caecocolicRegion):tissue?.name??key;
  return {name,ids:cells.map(p=>p.id)};
 });
 return steps.every(step=>step.ids.length)&&steps.flatMap(step=>step.ids).length===parts.length?steps:null;
}
export function wallPeelDepth(steps:WallPeelStep[],hidden:readonly string[]):number{
 const ids=new Set(hidden);let depth=0;for(const step of steps){if(!step.ids.every(id=>ids.has(id)))break;depth++;}return depth;
}
