import {isSurfaceSystem,type Atlas,type Part,type SceneState,type SystemId,type View} from './anatomy';
import {DEPTH_LAYERS,depthLayerFor} from './depth-layers';
export const REGIONS=['all','head-neck','torso','upper-right','upper-left','lower-right','lower-left'] as const;
export const REGION_NAMES=['Whole body','Head & neck','Torso & pelvis','Right arm','Left arm','Right leg','Left leg'];
export const VIEWS:View[]=['three-quarter','front','back','side','right','superior','inferior'];
const LEGACY_PEEL_LAYERS=[['skin'],['investing-fascia'],['superficial-muscles','second-muscles','intermediate-muscles','deep-muscles','deepest-muscles'],['visceral-coverings','ligaments'],['anterior-organs','deep-organs'],['superficial-veins','deep-vessels','lymphatic'],['thoracic-bones','limb-bones','other-bones','spine','skull']];
const visibilityCache=new WeakMap<SceneState,{selected:Set<string>;hidden:Set<string>;visible:Set<SystemId>;depthHidden:Set<string>;surface:boolean}>();
function visibilitySets(s:SceneState){
 let sets=visibilityCache.get(s);
 if(!sets){sets={selected:new Set(s.selected),hidden:new Set(s.hidden??[]),visible:new Set(s.visible),depthHidden:new Set(s.depthHidden??[]),surface:s.visible.some(isSurfaceSystem)};visibilityCache.set(s,sets);}
 return sets;
}
export function inRegion(p:Part,region='all'){
 if(region==='all')return true;if(p.system==='integumentary')return false;
 if(p.regions?.length)return p.regions.includes(region);
 const x=(p.bounds[0][0]+p.bounds[1][0])/2,y=(p.bounds[0][1]+p.bounds[1][1])/2;
 return region==='head-neck'?y>1.35:region==='torso'?y>=.75&&y<=1.4&&Math.abs(x)<.23:region.startsWith('upper')?y>.7&&Math.abs(x)>.18&&(region.endsWith('left')?x>0:x<0):y<.85&&(region.endsWith('left')?x>0:x<0);
}
export function partVisible(p:Part,s:SceneState){
 const sets=visibilitySets(s);
 if(sets.selected.has(p.id))return true;
 if(p.suppressed||s.isolate||sets.hidden.has(p.id)||sets.depthHidden.has(depthLayerFor(p))||!sets.visible.has(p.system)||!inRegion(p,s.region))return false;
 // Some primary-source internal meshes protrude beyond the outer reference.
 // An intact, fully opaque surface should conceal unselected internal context.
 if(!isSurfaceSystem(p.system)&&sets.surface&&!sets.depthHidden.has('skin')&&(s.skinOpacity??0)>=.999&&(s.explode??0)<.001)return false;
 return true;
}
export function resolveSelection(atlas:Atlas,terms:string[]){
 return [...new Set(terms.flatMap(term=>{
  const sourceName=term.startsWith('DETAIL:')?term.slice(7):term;
  const concept=atlas.concepts.find(c=>c.id===term||c.id==='DETAIL:'+term||c.id==='ZA:'+sourceName||c.name.toLowerCase()===sourceName.toLowerCase());
  return concept?.elements??atlas.parts.filter(p=>p.id===term||p.sourceId===sourceName).map(p=>p.id);
 }))];
}
const number=(text:string|null,min:number,max:number,fallback:number)=>{if(text===null)return fallback;const n=Number(text);return Number.isFinite(n)?Math.min(max,Math.max(min,n)):fallback;};
export function readViewUrl(search:string,atlas:Atlas,base:SceneState):SceneState{
 const q=new URLSearchParams(search),view=q.get('view') as View,region=q.get('region')??'all',axis=q.get('cut'),cam=q.get('camera')?.split(',').map(Number);
 const legacyDepth=LEGACY_PEEL_LAYERS.slice(0,Math.round(number(q.get('peel'),0,7,0))).flat();
 const depthHidden=atlas.scope==='cell'?[]:(q.has('depth')?q.get('depth')!.split(','):legacyDepth).filter(id=>DEPTH_LAYERS.some(layer=>layer.id===id));
 return {...base,focus:number(q.get('focus'),0,1,0),selected:resolveSelection(atlas,q.getAll('select')),view:VIEWS.includes(view)?view:base.view,region:atlas.scope==='cell'?'all':REGIONS.includes(region as never)?region:'all',contextOpacity:number(q.get('context'),0,1,1),skinOpacity:number(q.get('skin'),0,1,base.skinOpacity??.1),peel:0,depthHidden,hidden:resolveSelection(atlas,q.getAll('hide')),isolate:q.get('isolate')==='1',explode:number(q.get('explode'),0,1,0),rotate:q.get('rotate')==='1',labels:q.get('labels')!=='0',visible:q.has('layers')?(q.get('layers')!.split(',').filter(id=>atlas.parts.some(p=>p.system===id)) as SystemId[]):base.visible,section:{enabled:['axial','sagittal','coronal'].includes(axis??''),axis:axis==='sagittal'||axis==='coronal'?axis:'axial',position:number(q.get('slice'),0,1,.5),flip:q.get('flip')==='1'},camera:cam&&[6,8].includes(cam.length)&&cam.every(n=>Number.isFinite(n)&&Math.abs(n)<1000)&&Math.hypot(cam[0]-cam[3],cam[1]-cam[4],cam[2]-cam[5])>.001?cam:undefined};
}
export function viewUrl(base:string,model:string,s:SceneState,camera?:number[]){
 const url=new URL(base);url.search='';url.hash='';const q=url.searchParams;q.set('model',model);q.set('view',s.view);s.selected.forEach(id=>q.append('select',id));q.set('layers',s.visible.join(','));q.set('context',String(s.contextOpacity??1));q.set('skin',String(s.skinOpacity??.1));q.set('region',s.region??'all');if(s.depthHidden?.length)q.set('depth',s.depthHidden.join(','));(s.hidden??[]).forEach(id=>q.append('hide',id));if(s.isolate)q.set('isolate','1');if(s.explode)q.set('explode',String(s.explode));if(s.rotate)q.set('rotate','1');if(s.labels===false)q.set('labels','0');if(s.section?.enabled){q.set('cut',s.section.axis);q.set('slice',String(s.section.position));if(s.section.flip)q.set('flip','1');}if(camera)q.set('camera',camera.map(n=>Number(n.toFixed(6))).join(','));return url.href;
}
