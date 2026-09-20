import type {Atlas,Part,SceneState,SystemId,View} from './anatomy';
export const REGIONS=['all','head-neck','torso','upper-right','upper-left','lower-right','lower-left'] as const;
export const REGION_NAMES=['Whole body','Head & neck','Torso & pelvis','Right arm','Left arm','Right leg','Left leg'];
export const VIEWS:View[]=['three-quarter','front','back','side','right','superior','inferior'];
export const PEEL_LAYERS:SystemId[][]=[['integumentary','regions'],['fascia'],['muscular'],['connective'],['digestive','respiratory','urinary','reproductive','endocrine'],['arterial','venous','lymphatic'],['skeletal']];
export const PEEL_NAMES=['Intact','Skin removed','Fascia removed','Muscles removed','Connective tissue removed','Organs removed','Vessels removed','Skeleton removed'];
export function inRegion(p:Part,region='all'){
 if(region==='all')return true;if(p.system==='integumentary')return false;
 if(p.regions?.length)return p.regions.includes(region);
 const x=(p.bounds[0][0]+p.bounds[1][0])/2,y=(p.bounds[0][1]+p.bounds[1][1])/2;
 return region==='head-neck'?y>1.35:region==='torso'?y>=.75&&y<=1.4&&Math.abs(x)<.23:region.startsWith('upper')?y>.7&&Math.abs(x)>.18&&(region.endsWith('left')?x>0:x<0):y<.85&&(region.endsWith('left')?x>0:x<0);
}
export function partVisible(p:Part,s:SceneState){
 if(s.selected.includes(p.id))return true;
 if(p.suppressed||s.isolate||s.hidden?.includes(p.id)||!s.visible.includes(p.system)||!inRegion(p,s.region))return false;
 return !PEEL_LAYERS.slice(0,s.peel??0).flat().includes(p.system);
}
export function resolveSelection(atlas:Atlas,terms:string[]){
 return [...new Set(terms.flatMap(term=>atlas.concepts.find(c=>c.id===term||c.id==='DETAIL:'+term||c.name.toLowerCase()===term.toLowerCase())?.elements??(atlas.parts.some(p=>p.id===term)?[term]:[])))];
}
const number=(text:string|null,min:number,max:number,fallback:number)=>{if(text===null)return fallback;const n=Number(text);return Number.isFinite(n)?Math.min(max,Math.max(min,n)):fallback;};
export function readViewUrl(search:string,atlas:Atlas,base:SceneState):SceneState{
 const q=new URLSearchParams(search),view=q.get('view') as View,region=q.get('region')??'all',axis=q.get('cut'),cam=q.get('camera')?.split(',').map(Number);
 return {...base,focus:number(q.get('focus'),0,1,0),selected:resolveSelection(atlas,q.getAll('select')),view:VIEWS.includes(view)?view:base.view,region:REGIONS.includes(region as never)?region:'all',contextOpacity:number(q.get('context'),0,1,1),skinOpacity:number(q.get('skin'),0,1,base.skinOpacity??.1),peel:Math.round(number(q.get('peel'),0,7,0)),hidden:resolveSelection(atlas,q.getAll('hide')),isolate:q.get('isolate')==='1',explode:number(q.get('explode'),0,1,0),labels:q.get('labels')!=='0',visible:q.has('layers')?(q.get('layers')!.split(',').filter(id=>atlas.parts.some(p=>p.system===id)) as SystemId[]):base.visible,section:{enabled:['axial','sagittal','coronal'].includes(axis??''),axis:axis==='sagittal'||axis==='coronal'?axis:'axial',position:number(q.get('slice'),0,1,.5),flip:q.get('flip')==='1'},camera:cam&&[6,8].includes(cam.length)&&cam.every(n=>Number.isFinite(n)&&Math.abs(n)<1000)&&Math.hypot(cam[0]-cam[3],cam[1]-cam[4],cam[2]-cam[5])>.001?cam:undefined};
}
export function viewUrl(base:string,model:string,s:SceneState,camera?:number[]){
 const url=new URL(base);url.search='';url.hash='';const q=url.searchParams;q.set('model',model);q.set('view',s.view);s.selected.forEach(id=>q.append('select',id));q.set('layers',s.visible.join(','));q.set('context',String(s.contextOpacity??1));q.set('skin',String(s.skinOpacity??.1));q.set('region',s.region??'all');q.set('peel',String(s.peel??0));(s.hidden??[]).forEach(id=>q.append('hide',id));if(s.isolate)q.set('isolate','1');if(s.explode)q.set('explode',String(s.explode));if(s.labels===false)q.set('labels','0');if(s.section?.enabled){q.set('cut',s.section.axis);q.set('slice',String(s.section.position));if(s.section.flip)q.set('flip','1');}if(camera)q.set('camera',camera.map(n=>Number(n.toFixed(6))).join(','));return url.href;
}
