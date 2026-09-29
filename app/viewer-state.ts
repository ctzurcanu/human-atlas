import {isSurfaceSystem,type Atlas,type Part,type SceneState,type SystemId,type View} from './anatomy';
import {packPartIds,partFingerprint,unpackPartIds} from './view-part-codes';
import {DEPTH_LAYERS,depthLayerFor,isSkinPart} from './depth-layers';
import {completedDepthLayers,depthLayerOpacity,partLayerOpacity} from './depth-control';
import type {SectionState} from './section-plane';
import {validCamera} from '../shared/camera-frame.mjs';
import {compactViewUrl,viewParameters} from './view-url';
export const REGIONS=['all','head-neck','torso','upper-right','upper-left','lower-right','lower-left'] as const;
export const REGION_NAMES=['Whole body','Head & neck','Torso & pelvis','Right arm','Left arm','Right leg','Left leg'];
export const VIEWS:View[]=['three-quarter','front','back','side','right','superior','inferior'];
// Fixed URL bit positions. Append new systems; never reorder existing entries.
const URL_SYSTEMS:readonly SystemId[]=['skeletal','muscular','cardiac','sensory','arterial','venous','nervous','respiratory','digestive','urinary','lymphatic','endocrine','reproductive','integumentary','pregnancy','fascia','schematic','attachments','regions','connective','cell-boundary','cell-nucleus','cell-endomembrane','cell-mitochondria','cell-cytoskeleton','cell-vesicles','cell-protein-synthesis','cell-division'];
const LEGACY_PEEL_LAYERS=[['skin'],['investing-fascia'],['superficial-muscles','second-muscles','intermediate-muscles','deep-muscles','deepest-muscles'],['visceral-coverings','ligaments'],['anterior-organs','deep-organs'],['superficial-veins','deep-vessels','lymphatic'],['thoracic-bones','limb-bones','other-bones','spine','skull']];
const visibilityCache=new WeakMap<SceneState,{selected:Set<string>;hidden:Set<string>;visible:Set<SystemId>;depthHidden:Set<string>;surface:boolean}>();
function visibilitySets(s:SceneState){
 let sets=visibilityCache.get(s);
 if(!sets){sets={selected:new Set(s.selected),hidden:new Set(s.hidden??[]),visible:new Set(s.visible),depthHidden:new Set(s.depthHidden??[]),surface:s.visible.some(isSurfaceSystem)};visibilityCache.set(s,sets);}
 return sets;
}
export function inRegion(p:Part,region='all'){
 if(region==='all')return true;
 // Imported reference regions have separately tagged skin. Keep the single
 // whole-body derived surface out of regional views, but retain local skin.
 if(isSkinPart(p)&&/^(skin(?: of body)?|body surface)$/i.test(p.name)&&!p.regions?.length)return false;
 if(p.regions?.length)return p.regions.includes(region);
 const x=(p.bounds[0][0]+p.bounds[1][0])/2,y=(p.bounds[0][1]+p.bounds[1][1])/2;
 return region==='head-neck'?y>1.35:region==='torso'?y>=.75&&y<=1.4&&Math.abs(x)<.23:region.startsWith('upper')?y>.7&&Math.abs(x)>.18&&(region.endsWith('left')?x>0:x<0):y<.85&&(region.endsWith('left')?x>0:x<0);
}
function inSectionRegion(p:Part,region='all'){
 if(inRegion(p,region))return true;
 // A regional cut must include a neighboring structure where its actual
 // geometry enters the cut volume (for example gluteal muscle above the hip).
 // Source chapter tags alone omit that tissue and leave a false empty space.
 if(region==='torso'&&!isSkinPart(p)&&p.regions?.some(item=>item.startsWith('lower-'))){
  const [min,max]=p.bounds;
  return min[0]<=.23&&max[0]>=-.23&&min[1]<=1.4&&max[1]>=.75&&min[2]<=.17&&max[2]>=-.17;
 }
 return false;
}
export function partVisible(p:Part,s:SceneState){
 if(p.suppressed)return false;
 const sets=visibilitySets(s);
 if(s.depth!==undefined&&sets.depthHidden.has(depthLayerFor(p)))return false;
 if(sets.selected.has(p.id))return true;
 if(s.isolate||sets.hidden.has(p.id)||sets.depthHidden.has(depthLayerFor(p))||!sets.visible.has(p.system)||!inRegion(p,s.region))return false;
 if(partLayerOpacity(p,s)<=0)return false;
 // Some primary-source internal meshes protrude beyond the outer reference.
 // An intact, fully opaque surface should conceal unselected internal context.
 if(!isSkinPart(p)&&sets.surface&&!sets.depthHidden.has('skin')&&depthLayerOpacity('skin',s)>=.999&&(s.explode??0)<.001&&!s.sections?.some(section=>section.enabled)&&!s.section?.enabled)return false;
 return true;
}
// Sections use the layer checkboxes as their source set, including pieces that
// lie entirely on either side of a cut. Selection alone must not add a piece.
export function sectionPartVisible(p:Part,s:SceneState){
 const sets=visibilitySets(s);
 if(p.suppressed||s.depth!==undefined&&sets.depthHidden.has(depthLayerFor(p)))return false;
 if(s.isolate)return sets.selected.has(p.id);
 // A single whole-body skin mesh is still one checked item in a section.
 // Excluding it by region made 100% skin appear absent in every regional cut.
 const wholeSkin=isSkinPart(p)&&/^(skin(?: of body)?|body surface)$/i.test(p.name)&&!p.regions?.length;
 return sets.visible.has(p.system)&&!sets.hidden.has(p.id)&&!sets.depthHidden.has(depthLayerFor(p))&&(inSectionRegion(p,s.region)||wholeSkin)&&partLayerOpacity(p,s)>0;
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
 const q=viewParameters(search),view=q.get('view') as View,region=q.get('region')??'all',cam=(q.get('frame')??q.get('camera'))?.split(',').map(Number);
 const match=q.has('p')&&q.get('p')===partFingerprint(atlas);
 const selected=q.has('select')?resolveSelection(atlas,q.getAll('select')):match&&q.has('s')?unpackPartIds(q.get('s')!,atlas,false):[];
 const hidden=q.has('hide')?resolveSelection(atlas,q.getAll('hide')):match&&q.has('h')?unpackPartIds(q.get('h')!,atlas,true):[];
 const mask=q.get('l'),maskNumber=mask&&/^[0-9a-z]{1,7}$/.test(mask)?parseInt(mask,36):NaN;
 const visible=q.has('layers')?(q.get('layers')!.split(',').filter(id=>atlas.parts.some(p=>p.system===id)) as SystemId[]):Number.isSafeInteger(maskNumber)&&maskNumber>=0&&maskNumber<2**URL_SYSTEMS.length?URL_SYSTEMS.filter((_,index)=>(maskNumber&2**index)!==0):base.visible;
 const readSection=(suffix:string):SectionState=>{const axis=q.get('cut'+suffix);return {enabled:['axial','sagittal','coronal','oblique'].includes(axis??''),axis:axis==='sagittal'||axis==='coronal'||axis==='oblique'?axis:'axial',position:number(q.get('slice'+suffix),0,1,.38),flip:axis?q.get('flip'+suffix)==='1':true,...(axis==='oblique'?{azimuth:number(q.get('azimuth'+suffix),-180,180,35),elevation:number(q.get('elevation'+suffix),-90,90,30)}:{})};};
 const sections=[readSection('')];if(q.get('sections')==='2'||q.has('cut2'))sections.push(readSection('2'));
 const activeSection=sections.length===2&&q.get('sectionTab')==='2'?1:0;
 const legacyDepth=LEGACY_PEEL_LAYERS.slice(0,Math.round(number(q.get('peel'),0,7,0))).flat();
 const depth=atlas.scope!=='cell'&&q.has('dp')?number(q.get('dp'),0,1,0):undefined;
 const depthHidden=atlas.scope==='cell'?[]:(q.has('depth')?q.get('depth')!.split(','):depth!==undefined?completedDepthLayers(depth):legacyDepth).filter(id=>DEPTH_LAYERS.some(layer=>layer.id===id));
 return {...base,guestExtensions:q.get('bioExpand')?.split(',').filter(id=>/^HGNC:\d+$/.test(id)),guestView:q.get('bioView')?.slice(0,80)??undefined,guestQuery:q.get('bio')?.slice(0,200)??undefined,focus:number(q.get('focus'),0,1,0),selected,view:VIEWS.includes(view)?view:base.view,region:atlas.scope==='cell'?'all':REGIONS.includes(region as never)?region:'all',contextOpacity:number(q.get('context'),0,1,1),skinOpacity:number(q.get('skin'),0,1,base.skinOpacity??.1),peel:0,depth,depthHidden,hidden,isolate:q.get('isolate')==='1',explode:number(q.get('explode'),0,1,0),rotate:q.get('rotate')==='1',labels:q.get('labels')!=='0',visible,sections,activeSection,section:sections[activeSection],camera:cam&&validCamera(cam)?cam:undefined};
}
export function viewUrl(base:string,model:string,s:SceneState,camera?:number[],atlas?:Atlas){
 const url=new URL(base);url.search='';url.hash='';const q=url.searchParams;q.set('model',model);if(s.view!=='three-quarter')q.set('view',s.view);if(s.guestView)q.set('bioView',s.guestView);if(s.guestExtensions?.length)q.set('bioExpand',s.guestExtensions.join(','));if(s.guestQuery)q.set('bio',s.guestQuery);
 const packedSelection=atlas&&s.selected.length?packPartIds(s.selected,atlas,false):undefined;
 const packedHidden=atlas&&s.hidden?.length?packPartIds(s.hidden,atlas,true):undefined;
 if(packedSelection||packedHidden)q.set('p',partFingerprint(atlas!));
 if(packedSelection)q.set('s',packedSelection);else s.selected.forEach(id=>q.append('select',id));
 if(packedHidden)q.set('h',packedHidden);else (s.hidden??[]).forEach(id=>q.append('hide',id));
 if(atlas&&s.visible.every(id=>URL_SYSTEMS.includes(id))){const mask=URL_SYSTEMS.reduce((value,id,index)=>s.visible.includes(id)?value+2**index:value,0);q.set('l',mask.toString(36));}else q.set('layers',s.visible.join(','));
 if((s.contextOpacity??1)!==1)q.set('context',String(s.contextOpacity));if((s.skinOpacity??.1)!==.1)q.set('skin',String(s.skinOpacity));if(s.region&&s.region!=='all')q.set('region',s.region);if(s.depth!==undefined){q.set('dp',String(Number(s.depth.toFixed(6))));const expected=new Set<string>(completedDepthLayers(s.depth)),hidden=new Set(s.depthHidden??[]);if(expected.size!==hidden.size||[...expected].some(id=>!hidden.has(id)))q.set('depth',[...hidden].join(','));}else if(s.depthHidden?.length)q.set('depth',s.depthHidden.join(','));if(s.isolate)q.set('isolate','1');if(s.explode)q.set('explode',String(s.explode));if(s.rotate)q.set('rotate','1');if(s.labels===false)q.set('labels','0');
 const sections=s.sections?.length?s.sections.slice(0,2):[s.section].filter((section):section is NonNullable<typeof section>=>!!section);
 const writeSection=(section:SectionState|undefined,suffix:string)=>{if(!section?.enabled)return;q.set('cut'+suffix,section.axis);q.set('slice'+suffix,String(section.position));if(section.flip)q.set('flip'+suffix,'1');if(section.axis==='oblique'){q.set('azimuth'+suffix,String(section.azimuth??35));q.set('elevation'+suffix,String(section.elevation??30));}};
 writeSection(sections[0],'');if(sections.length===2){q.set('sections','2');writeSection(sections[1],'2');if(s.activeSection===1)q.set('sectionTab','2');}
 if(camera)q.set(camera.length===12?'frame':'camera',camera.map(n=>Number(n.toFixed(6))).join(','));return compactViewUrl(url.href);
}
