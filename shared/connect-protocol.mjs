export const CONNECT_PATH = '/atlas-connect';
import {validCamera} from './camera-frame.mjs';
export const MAX_MESSAGE_BYTES = 1_048_576;
export const MODELS = ['male-full','male-detail','male','female','embryo','cell','local-male','local-female','local-reference'];
const object = value => !!value && typeof value === 'object' && !Array.isArray(value);
const text = (value,max=2048) => typeof value === 'string' && value.length <= max;
const strings = (value,max=10000) => Array.isArray(value) && value.length <= max && value.every(item=>text(item,256));
const number = (value,min=-10000,max=10000) => typeof value === 'number' && Number.isFinite(value) && value >= min && value <= max;
const vector = (value,length) => Array.isArray(value) && value.length === length && value.every(item=>number(item));
const section = value => object(value) && typeof value.enabled === 'boolean' && ['axial','sagittal','coronal','oblique'].includes(value.axis) && number(value.position,0,1) && typeof value.flip === 'boolean' && (value.azimuth===undefined||number(value.azimuth,-180,180)) && (value.elevation===undefined||number(value.elevation,-90,90));
const terminology = value => object(value) && ['ta98','tha','fma','latin','ontology'].every(key=>value[key]===undefined||value[key]===null||text(value[key],512));
const choice = value => object(value) && text(value.id,256) && text(value.name,512) && strings(value.elements) && typeof value.group==='boolean' && (value.terminology===undefined||terminology(value.terminology));
export function validPose(value){
 return object(value) && MODELS.includes(value.model) && validCamera(value.camera) && vector(value.up,3) && vector(value.rotation,4);
}
export function validSnapshot(value){
 if(!object(value)||value.version!==1||!MODELS.includes(value.model)||!text(value.hierarchy,256)||!strings(value.guestSources,16))return false;
 if(!value.guestSources.every(source=>{try{return ['https:','http:'].includes(new URL(source).protocol);}catch{return false;}}))return false;
 const state=value.state,ui=value.ui;
 if(!object(state)||!strings(state.selected)||!strings(state.visible,64)||!number(state.explode,0,1)||typeof state.isolate!=='boolean'||!['three-quarter','front','back','side','right','superior','inferior'].includes(state.view)||typeof state.rotate!=='boolean'||!number(state.reset,0,Number.MAX_SAFE_INTEGER))return false;
 if(['labels','selectionGroup','inspectorOpen'].some(key=>state[key]!==undefined&&typeof state[key]!=='boolean'))return false;
 if(['focus','peel'].some(key=>state[key]!==undefined&&!number(state[key],0,10000)))return false;
 if(['contextOpacity','skinOpacity'].some(key=>state[key]!==undefined&&!number(state[key],0,1)))return false;
 if(state.region!==undefined&&!['all','head-neck','torso','upper-right','upper-left','lower-right','lower-left'].includes(state.region))return false;
 if(['hidden','depthHidden'].some(key=>state[key]!==undefined&&!strings(state[key])))return false;
 if(state.camera!==undefined&&!validCamera(state.camera))return false;
 if(state.section!==undefined&&!section(state.section))return false;
 if(state.sections!==undefined&&(!Array.isArray(state.sections)||state.sections.length<1||state.sections.length>2||!state.sections.every(section)))return false;
 if(state.activeSection!==undefined&&!number(state.activeSection,0,1))return false;
 if(!object(ui)||!['search','advanced','sections',null].includes(ui.panel)||!['layers','search','advanced','sections','details'].includes(ui.frontPanel)||['details','layersVisible','mobileLayersOpen','addSelection'].some(key=>typeof ui[key]!=='boolean')||!text(ui.query,2000))return false;
 if(ui.expanded!==undefined&&(!Array.isArray(ui.expanded)||ui.expanded.length>1000||!ui.expanded.every(item=>object(item)&&text(item.key,512)&&typeof item.expanded==='boolean')))return false;
 if(value.choice!==null&&(!choice(value.choice)||value.choice.children!==undefined&&(!Array.isArray(value.choice.children)||value.choice.children.length>1000||!value.choice.children.every(choice))))return false;
 if(value.relationshipExpansion!==undefined&&value.relationshipExpansion!==null){const expansion=value.relationshipExpansion;if(!object(expansion)||!text(expansion.anchorId,256)||!strings(expansion.anchorIds)||!text(expansion.selectedKey,500000)||!strings(expansion.frontier)||typeof expansion.exhausted!=='boolean'||!number(expansion.depth,0,10000))return false;}
 if(!strings(value.covering)||value.pose!==null&&!validPose(value.pose))return false;
 if(value.advanced!==null){
  const advanced=value.advanced,doc=advanced?.document;
  if(!object(advanced)||!object(doc)||doc.schemaVersion!==1||!['slides','quiz','layers','connect'].includes(doc.type)||!Array.isArray(doc.views)||doc.views.length>500)return false;
  if(!doc.views.every(view=>object(view)&&text(view.id,256)&&text(view.name,512)&&text(view.url,100000)&&text(view.model,64)&&text(view.addedAt,64)))return false;
  if(advanced.slideIndex!==null&&(!Number.isInteger(advanced.slideIndex)||advanced.slideIndex<0||advanced.slideIndex>=doc.views.length))return false;
 }
 return true;
}
