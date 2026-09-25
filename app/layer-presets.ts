import type {SceneState,SystemId} from './anatomy';

export const LAYER_PRESETS=[
 {id:'all',label:'All'},
 {id:'skin',label:'Skin'},
 {id:'muscles',label:'Muscles'},
 {id:'vessels',label:'Vessels'},
 {id:'skeleton',label:'Skeleton'},
 {id:'organs',label:'Organs'},
] as const;

export type LayerPreset=typeof LAYER_PRESETS[number]['id'];
const presetSystems:Record<Exclude<LayerPreset,'all'>,SystemId[]>={
 skin:['integumentary'],
 muscles:['muscular'],
 vessels:['arterial','venous'],
 skeleton:['skeletal'],
 organs:['cardiac','respiratory','digestive','urinary','endocrine','reproductive'],
};

export function presetLayers(preset:LayerPreset,available:SystemId[]):SystemId[]{
 const desired=preset==='all'?available.filter(id=>!['pregnancy','attachments','regions','schematic'].includes(id)):presetSystems[preset];
 return desired.filter(id=>available.includes(id));
}

export function applyLayerPreset(state:SceneState,preset:LayerPreset,available:SystemId[]):SceneState{
 return {...state,visible:presetLayers(preset,available),skinOpacity:preset==='skin'?1:preset==='all'?.1:0,
  selected:[],isolate:false,contextOpacity:1,peel:0,depthHidden:[],explode:0,rotate:false,camera:undefined};
}
