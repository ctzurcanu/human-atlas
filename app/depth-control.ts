import {isSurfaceSystem,type Part,type SceneState} from './anatomy';
import {DEPTH_LAYERS,depthLayerFor,isSkinPart,type DepthLayerId} from './depth-layers';

type DepthState=Pick<SceneState,'depth'|'skinOpacity'|'depthHidden'|'visible'>;
type OpacityState=Pick<SceneState,'depth'|'skinOpacity'>;
const layerIndexes=new Map<string,number>(DEPTH_LAYERS.map((layer,index)=>[layer.id,index]));
const clamp=(value:number)=>Math.max(0,Math.min(1,value));
const layerProgress=(amount:number)=>{
 const value=clamp(amount)*DEPTH_LAYERS.length,nearest=Math.round(value);
 return Math.abs(value-nearest)<2e-5?nearest:value;
};
export function completedDepthLayers(amount:number){return DEPTH_LAYERS.slice(0,Math.floor(layerProgress(amount))).map(layer=>layer.id);}

/** Legacy views keep their existing surface opacity and manual layer choices. */
export function depthPosition(state:DepthState){
 if(state.depth!==undefined)return clamp(state.depth);
 let complete=0;
 while(complete<DEPTH_LAYERS.length&&(state.depthHidden?.includes(DEPTH_LAYERS[complete].id)||complete===0&&(!state.visible.some(system=>isSurfaceSystem(system)||system==='regions')||(state.skinOpacity??.1)<=0)))complete++;
 return (complete===0?1-clamp(state.skinOpacity??.1):complete)/DEPTH_LAYERS.length;
}
/** Completed layers use depthHidden; only the layer between the ticks fades. */
export function depthLayerOpacity(id:DepthLayerId,state:OpacityState){
 if(state.depth===undefined)return id==='skin'?state.skinOpacity??.1:1;
 const progress=layerProgress(state.depth),index=layerIndexes.get(id)!;
 return index===Math.floor(progress)?1-(progress-index):1;
}
export function partLayerOpacity(part:Part,state:OpacityState){
 return state.depth===undefined?(isSkinPart(part)?state.skinOpacity??.1:1):depthLayerOpacity(depthLayerFor(part),state);
}
/** Replace the previous slider's hidden prefix, retaining independent layer hides. */
export function setDepthPosition(state:SceneState,amount:number,parts:readonly Part[]):SceneState{
 const depth=clamp(amount),previous=new Set<string>(completedDepthLayers(depthPosition(state)));
 const depthHidden=[...new Set([...(state.depthHidden??[]).filter(id=>!previous.has(id)),...completedDepthLayers(depth)])];
 const surfaces=parts.filter(part=>!part.suppressed&&isSkinPart(part)).map(part=>part.system);
 return {...state,depth,depthHidden,skinOpacity:1,visible:[...new Set([...state.visible,...surfaces])]};
}
