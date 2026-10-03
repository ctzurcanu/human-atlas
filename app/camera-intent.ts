import type {SceneState} from './anatomy';

/** Display/selection changes must never imply permission to move the camera. */
export function explicitCameraAction(state:SceneState,previous:{view:string;reset:number;focus:number},newCamera:boolean){
 return state.view!==previous.view||state.reset!==previous.reset||!!state.focus&&state.focus!==previous.focus||newCamera;
}
