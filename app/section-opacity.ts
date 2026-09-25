import type {Part,SceneState} from './anatomy';

// Match the opacity expression used by the main mesh shader in scene.tsx.
export function sectionCapOpacity(part:Pick<Part,'id'|'system'>,state:Pick<SceneState,'selected'|'contextOpacity'|'skinOpacity'>,sourceOpacity=1){
 if(state.selected.includes(part.id))return 1;
 const skin=part.system==='integumentary'||part.system==='cell-boundary'?state.skinOpacity??.1:1;
 return Math.max(0,Math.min(1,(state.contextOpacity??1)*skin*sourceOpacity));
}
