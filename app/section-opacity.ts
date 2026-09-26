import type {Part,SceneState} from './anatomy';
import {isSkinPart} from './depth-layers';

// Match the opacity expression used by the main mesh shader in scene.tsx.
export function sectionCapOpacity(part:Part,state:Pick<SceneState,'selected'|'contextOpacity'|'skinOpacity'> & Partial<Pick<SceneState,'section'|'sections'>>,sourceOpacity=1){
 const skin=isSkinPart(part)?state.skinOpacity??.1:1;
 const sectionOn=state.sections?.length?state.sections.some(section=>section.enabled):state.section?.enabled;
 const context=sectionOn||!state.selected.length||state.selected.includes(part.id)?1:state.contextOpacity??1;
 return Math.max(0,Math.min(1,context*skin*sourceOpacity));
}
