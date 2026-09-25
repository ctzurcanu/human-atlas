import type {SceneState} from './anatomy';
import type {SectionState} from './section-plane';

export const defaultSection:SectionState={enabled:false,axis:'axial',position:.38,flip:true,azimuth:35,elevation:30};
export function sectionStack(state:SceneState){return state.sections?.length?state.sections.slice(0,2):[state.section??defaultSection];}
export function enabledSections(state:SceneState){return sectionStack(state).map((section,index)=>({section,index})).filter(({section})=>section.enabled);}
