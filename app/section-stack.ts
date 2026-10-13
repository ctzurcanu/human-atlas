import type {SceneState} from './anatomy';
import type {SectionState} from './section-plane';

export const defaultSection:SectionState={enabled:false,axis:'axial',position:.38,flip:true,azimuth:35,elevation:30};
export function sectionStack(state:SceneState){return state.sections?.length?state.sections.slice(0,2):[state.section??defaultSection];}
export function enabledSections(state:SceneState){return sectionStack(state).map((section,index)=>({section,index})).filter(({section})=>section.enabled);}

/** Close every cut without changing its settings or the current camera. */
export function closeSections(state:SceneState):SceneState{
 const sections=sectionStack(state).map(section=>({...section,enabled:false}));
 const activeSection=Math.max(0,Math.min(state.activeSection??0,sections.length-1));
 return {...state,sections,activeSection,section:sections[activeSection]};
}

/** Moving a live plane leaves anatomy, opacity and camera state unchanged. */
export function sectionPlaneOnlyChange(previous:SceneState|null,next:SceneState){
 if(!previous)return false;
 const ignored=new Set(['section','sections','activeSection']);
 if([...new Set([...Object.keys(previous),...Object.keys(next)])].some(key=>!ignored.has(key)&&previous[key as keyof SceneState]!==next[key as keyof SceneState]))return false;
 return enabledSections(previous).map(cut=>cut.index).join(',')===enabledSections(next).map(cut=>cut.index).join(',');
}
