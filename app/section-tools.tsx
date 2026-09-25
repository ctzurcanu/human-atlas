import {useMemo} from 'react';
import {Plus,X} from 'lucide-react';
import {Button} from '@/components/ui/button';
import {type Atlas,type SceneState} from './anatomy';
import {sectionSetBounds} from './section-set';
import {sectionAngles,sectionFraction,sectionPoint,type SectionAxis,type SectionState} from './section-plane';
import {sectionStack} from './section-stack';

type Props={atlas:Atlas;state:SceneState;setState:(update:(state:SceneState)=>SceneState)=>void;close:()=>void};
export default function SectionTools({atlas,state,setState,close}:Props){
 const sections=sectionStack(state),active=Math.min(state.activeSection??0,sections.length-1),section=sections[active];
 const bounds=useMemo(()=>sectionSetBounds(atlas,state),[atlas,state]);
 const angles=sectionAngles(section);
 const setStack=(next:SectionState[],index:number)=>setState(current=>({...current,sections:next,activeSection:index,section:next[index],camera:undefined}));
 const update=(change:Partial<SectionState>)=>setState(current=>{
  const stack=sectionStack(current),index=Math.min(current.activeSection??0,stack.length-1),previous=stack[index],next={...previous,...change};stack[index]=next;
  return {...current,sections:stack,activeSection:index,section:next,...(change.enabled===true&&!previous.enabled?{camera:undefined}:{})};
 });
 const chooseAxis=(axis:SectionAxis)=>setState(current=>{
  const stack=sectionStack(current),index=Math.min(current.activeSection??0,stack.length-1),previous=stack[index],anchor=sectionPoint(bounds,previous),presetAngles=axis==='oblique'&&previous.axis!=='oblique'?sectionAngles(previous):{};
  const next={...previous,...presetAngles,axis,enabled:true,flip:previous.enabled?previous.flip:true};next.position=sectionFraction(bounds,next,anchor);stack[index]=next;
  return {...current,sections:stack,activeSection:index,section:next,camera:undefined};
 });
 const changeAngle=(key:'azimuth'|'elevation',value:number)=>setState(current=>{
  const stack=sectionStack(current),index=Math.min(current.activeSection??0,stack.length-1),previous=stack[index],anchor=sectionPoint(bounds,previous);
  const next={...previous,...sectionAngles(previous),axis:'oblique' as const,[key]:value,enabled:true};next.position=sectionFraction(bounds,next,anchor);stack[index]=next;
  return {...current,sections:stack,activeSection:index,section:next,camera:undefined};
 });
 const add=()=>{if(sections.length===2)return;const next=[...sections,{...section,enabled:true,position:Math.min(.95,section.position+.12)}];setStack(next,1);};
 const remove=()=>{if(sections.length===1)return;const next=sections.filter((_,index)=>index!==active);setStack(next,0);};
 const positionLabel=section.axis==='axial'?['Superior','Inferior']:section.axis==='sagittal'?['Right','Left']:section.axis==='coronal'?['Anterior','Posterior']:['Start','End'];
 return <section className="section-panel glass" aria-label="Sections">
  <div className="panel-heading"><strong>Sections</strong><Button variant="ghost" className="icon-button" onClick={close} aria-label="Close sections"><X size={17}/></Button></div>
  <div className="section-preview"><canvas aria-label="3D model with section planes"/></div>
  <div className="section-tabs" role="tablist" aria-label="Sections">
   {sections.map((item,index)=><button type="button" role="tab" id={'section-tab-'+index} aria-controls="section-settings" aria-selected={active===index} className={active===index?'active':''} key={index} onClick={()=>setStack(sections,index)}>Section {index+1}{item.enabled?'':' · closed'}</button>)}
   {sections.length<2&&<button type="button" className="section-add icon-button" aria-label="Add second section" title="Add second section" onClick={add}><Plus size={16}/></button>}
   {sections.length===2&&<button type="button" className="section-remove icon-button" aria-label={'Remove section '+(active+1)} title="Remove this section" onClick={remove}><X size={15}/></button>}
  </div>
  <div className="section-scroll" role="tabpanel" id="section-settings" aria-labelledby={'section-tab-'+active}>
   <fieldset className="section-presets"><legend>Plane</legend>{(['axial','sagittal','coronal'] as const).map(axis=><button type="button" key={axis} aria-pressed={section.axis===axis} onClick={()=>chooseAxis(axis)}><strong>{axis[0].toUpperCase()+axis.slice(1)}</strong><span>{axis==='axial'?'Across':axis==='sagittal'?'Side to side':'Front to back'}</span></button>)}<button type="button" aria-pressed={section.axis==='oblique'} onClick={()=>chooseAxis('oblique')}><strong>Oblique</strong><span>Any angle</span></button></fieldset>
   {section.axis==='oblique'&&<div className="section-angle-controls"><label>Turn <output>{Math.round(angles.azimuth)}°</output><input type="range" min="-180" max="180" step="1" value={angles.azimuth} onChange={event=>changeAngle('azimuth',Number(event.target.value))}/></label><label>Tilt <output>{Math.round(angles.elevation)}°</output><input type="range" min="-90" max="90" step="1" value={angles.elevation} onChange={event=>changeAngle('elevation',Number(event.target.value))}/></label></div>}
   <div className="section-position"><label htmlFor="section-position">Position <output>{Math.round(section.position*100)}%</output></label><input id="section-position" type="range" min="0" max="1000" step="1" value={Math.round(section.position*1000)} onChange={event=>update({enabled:true,position:Number(event.target.value)/1000})}/><div><span>{positionLabel[0]}</span><span>{positionLabel[1]}</span></div></div>
   <label className="check-label"><input type="checkbox" checked={section.flip} onChange={event=>update({flip:event.target.checked})}/>Keep opposite side</label>
   <button type="button" className="section-toggle" onClick={()=>update({enabled:!section.enabled})}>{section.enabled?'Close section':'Open section'}</button>
  </div>
 </section>;
}
