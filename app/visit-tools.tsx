import {useEffect,useMemo,useState,useRef} from 'react';
import * as T from 'three';
import {Button} from '@/components/ui/button';
import {X,Play,Pause,RotateCcw} from 'lucide-react';
import type {Atlas} from './anatomy';
import type {AlimentaryVisit} from './alimentary-visit';
import type {VisitState} from './model-visit';
import {visitRegions} from './visit-map';
type Props={atlas:Atlas;progressRef:{current:number};alimentary?:AlimentaryVisit;value:VisitState;change:(update:(current:VisitState)=>VisitState)=>void;close:()=>void};
export default function VisitTools({atlas,progressRef,alimentary,value,change,close}:Props){
 const listRef=useRef<HTMLOListElement>(null);
 const bounds=useMemo(()=>{const b=new T.Box3();for(const p of atlas.parts)if(!p.suppressed)b.union(new T.Box3(new T.Vector3().fromArray(p.bounds[0]),new T.Vector3().fromArray(p.bounds[1])));return b;},[atlas]);
 const [progress,setProgress]=useState(progressRef.current);
 useEffect(()=>{const timer=setInterval(()=>setProgress(progressRef.current),80);return()=>clearInterval(timer);},[progressRef]);
 const regions=useMemo(()=>visitRegions(bounds,value.path,alimentary),[alimentary,value.path,bounds]);
 const seek=(position:number,reorient=false)=>{progressRef.current=position;setProgress(position);change(s=>({...s,seek:position,seekId:s.seekId+1,...(reorient?{restart:s.restart+1,lookAlong:value.path==='alimentary'?true:s.lookAlong}:{})}));};
 const gradient=regions.map((r,i)=>`${r.color} ${r.start*100}%, ${r.color} ${(regions[i+1]?.start??1)*100}%`).join(',');
 const current=Math.max(0,regions.reduce((index,r,i)=>r.start<=progress?i:index,-1));
 useEffect(()=>{const list=listRef.current,item=list?.children[current] as HTMLElement|undefined;if(list&&item)list.scrollTo({top:item.offsetTop-list.offsetTop-(list.clientHeight-item.offsetHeight)/2,behavior:'smooth'});},[current,value.path]);
 const positionText=value.path==='orbit'?`${Math.round(progress*360)}°`:value.path==='alimentary'?`${Math.round(progress*100)}%`:`${Math.round(progress*bounds.getSize(new T.Vector3()).y*100)} cm`;
 return <section className="detail-sheet visit-panel glass" aria-label="Visit model"><div className="panel-heading"><strong>Visit model</strong><Button variant="ghost" className="icon-button" aria-label="Close visit" onClick={close}><X size={18}/></Button></div>
 <div className="section-preview visit-preview"><canvas aria-label="3D model with visit path"/></div>
 <label className="visit-path"><select aria-label="Visit path" value={value.path} onChange={e=>{progressRef.current=0;setProgress(0);change(s=>({...s,path:e.target.value as VisitState['path'],playing:false,seekId:0,lookAlong:e.target.value==='alimentary',restart:s.restart+1}));}}><option value="orbit">Around the model</option><option value="front">Front · head to feet</option><option value="back">Back · head to feet</option>{alimentary&&<option value="alimentary">Alimentary tube · mouth to anus</option>}</select></label>
 <div className="visit-play"><Button onClick={()=>change(s=>({...s,playing:true,direction:1,seekId:0,restart:s.restart+1}))}><Play size={16}/>Forward</Button><Button onClick={()=>change(s=>({...s,playing:true,direction:-1,seekId:0,restart:s.restart+1}))}><Play size={16}/>Reverse</Button><Button aria-label={value.playing?'Pause visit':'Resume visit'} onClick={()=>change(s=>({...s,playing:!s.playing,...(!s.playing&&(s.direction===1?progressRef.current>=1:progressRef.current<=0)?{restart:s.restart+1,seekId:0}:{})}))}>{value.playing?<Pause size={16}/>:<Play size={16}/>}</Button></div>
 <label className="check-label"><input type="checkbox" checked={value.lookAlong} onChange={e=>change(s=>({...s,lookAlong:e.target.checked}))}/>Look along travel direction</label>
 <label>Speed <output>{value.speed.toFixed(2)}×</output><input aria-label="Visit speed" type="range" min=".1" max="4" step=".01" value={value.speed} onChange={e=>change(s=>({...s,speed:Number(e.target.value)}))}/></label>
 <label>Visit position <output>{positionText}</output><input className="visit-position" style={{background:gradient?`linear-gradient(to right,${gradient})`:undefined}} aria-label="Visit position" aria-valuetext={regions.length?`${regions[current].name}, ${Math.round(progress*100)} percent`:undefined} type="range" min="0" max="1" step=".0001" value={progress} onChange={e=>seek(Number(e.target.value))}/></label>
 {!!regions.length&&<ol ref={listRef} className="visit-organs" aria-label="Organs visited">{regions.map((r,i)=><li key={r.name}><button type="button" aria-current={current===i?'step':undefined} onClick={()=>seek(r.start,true)}><span className="visit-bullet" style={{background:r.color}}/>{r.name}</button></li>)}</ol>}
 <div className="visit-footer"><Button className="visit-restart" variant="ghost" onClick={()=>change(s=>({...s,seekId:0,restart:s.restart+1}))}><RotateCcw size={16}/>Restart path</Button><label className="check-label"><input type="checkbox" checked={value.loop} onChange={e=>change(s=>({...s,loop:e.target.checked}))}/>Loop</label></div></section>;
}
