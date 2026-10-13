import {lateralityClass} from './laterality';
import type {LabelPlacement} from './label-layout';

/** Shared selected-item callouts for atlas tissues and sensory overlays. */
export function appendAnatomyLabelCallouts(layer:SVGSVGElement,placements:LabelPlacement[],onShowComponents?:(id:string,components:string[])=>void,onSelectLabel?:(id:string)=>void){
 const ns='http://www.w3.org/2000/svg';
 const leaders=document.createElementNS(ns,'g'),pills=document.createElementNS(ns,'g');
 leaders.setAttribute('class','selected-anatomy-label');pills.setAttribute('class','selected-anatomy-label');layer.appendChild(leaders);layer.appendChild(pills);
 const make=(parent:SVGElement,name:string,attrs:Record<string,string>)=>{const element=document.createElementNS(ns,name);Object.entries(attrs).forEach(([key,value])=>element.setAttribute(key,value));parent.appendChild(element);return element;};
 for(const placed of placements){
  const cls=lateralityClass(placed.laterality??null);
  make(leaders,'line',{class:cls,x1:String(placed.x),y1:String(placed.y),x2:String(placed.leaderX),y2:String(placed.leaderY)});
  make(leaders,'circle',{class:cls,cx:String(placed.x),cy:String(placed.y),r:'3'});
 }
 for(const placed of placements){
  const cls=lateralityClass(placed.laterality??null);
  const item=make(pills,'g',onSelectLabel?{class:'anatomy-label-item',role:'button',tabindex:'0','aria-label':`Select ${placed.fullText??placed.text}`}:{ });
  make(item,'rect',{class:cls,x:String(placed.left),y:String(placed.top),width:String(placed.right-placed.left),height:String(placed.bottom-placed.top),rx:'5'});
  const label=make(item,'text',{class:cls,x:String(placed.left+8),y:String(placed.leaderY+3.5)});label.textContent=placed.label;
  const title=make(item,'title',{});title.textContent=placed.fullText??placed.text;
  if(onSelectLabel){
   const activate=(event:Event)=>{event.preventDefault();event.stopPropagation();onSelectLabel(placed.id);};
   item.addEventListener('pointerdown',event=>event.stopPropagation());item.addEventListener('click',activate);
   item.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' ')activate(event);});
  }
  if(placed.components?.length&&onShowComponents){
   const button=make(pills,'g',{class:'anatomy-label-components',role:'button',tabindex:'0','aria-label':`Show components of ${placed.fullText??placed.text}`});
   make(button,'rect',{x:String(placed.right-24),y:String(placed.top),width:'24',height:String(placed.bottom-placed.top),rx:'4'});
   make(button,'path',{d:`M ${placed.right-17} ${placed.leaderY-2} l 5 5 l 5 -5`,fill:'none'});
   const activate=(event:Event)=>{event.preventDefault();event.stopPropagation();onShowComponents(placed.id,placed.components!);};
   button.addEventListener('pointerdown',event=>event.stopPropagation());button.addEventListener('click',activate);
   button.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' ')activate(event);});
   const hint=make(button,'title',{});hint.textContent='Show components';
  }
 }
}
