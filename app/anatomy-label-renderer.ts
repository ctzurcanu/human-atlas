import {lateralityClass} from './laterality';
import type {LabelPlacement} from './label-layout';

/** Shared selected-item callouts for atlas tissues and sensory overlays. */
export function appendAnatomyLabelCallouts(layer:SVGSVGElement,placements:LabelPlacement[]){
 const ns='http://www.w3.org/2000/svg';
 const leaders=document.createElementNS(ns,'g'),pills=document.createElementNS(ns,'g');
 leaders.setAttribute('class','selected-anatomy-label');pills.setAttribute('class','selected-anatomy-label');layer.appendChild(leaders);layer.appendChild(pills);
 const make=(parent:SVGGElement,name:string,attrs:Record<string,string>)=>{const element=document.createElementNS(ns,name);Object.entries(attrs).forEach(([key,value])=>element.setAttribute(key,value));parent.appendChild(element);return element;};
 for(const placed of placements){
  const cls=lateralityClass(placed.laterality??null);
  make(leaders,'line',{class:cls,x1:String(placed.x),y1:String(placed.y),x2:String(placed.leaderX),y2:String(placed.leaderY)});
  make(leaders,'circle',{class:cls,cx:String(placed.x),cy:String(placed.y),r:'3'});
 }
 for(const placed of placements){
  const cls=lateralityClass(placed.laterality??null);
  make(pills,'rect',{class:cls,x:String(placed.left),y:String(placed.top),width:String(placed.right-placed.left),height:String(placed.bottom-placed.top),rx:'5'});
  const label=make(pills,'text',{class:cls,x:String(placed.left+8),y:String(placed.leaderY+3.5)});label.textContent=placed.label;
  const title=make(pills,'title',{});title.textContent=placed.fullText??placed.text;
 }
}
