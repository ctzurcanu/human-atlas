import type {Laterality} from './laterality';

export interface LabelAnchor {id:string;x:number;y:number;text:string;fullText?:string;laterality?:Laterality}
export interface LabelPlacement extends LabelAnchor {side:'left'|'right';label:string;left:number;right:number;top:number;bottom:number;leaderX:number;leaderY:number}
export interface LabelArea {left:number;right:number;top:number;bottom:number}
export interface LabelPanel extends LabelArea {kind:'top'|'bottom'|'panel'}
export const LABEL_PANEL_GAP=5;

export function labelAreaForPanels(width:number,height:number,panels:LabelPanel[]):LabelArea{
 const area={left:LABEL_PANEL_GAP,right:width-LABEL_PANEL_GAP,top:LABEL_PANEL_GAP,bottom:height-LABEL_PANEL_GAP};
 for(const panel of panels){
  if(panel.right<=0||panel.left>=width||panel.bottom<=0||panel.top>=height)continue;
  if(panel.kind==='top')area.top=Math.max(area.top,panel.bottom+LABEL_PANEL_GAP);
  else if(panel.kind==='bottom')area.bottom=Math.min(area.bottom,panel.top-LABEL_PANEL_GAP);
  else if(panel.right-panel.left<width/2){
   if((panel.left+panel.right)/2<width/2)area.left=Math.max(area.left,panel.right+LABEL_PANEL_GAP);
   else area.right=Math.min(area.right,panel.left-LABEL_PANEL_GAP);
  }else if(panel.top>=height-panel.bottom)area.bottom=Math.min(area.bottom,panel.top-LABEL_PANEL_GAP);
  else area.top=Math.max(area.top,panel.bottom+LABEL_PANEL_GAP);
 }
 return area;
}

const clamp=(value:number,min:number,max:number)=>Math.max(min,Math.min(max,value));

// Place selected labels in two columns outside the projected selection. Each
// column reserves one row per label, so text and leaders remain legible as the
// model rotates or the selected relationship set grows.
export function layoutAnatomyLabels(anchors:LabelAnchor[],area:LabelArea,measure:(text:string)=>number):LabelPlacement[]{
 const available=Math.max(0,area.right-area.left),height=Math.max(0,area.bottom-area.top);
 if(!anchors.length||available<130||height<28)return [];
 const compact=available<520;
 const rows=Math.ceil(anchors.length/2),defaultRowHeight=compact?21:24;
 // Dense relationship sets (such as a hip bone) need tighter rows rather than
 // disappearing at a fixed selection count. Keep a readable minimum spacing.
 const gap=Math.min(compact?25:30,Math.max(18,(height-defaultRowHeight)/Math.max(1,rows-1)));
 const rowHeight=Math.min(defaultRowHeight,gap-3);
 const maxPerSide=Math.max(1,Math.floor((height-rowHeight)/gap)+1);
 const maxWidth=Math.max(70,Math.min(compact?135:230,(available-50)/2));
 const midpoint=(area.left+area.right)/2;
 const selected=anchors.filter(anchor=>Number.isFinite(anchor.x)&&Number.isFinite(anchor.y)).slice(0,maxPerSide*2);
 const left:LabelAnchor[]=[],right:LabelAnchor[]=[];
 for(const anchor of selected){
  const preferred=anchor.x<midpoint?left:right,alternate=preferred===left?right:left;
  if(preferred.length<maxPerSide)preferred.push(anchor);
  else if(alternate.length<maxPerSide)alternate.push(anchor);
 }
 const fit=(text:string)=>{
  if(measure(text)<=maxWidth-16)return text;
  let lo=0,hi=text.length;
  while(lo<hi){const mid=Math.ceil((lo+hi)/2);if(measure(`${text.slice(0,mid)}…`)<=maxWidth-16)lo=mid;else hi=mid-1;}
  return `${text.slice(0,lo)}…`;
 };
 const layout=(items:LabelAnchor[],side:'left'|'right')=>{
  items.sort((a,b)=>a.y-b.y||a.x-b.x||a.id.localeCompare(b.id));
  const output:LabelPlacement[]=[];
  const rowTop=area.top+rowHeight/2,rowBottom=area.bottom-rowHeight/2;
  for(let index=0;index<items.length;index++){
   const anchor=items[index],label=fit(anchor.text),width=Math.min(maxWidth,measure(label)+16);
   const minY=Math.max(rowTop+index*gap,index?output[index-1].leaderY+gap:rowTop),maxY=rowBottom-(items.length-1-index)*gap;
   const cy=clamp(anchor.y,minY,maxY);
   const leftX=side==='left'?area.left:area.right-width;
   output.push({...anchor,side,label,left:leftX,right:leftX+width,top:cy-rowHeight/2,bottom:cy+rowHeight/2,leaderX:side==='left'?leftX+width:leftX,leaderY:cy});
  }
  return output;
 };
 return [...layout(left,'left'),...layout(right,'right')];
}
