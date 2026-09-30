export interface SceneFrame {canvas:HTMLCanvasElement;labels:SVGSVGElement;transitioning?:boolean}
export interface CaptureSize {width:number;height:number}
export interface CapturedOverlay {image:HTMLImageElement;x:number;y:number;width:number;height:number}
export type SceneFrameListener=(frame:SceneFrame,time:number)=>void;
export type SceneFrameCapture={(size?:CaptureSize):SceneFrame;release:()=>void;invalidate:()=>void;subscribe:(listener:SceneFrameListener,size:CaptureSize)=>()=>void};
export const VIDEO_FPS=60;
export const captureScale=(size:CaptureSize,width:number,height:number)=>Math.max(1,Math.min(size.width/width,size.height/height));
export const videoBitrate=(size:CaptureSize)=>Math.round(Math.min(100_000_000,Math.max(12_000_000,size.width*size.height*VIDEO_FPS*.4)));
let fontCss:Promise<string>|undefined;
const captureIncluded=(element:HTMLElement)=>!element.classList?.contains('anatomy-labels')&&!element.classList?.contains('recording-outline')&&!element.hasAttribute?.('data-scene-canvas')&&!['SCRIPT','STYLE'].includes(element.tagName);
// Copy the properties used by the atlas rather than hundreds of unrelated
// browser defaults per node. Large computed-style snapshots block animation.
const captureStyleProperties=[
 'display','position','top','right','bottom','left','z-index','box-sizing','width','height','min-width','min-height','max-width','max-height',
 'margin-top','margin-right','margin-bottom','margin-left','padding-top','padding-right','padding-bottom','padding-left',
 'flex-direction','flex-wrap','flex-grow','flex-shrink','flex-basis','order','align-items','align-self','align-content','justify-content','justify-items','justify-self','row-gap','column-gap',
 'grid-template-columns','grid-template-rows','grid-auto-columns','grid-auto-rows','grid-auto-flow','grid-column','grid-row',
 'overflow-x','overflow-y','scrollbar-width','scrollbar-color','scrollbar-gutter','visibility','opacity','transform','transform-origin','translate','scale','rotate','isolation',
 'color','background-color','background-image','background-size','background-position','background-repeat','background-clip',
 'border-top','border-right','border-bottom','border-left','border-radius','box-shadow','outline','outline-offset',
 'font-family','font-size','font-weight','font-style','font-stretch','font-variant','font-feature-settings','font-variation-settings','font-kerning','font-synthesis','-webkit-font-smoothing','line-height','letter-spacing',
 'text-align','text-decoration','text-transform','text-overflow','white-space','word-break','overflow-wrap','vertical-align','direction',
 'appearance','accent-color','list-style-type','list-style-position','content','filter','backdrop-filter','clip-path',
 'fill','stroke','stroke-width','stroke-linecap','stroke-linejoin','stroke-dasharray','stroke-dashoffset','fill-rule','paint-order','text-anchor',
];

function imageFromSvg(markup:string){
 return new Promise<HTMLImageElement>((resolve,reject)=>{
  const image=new Image();
  image.onload=()=>resolve(image);
  image.onerror=()=>reject(new Error('The page overlay could not be rendered.'));
  image.src=`data:image/svg+xml;charset=utf-8,${encodeURIComponent(markup)}`;
 });
}

function preserveScroll(source:Element,clone:Element,included=captureIncluded){
 const children=Array.from(source.children).filter(element=>included(element as HTMLElement)),copies=Array.from(clone.children).filter(element=>element.tagName.toUpperCase()!=='STYLE');
 children.forEach((child,index)=>{if(copies[index])preserveScroll(child,copies[index],included);});
 if(!(source instanceof HTMLElement)||!(clone instanceof HTMLElement)||!(source.scrollTop||source.scrollLeft))return;
 clone.style.overflow='hidden';
 for(const child of Array.from(clone.children)){
  if(!(child instanceof HTMLElement))continue;
  const position=child.style.position;
  if(position==='sticky'){const rect=children[copies.indexOf(child)]?.getBoundingClientRect(),container=source.getBoundingClientRect();if(rect){child.style.position='absolute';child.style.top=`${rect.top-container.top}px`;child.style.left=`${rect.left-container.left}px`;}}else child.style.transform=`translate(${-source.scrollLeft}px,${-source.scrollTop}px) ${child.style.transform==='none'?'':child.style.transform}`;
 }
}

async function renderOverlay(root:HTMLElement,width:number,height:number,scale:number,included=captureIncluded,positioned=false){
 const {toSvg,getFontEmbedCSS}=await import('html-to-image');
 fontCss??=getFontEmbedCSS(root.ownerDocument.body).catch(error=>{fontCss=undefined;throw error;});
 const url=await toSvg(root,{
  fontEmbedCSS:await fontCss,
  width,height,
  style:positioned?{position:'relative',top:'0',left:'0',right:'auto',bottom:'auto',transform:'none',margin:'0'}:{background:'transparent'},
  filter:included,
  includeStyleProperties:captureStyleProperties,
 });
 const markup=decodeURIComponent(url.slice(url.indexOf(',')+1)),documentClone=new DOMParser().parseFromString(markup,'image/svg+xml');
 const clone=documentClone.querySelector('foreignObject')?.firstElementChild;
 if(!clone)throw new Error('The page overlay could not be created.');
 preserveScroll(root,clone,included);
 documentClone.documentElement.setAttribute('width',String(Math.round(width*scale)));
 documentClone.documentElement.setAttribute('height',String(Math.round(height*scale)));
 return imageFromSvg(new XMLSerializer().serializeToString(documentClone));
}

export const renderPageOverlay=(studio:HTMLElement,scale=1,separateAdvanced=false)=>renderOverlay(studio.ownerDocument.body,studio.clientWidth,studio.clientHeight,scale,element=>captureIncluded(element)&&(!separateAdvanced||!element.classList?.contains('advanced-panel')));
export async function renderAdvancedOverlay(studio:HTMLElement,scale=1):Promise<CapturedOverlay|null>{
 const panel=studio.querySelector<HTMLElement>('.advanced-panel');if(!panel)return null;
 const rect=panel.getBoundingClientRect(),base=studio.getBoundingClientRect();if(!rect.width||!rect.height)return null;
 return {image:await renderOverlay(panel,rect.width,rect.height,scale,captureIncluded,true),x:rect.left-base.left,y:rect.top-base.top,width:rect.width,height:rect.height};
}

export async function renderSceneLabels(labels:SVGSVGElement,width:number,height:number,scale=1){
 if(!labels.childElementCount)return null;
 const clone=labels.cloneNode(true) as SVGSVGElement;
 const source=[labels,...labels.querySelectorAll('*')],target=[clone,...clone.querySelectorAll('*')];
 source.forEach((element,index)=>{
  const style=getComputedStyle(element);
  for(const property of ['fill','stroke','stroke-width','font-family','font-size','font-weight','opacity','paint-order','stroke-linejoin'])target[index].setAttribute(property,style.getPropertyValue(property));
 });
 clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('viewBox',`0 0 ${width} ${height}`);clone.setAttribute('width',String(Math.round(width*scale)));clone.setAttribute('height',String(Math.round(height*scale)));
 return imageFromSvg(new XMLSerializer().serializeToString(clone));
}

type LabelStyle={fill:string;stroke:string;strokeWidth:number;opacity:number;font:string;paintOrder:string;lineJoin:CanvasLineJoin;textAlign:CanvasTextAlign};
const labelStyles=new WeakMap<SVGSVGElement,{theme:string;styles:Map<string,LabelStyle>}>();

// Atlas labels use these SVG primitives. Painting them directly avoids image
// decoding between a camera frame and its corresponding label positions.
function drawLiveLabels(context:CanvasRenderingContext2D,labels:SVGSVGElement,width:number){
 const theme=`${document.documentElement.dataset.theme!=='light'}:${width<768}`;
 let cached=labelStyles.get(labels);if(!cached||cached.theme!==theme){cached={theme,styles:new Map()};labelStyles.set(labels,cached);}
 const visit=(element:Element,parentKey:string)=>{
  const tag=element.tagName.toLowerCase();if(tag==='title')return;
  const key=`${parentKey}/${tag}.${element.getAttribute('class')??''}.${element.getAttribute('style')??''}`;
  let style=cached.styles.get(key);
  if(!style){const computed=getComputedStyle(element),anchor=element.getAttribute('text-anchor');style={fill:computed.fill,stroke:computed.stroke,strokeWidth:parseFloat(computed.strokeWidth)||0,opacity:Number(computed.opacity),font:`${computed.fontStyle} ${computed.fontWeight} ${computed.fontSize} ${computed.fontFamily}`,paintOrder:computed.paintOrder,lineJoin:computed.strokeLinejoin as CanvasLineJoin,textAlign:anchor==='middle'?'center':anchor==='end'?'right':'left'};cached.styles.set(key,style);}
  context.save();context.globalAlpha*=style.opacity;
  if(tag==='g'){for(const child of Array.from(element.children))visit(child,key);context.restore();return;}
  const value=(name:string)=>Number(element.getAttribute(name)??0);
  context.fillStyle=style.fill;context.strokeStyle=style.stroke;context.lineWidth=style.strokeWidth;context.lineJoin=style.lineJoin;
  const fill=style.fill!=='none',stroke=style.stroke!=='none'&&style.strokeWidth>0;
  if(tag==='text'){
   context.font=style.font;context.textAlign=style.textAlign;context.textBaseline='alphabetic';
   const text=element.textContent??'',x=value('x'),y=value('y'),paintFill=()=>{if(fill)context.fillText(text,x,y);},paintStroke=()=>{if(stroke)context.strokeText(text,x,y);};
   if(style.paintOrder.startsWith('stroke')){paintStroke();paintFill();}else{paintFill();paintStroke();}
  }else{
   context.beginPath();
   if(tag==='line'){context.moveTo(value('x1'),value('y1'));context.lineTo(value('x2'),value('y2'));}
   else if(tag==='circle')context.arc(value('cx'),value('cy'),value('r'),0,Math.PI*2);
   else if(tag==='rect')context.roundRect(value('x'),value('y'),value('width'),value('height'),value('rx'));
   if(fill&&tag!=='line')context.fill();if(stroke)context.stroke();
  }
  context.restore();
 };
 if(getComputedStyle(labels).display!=='none')for(const child of Array.from(labels.children))visit(child,'');
}

export function drawPageFrame(output:HTMLCanvasElement,scene:HTMLCanvasElement,labels:HTMLImageElement|SVGSVGElement|null,overlay:HTMLImageElement|null,sourceWidth:number,sourceHeight:number,advanced:CapturedOverlay|null=null){
 const context=output.getContext('2d');if(!context)throw new Error('The page frame could not be created.');
 const scale=Math.min(output.width/sourceWidth,output.height/sourceHeight),width=sourceWidth*scale,height=sourceHeight*scale,x=(output.width-width)/2,y=(output.height-height)/2;
 context.imageSmoothingEnabled=true;context.imageSmoothingQuality='high';
 context.fillStyle=document.documentElement.dataset.theme!=='light'?'#141b23':'#f2f3f3';context.fillRect(0,0,output.width,output.height);
 context.drawImage(scene,x,y,width,height);
 if(labels instanceof SVGSVGElement){context.save();context.translate(x,y);context.scale(scale,scale);context.beginPath();context.rect(0,0,sourceWidth,sourceHeight);context.clip();drawLiveLabels(context,labels,sourceWidth);context.restore();}
 else if(labels)context.drawImage(labels,x,y,width,height);
 if(overlay)context.drawImage(overlay,x,y,width,height);
 if(advanced)context.drawImage(advanced.image,x+advanced.x*scale,y+advanced.y*scale,advanced.width*scale,advanced.height*scale);
}
