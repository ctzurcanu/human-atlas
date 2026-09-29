export interface SceneFrame {canvas:HTMLCanvasElement;labels:SVGSVGElement}
export interface CaptureSize {width:number;height:number}
export type SceneFrameCapture={(size?:CaptureSize):SceneFrame;release:()=>void};
export const captureScale=(size:CaptureSize,width:number,height:number)=>Math.max(1,Math.min(size.width/width,size.height/height));
export const videoBitrate=(size:CaptureSize)=>Math.round(Math.min(80_000_000,Math.max(8_000_000,size.width*size.height*30*.4)));
let fontCss:Promise<string>|undefined;
const captureIncluded=(element:HTMLElement)=>!element.classList?.contains('anatomy-labels')&&!element.classList?.contains('recording-outline')&&!element.hasAttribute?.('data-scene-canvas')&&!['SCRIPT','STYLE'].includes(element.tagName);

function imageFromSvg(markup:string){
 return new Promise<HTMLImageElement>((resolve,reject)=>{
  const image=new Image();
  image.onload=()=>resolve(image);
  image.onerror=()=>reject(new Error('The page overlay could not be rendered.'));
  image.src=`data:image/svg+xml;charset=utf-8,${encodeURIComponent(markup)}`;
 });
}

function preserveScroll(source:Element,clone:Element){
 const children=Array.from(source.children).filter(element=>captureIncluded(element as HTMLElement)),copies=Array.from(clone.children).filter(element=>element.tagName.toUpperCase()!=='STYLE');
 children.forEach((child,index)=>{if(copies[index])preserveScroll(child,copies[index]);});
 if(!(source instanceof HTMLElement)||!(clone instanceof HTMLElement)||!(source.scrollTop||source.scrollLeft))return;
 clone.style.overflow='hidden';
 for(const child of Array.from(clone.children)){
  if(!(child instanceof HTMLElement))continue;
  const position=child.style.position;
  if(position==='sticky'){const rect=children[copies.indexOf(child)]?.getBoundingClientRect(),container=source.getBoundingClientRect();if(rect){child.style.position='absolute';child.style.top=`${rect.top-container.top}px`;child.style.left=`${rect.left-container.left}px`;}}else child.style.transform=`translate(${-source.scrollLeft}px,${-source.scrollTop}px) ${child.style.transform==='none'?'':child.style.transform}`;
 }
}

export async function renderPageOverlay(studio:HTMLElement,scale=1){
 const {toSvg,getFontEmbedCSS}=await import('html-to-image'),root=studio.ownerDocument.body;
 fontCss??=getFontEmbedCSS(root).catch(error=>{fontCss=undefined;throw error;});
 const url=await toSvg(root,{
  fontEmbedCSS:await fontCss,
  width:studio.clientWidth,
  height:studio.clientHeight,
  style:{background:'transparent'},
  filter:captureIncluded,
 });
 const markup=decodeURIComponent(url.slice(url.indexOf(',')+1)),documentClone=new DOMParser().parseFromString(markup,'image/svg+xml');
 const clone=documentClone.querySelector('foreignObject')?.firstElementChild;
 if(!clone)throw new Error('The page overlay could not be created.');
 preserveScroll(root,clone);
 documentClone.documentElement.setAttribute('width',String(Math.round(studio.clientWidth*scale)));
 documentClone.documentElement.setAttribute('height',String(Math.round(studio.clientHeight*scale)));
 return imageFromSvg(new XMLSerializer().serializeToString(documentClone));
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

export function drawPageFrame(output:HTMLCanvasElement,scene:HTMLCanvasElement,labels:HTMLImageElement|null,overlay:HTMLImageElement|null,sourceWidth:number,sourceHeight:number){
 const context=output.getContext('2d');if(!context)throw new Error('The page frame could not be created.');
 const scale=Math.min(output.width/sourceWidth,output.height/sourceHeight),width=sourceWidth*scale,height=sourceHeight*scale,x=(output.width-width)/2,y=(output.height-height)/2;
 context.imageSmoothingEnabled=true;context.imageSmoothingQuality='high';
 context.fillStyle=matchMedia('(prefers-color-scheme: dark)').matches?'#141b23':'#f2f3f3';context.fillRect(0,0,output.width,output.height);
 context.drawImage(scene,x,y,width,height);
 if(labels)context.drawImage(labels,x,y,width,height);
 if(overlay)context.drawImage(overlay,x,y,width,height);
}
