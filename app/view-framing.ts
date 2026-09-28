import {Box3,MathUtils,PerspectiveCamera,Quaternion,Vector3} from 'three';

export type ViewArea={left:number;right:number;top:number;bottom:number};
export type ViewRect=ViewArea&{kind:'top'|'bottom'|'panel'};
const clamp=(value:number)=>MathUtils.clamp(value,0,1);
const span=(a:number,b:number)=>Math.max(.000001,b-a);
export const immersiveArea:ViewArea={left:.04,right:.96,top:.04,bottom:.86};

/** Fractions of the local viewport; panel dimensions never travel with a view. */
export function usableViewArea(rects:ViewRect[]):ViewArea{
 let left=.02,right=.98,top=.04,bottom=.94;
 for(const rect of rects){
  if(rect.kind==='top')top=Math.max(top,rect.bottom+.015);
  else if(rect.kind==='bottom')bottom=Math.min(bottom,rect.top-.015);
  else if(rect.right-rect.left<.5&&rect.bottom-rect.top>.35){
   if((rect.left+rect.right)/2<.5)left=Math.max(left,rect.right+.02);
   else right=Math.min(right,rect.left-.02);
  }else if(rect.top>.15)bottom=Math.min(bottom,rect.top-.02);
 }
 left=clamp(left);right=clamp(right);top=clamp(top);bottom=clamp(bottom);
 if(right-left<.08){const center=clamp((left+right)/2);left=Math.max(0,center-.04);right=Math.min(1,center+.04);}
 if(bottom-top<.08){const center=clamp((top+bottom)/2);top=Math.max(0,center-.04);bottom=Math.min(1,center+.04);}
 return {left,right,top,bottom};
}
export function preferredAnchor(area:ViewArea):[number,number]{
 return [.5,MathUtils.clamp((Math.min(.46,(area.top+area.bottom)/2)-area.top)/span(area.top,area.bottom),0,1)];
}
export function screenAnchor(area:ViewArea,anchor:readonly number[]){
 // Even a saved screen center is kept above the physical viewport midpoint.
 return [area.left+clamp(anchor[0])*span(area.left,area.right),Math.min(.46,area.top+clamp(anchor[1])*span(area.top,area.bottom))];
}
export function setFrameOffset(camera:PerspectiveCamera,area:ViewArea,anchor:readonly number[]=preferredAnchor(area)){
 const [x,y]=screenAnchor(area,anchor);
 // Three's projection API wants pixels, but ratios are enough (height = 1).
 camera.setViewOffset(camera.aspect,1,(.5-x)*camera.aspect,.5-y,camera.aspect,1);
}
export function frameDistance(bounds:Box3,direction:Vector3,up:Vector3,aspect:number,area:ViewArea,fov:number){
 const size=bounds.getSize(new Vector3()),right=new Vector3().crossVectors(up,direction).normalize();
 if(right.lengthSq()<.001)right.set(1,0,0);
 const vertical=new Vector3().crossVectors(direction,right).normalize();
 const extent=(axis:Vector3)=>Math.abs(axis.x)*size.x+Math.abs(axis.y)*size.y+Math.abs(axis.z)*size.z;
 const tan=Math.tan(MathUtils.degToRad(fov/2));
 const [x,y]=screenAnchor(area,preferredAnchor(area)),availableHeight=Math.max(.02,2*Math.min(y-area.top,area.bottom-y)),availableWidth=Math.max(.02,2*Math.min(x-area.left,area.right-x));
 return Math.max(extent(vertical)/(2*tan*availableHeight),extent(right)/(2*tan*Math.max(.001,aspect)*availableWidth))+extent(direction)/2;
}
// Frame v1: unit direction [0..1]×3, model-relative target×3, usable-area
// anchor [0..1]×2, distance fraction [0..1], unit up [0..1]×3.
// Targets can leave [0..1] when the user pans beyond the model's bounds.
export function captureFrame(camera:PerspectiveCamera,target:Vector3,bounds:Box3,area:ViewArea):number[]{
 const size=bounds.getSize(new Vector3()),direction=camera.position.clone().sub(target),distance=direction.length();direction.normalize();
 const fit=frameDistance(bounds,direction,camera.up,camera.aspect,area,camera.fov);
 const x=.5-(camera.view?.enabled?camera.view.offsetX/camera.view.fullWidth:0),y=.5-(camera.view?.enabled?camera.view.offsetY/camera.view.fullHeight:0);
 const up=camera.up.clone().normalize();
 return [...direction.toArray().map(value=>(value+1)/2),...target.toArray().map((value,i)=>(value-bounds.min.getComponent(i))/Math.max(size.getComponent(i),.000001)),clamp((x-area.left)/span(area.left,area.right)),clamp((y-area.top)/span(area.top,area.bottom)),MathUtils.clamp(distance/(distance+fit),.000001,.999999),...up.toArray().map(value=>(value+1)/2)];
}
export function decodeFrame(values:number[],bounds:Box3,aspect:number,area:ViewArea,fov:number){
 if(values.length!==12){
  // Old URLs lack their source viewport dimensions. Use a canonical source
  // projection, then adapt its framing instead of copying its screen offset.
  const legacy=new PerspectiveCamera(34,16/9),target=new Vector3().fromArray(values,3),sourceArea=usableViewArea([]);
  legacy.position.fromArray(values);if(Math.abs(legacy.position.clone().sub(target).normalize().y)>.95)legacy.up.set(0,0,-1);legacy.lookAt(target);setFrameOffset(legacy,sourceArea);
  return decodeFrame(captureFrame(legacy,target,bounds,sourceArea),bounds,aspect,area,fov);
 }
 const direction=new Vector3().fromArray(values).multiplyScalar(2).subScalar(1).normalize(),up=new Vector3().fromArray(values,9).multiplyScalar(2).subScalar(1).normalize();
 const size=bounds.getSize(new Vector3()),target=new Vector3().fromArray(values,3).multiply(size).add(bounds.min);
 const zoom=MathUtils.clamp(values[8],.000001,.999999),distance=frameDistance(bounds,direction,up,aspect,area,fov)*zoom/(1-zoom);
 return {position:target.clone().addScaledVector(direction,distance),target,up,anchor:values.slice(6,8)};
}
export function applyFrame(camera:PerspectiveCamera,target:Vector3,values:number[],bounds:Box3,area:ViewArea){
 const frame=decodeFrame(values,bounds,camera.aspect,area,camera.fov);
 camera.position.copy(frame.position);target.copy(frame.target);camera.up.copy(frame.up);camera.lookAt(target);setFrameOffset(camera,area,frame.anchor);camera.updateMatrixWorld(true);
}
export function interpolateFrame(camera:PerspectiveCamera,target:Vector3,from:number[],to:number[],bounds:Box3,area:ViewArea,mix:number){
 const a=decodeFrame(from,bounds,camera.aspect,area,camera.fov),b=decodeFrame(to,bounds,camera.aspect,area,camera.fov);
 const start=a.position.sub(a.target),end=b.position.sub(b.target),startDistance=Math.max(.000001,start.length()),endDistance=Math.max(.000001,end.length());
 const rotation=new Quaternion().setFromUnitVectors(start.normalize(),end.normalize());
 const direction=start.applyQuaternion(new Quaternion().slerp(rotation,mix)),distance=Math.exp(MathUtils.lerp(Math.log(startDistance),Math.log(endDistance),mix));
 target.copy(a.target.lerp(b.target,mix));camera.position.copy(target).addScaledVector(direction,distance);camera.up.copy(a.up.lerp(b.up,mix).normalize());camera.lookAt(target);
 setFrameOffset(camera,area,a.anchor.map((value,i)=>MathUtils.lerp(value,b.anchor[i],mix)));camera.updateMatrixWorld(true);
}
