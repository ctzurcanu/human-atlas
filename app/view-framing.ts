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
/** Pull back only as far as needed; retain the current orientation and pan. */
export function zoomOutToFit(camera:PerspectiveCamera,bounds:Box3,area:ViewArea){
 if(bounds.isEmpty())return false;
 const paddingX=(area.right-area.left)*.02,paddingY=(area.bottom-area.top)*.02;
 const left=area.left+paddingX,right=area.right-paddingX,top=area.top+paddingY,bottom=area.bottom-paddingY;
 camera.updateMatrixWorld(true);
 const points=Array.from({length:8},(_,corner)=>new Vector3(bounds[(corner&1)?'max':'min'].x,bounds[(corner&2)?'max':'min'].y,bounds[(corner&4)?'max':'min'].z).applyMatrix4(camera.matrixWorldInverse));
 const fits=()=>points.every(point=>{const projected=point.clone().applyMatrix4(camera.projectionMatrix),x=(projected.x+1)/2,y=(1-projected.y)/2;return -point.z>=camera.near&&-point.z<=camera.far&&x>=left-1e-10&&x<=right+1e-10&&y>=top-1e-10&&y<=bottom+1e-10;});
 if(fits())return false;
 let projection=camera.projectionMatrix.elements,anchorX=(1-projection[8])/2,anchorY=(1+projection[9])/2;
 // A panel opened since the last fit can cover the projection center. Re-anchor
 // only during this explicit fit, otherwise pulling back cannot uncover it.
 if(anchorX<=left||anchorX>=right||anchorY<=top||anchorY>=bottom){setFrameOffset(camera,area);projection=camera.projectionMatrix.elements;anchorX=(1-projection[8])/2;anchorY=(1+projection[9])/2;}
 let pullback=0,maxDepth=0;
 for(const point of points){
  const horizontal=point.x*projection[0]/2,vertical=-point.y*projection[5]/2,depth=-point.z;
  const needed=Math.max(camera.near*1.01,Math.abs(horizontal)/(horizontal<0?anchorX-left:right-anchorX),Math.abs(vertical)/(vertical<0?anchorY-top:bottom-anchorY));
  pullback=Math.max(pullback,needed-depth);maxDepth=Math.max(maxDepth,depth);
 }
 if(pullback>0)camera.position.addScaledVector(new Vector3(0,0,1).applyQuaternion(camera.quaternion),pullback);
 camera.far=Math.max(camera.far,(maxDepth+pullback)*1.05);camera.updateProjectionMatrix();camera.updateMatrixWorld(true);
 return true;
}
/** Enlarge shrinking content to more than half of each usable dimension, up to its fit limit. */
export function zoomInToFill(camera:PerspectiveCamera,bounds:Box3,area:ViewArea,maxAdvance=Infinity){
 if(bounds.isEmpty()||maxAdvance<=0)return false;
 camera.updateMatrixWorld(true);
 const projection=camera.projectionMatrix.elements,anchorX=(1-projection[8])/2,anchorY=(1+projection[9])/2;
 const width=area.right-area.left,height=area.bottom-area.top,left=area.left+width*.02,right=area.right-width*.02,top=area.top+height*.02,bottom=area.bottom-height*.02;
 if(anchorX<=left||anchorX>=right||anchorY<=top||anchorY>=bottom)return false;
 const points=Array.from({length:8},(_,corner)=>{
  const point=new Vector3(bounds[(corner&1)?'max':'min'].x,bounds[(corner&2)?'max':'min'].y,bounds[(corner&4)?'max':'min'].z).applyMatrix4(camera.matrixWorldInverse);
  return {x:point.x*projection[0]/2,y:-point.y*projection[5]/2,depth:-point.z};
 });
 const fills=(advance:number)=>{
  let minX=Infinity,maxX=-Infinity,minY=Infinity,maxY=-Infinity;
  for(const point of points){const depth=point.depth-advance,x=point.x/depth,y=point.y/depth;minX=Math.min(minX,x);maxX=Math.max(maxX,x);minY=Math.min(minY,y);maxY=Math.max(maxY,y);}
  return maxX-minX>=width*.52-1e-10&&maxY-minY>=height*.52-1e-10;
 };
 if(fills(0))return false;
 let advance=maxAdvance;
 for(const point of points){
  const needed=Math.max(camera.near*1.01,Math.abs(point.x)/(point.x<0?anchorX-left:right-anchorX),Math.abs(point.y)/(point.y<0?anchorY-top:bottom-anchorY));
  advance=Math.min(advance,point.depth-needed);
 }
 if(advance<=1e-10)return false;
 if(fills(advance)){
  let low=0,high=advance;
  for(let iteration=0;iteration<40;iteration++){const middle=(low+high)/2;if(fills(middle))high=middle;else low=middle;}
  advance=high;
 }
 camera.position.addScaledVector(new Vector3(0,0,1).applyQuaternion(camera.quaternion),-advance);camera.updateMatrixWorld(true);
 return true;
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
