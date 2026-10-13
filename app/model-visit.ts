import * as T from 'three';
export type VisitState={path:'orbit'|'front'|'back'|'alimentary';playing:boolean;direction:1|-1;speed:number;restart:number;lookAlong:boolean;seek:number;seekId:number;loop:boolean};
export const initialVisit:VisitState={path:'orbit',playing:false,direction:1,speed:1,restart:0,lookAlong:false,seek:0,seekId:0,loop:false};
// Arc-length sampling keeps travel speed independent of control-point spacing.
export function visitCurve(bounds:T.Box3,path:VisitState['path'],waypoints?:[number,number,number][]){
 if(path==='alimentary'&&waypoints?.length){
  const points:T.Vector3[]=[],mapping:number[]=[];for(const p of waypoints){const v=new T.Vector3().fromArray(p);if(!points.length||points.at(-1)!.distanceToSquared(v)>1e-16)points.push(v);mapping.push(points.length-1);}
  // Smooth surface-centroid jitter over a physical 12 mm window, keeping endpoints.
  const distance=[0];for(let i=1;i<points.length;i++)distance.push(distance[i-1]+points[i].distanceTo(points[i-1]));
  const smooth=points.map((p,i)=>{if(i===0||i===points.length-1||i<=1||i>=points.length-2)return p;const v=new T.Vector3();let weight=0;for(let j=Math.max(0,i-8);j<=Math.min(points.length-1,i+8);j++){const w=Math.exp(-Math.pow((distance[j]-distance[i])/.012,2));v.addScaledVector(points[j],w);weight+=w;}return v.divideScalar(weight);});
  const curve=new T.CatmullRomCurve3(smooth,false,'centripetal') as T.CatmullRomCurve3&{waypointProgress:number[]};
  curve.arcLengthDivisions=Math.max(200,points.length*16);
  const lengths=curve.getLengths(),total=lengths.at(-1)!;
  curve.waypointProgress=mapping.map(i=>{const t=i/(points.length-1)*curve.arcLengthDivisions,k=Math.floor(t),f=t-k;return (lengths[k]*(1-f)+(lengths[Math.min(k+1,lengths.length-1)]*f))/total;});return curve;
 }
 const c=bounds.getCenter(new T.Vector3()),s=bounds.getSize(new T.Vector3()),radius=Math.max(s.x,s.z,s.y*.4,.001)*1.35;
 const points=path==='orbit'?Array.from({length:16},(_,i)=>{const a=i/16*Math.PI*2;return new T.Vector3(c.x+Math.sin(a)*radius,c.y+0,c.z+Math.cos(a)*radius);}):Array.from({length:7},(_,i)=>new T.Vector3(c.x,bounds.max.y-i/6*s.y,c.z+(path==='front'?1:-1)*Math.max(s.x*1.5,s.z*2,s.y*.18)));
 return new T.CatmullRomCurve3(points,path==='orbit','centripetal');
}
export function dorsalFrames(curve:T.Curve<T.Vector3>){
 const frames:T.Quaternion[]=[],up=new T.Vector3(0,0,1);let previous=new T.Quaternion();
 for(let i=0;i<=2048;i++){
  const h=.004,t=i/2048,tangent=curve.getPointAt(Math.min(1,t+h)).sub(curve.getPointAt(Math.max(0,t-h))).normalize();
  const dorsal=up.clone().addScaledVector(tangent,-up.dot(tangent));
  // At an axial pole retain the preceding frame's up instead of introducing a roll flip.
  if(dorsal.lengthSq()<.04)dorsal.set(0,1,0).applyQuaternion(previous).addScaledVector(tangent,-new T.Vector3(0,1,0).applyQuaternion(previous).dot(tangent));
  dorsal.normalize();const q=new T.Quaternion().setFromRotationMatrix(new T.Matrix4().lookAt(new T.Vector3(),tangent,dorsal));
  previous=q;frames.push(q);
 }
 return frames;
}
export class ModelVisit{
 progress=0; yaw=0;pitch=0; private command=-1;private seekId=-1;private path='';private curve?:T.Curve<T.Vector3>;private frames:T.Quaternion[]=[];
 update(camera:T.PerspectiveCamera,target:T.Vector3,bounds:T.Box3,state:VisitState,dt:number,waypoints?:[number,number,number][]){
  const starting=state.restart!==this.command||state.path!==this.path;
  if(starting){this.command=state.restart;this.path=state.path;this.curve=visitCurve(bounds,state.path,waypoints);this.progress=state.direction===1?0:1;this.yaw=0;this.pitch=0;this.frames=state.path==='alimentary'?dorsalFrames(this.curve):[];}
  const seeking=state.seekId>0&&state.seekId!==this.seekId;this.seekId=state.seekId;if(seeking&&state.seekId>0)this.progress=T.MathUtils.clamp(state.seek,0,1);
  if((!state.playing&&!seeking)||!this.curve||!state.playing&&state.seekId===0)return false;
  const step=(state.playing&&!seeking?dt:0)*state.speed/(state.path==='alimentary'?120:30)*state.direction,next=this.progress+step;
  this.progress=state.loop&&(next<0||next>1)?((next%1)+1)%1:T.MathUtils.clamp(next,0,1);
  camera.position.copy(this.curve.getPointAt(this.progress));
  const h=state.path==='alimentary'?Math.min(.012,.025/Math.max(this.curve.getLength(),.01)):.002;
  const tangent=state.lookAlong&&state.path!=='orbit'?this.curve.getPointAt(Math.min(1,this.progress+h)).sub(this.curve.getPointAt(Math.max(0,this.progress-h))).multiplyScalar(state.direction).normalize():(state.path==='front'||state.path==='back'?new T.Vector3(bounds.getCenter(new T.Vector3()).x,camera.position.y,bounds.getCenter(new T.Vector3()).z):bounds.getCenter(new T.Vector3())).sub(camera.position).normalize();
  let up=new T.Vector3(0,1,0);
  if(state.path==='alimentary'&&state.lookAlong&&this.frames.length){const t=this.progress*(this.frames.length-1),i=Math.floor(t),q=this.frames[i].clone().slerp(this.frames[Math.min(i+1,this.frames.length-1)],t-i);up.set(0,1,0).applyQuaternion(q);}
  const base=new T.Quaternion().setFromRotationMatrix(new T.Matrix4().lookAt(new T.Vector3(),tangent,up));
  if(state.path==='alimentary'&&state.lookAlong&&state.direction===1&&waypoints?.length){
   const curve=this.curve as T.Curve<T.Vector3>&{waypointProgress:number[]},exit=curve.waypointProgress[waypoints.length-2];
   if(this.progress>exit){const t=T.MathUtils.smoothstep(this.progress,exit,1),anus=new T.Vector3().fromArray(waypoints[waypoints.length-2]),toward=anus.sub(camera.position).normalize();
    const facing=new T.Quaternion().setFromRotationMatrix(new T.Matrix4().lookAt(new T.Vector3(),toward,new T.Vector3(0,1,0)));base.slerp(facing,t);}
  }
  if(state.path==='alimentary'&&state.lookAlong&&waypoints?.length){const curve=this.curve as T.Curve<T.Vector3>&{waypointProgress:number[]},mouth=curve.waypointProgress[1];if(this.progress<mouth&&this.frames.length){const entry=this.frames[Math.min(2048,Math.ceil(mouth*2048))];base.slerp(entry,T.MathUtils.smoothstep(this.progress,0,mouth));}}
  base.multiply(new T.Quaternion().setFromEuler(new T.Euler(this.pitch,this.yaw,0,'YXZ')));
  if(dt>0&&!seeking&&!starting)camera.quaternion.rotateTowards(base,dt*state.speed*1.25);else camera.quaternion.copy(base);camera.up.set(0,1,0).applyQuaternion(camera.quaternion);
  target.copy(camera.position).addScaledVector(new T.Vector3(0,0,-1).applyQuaternion(camera.quaternion),Math.max(bounds.getSize(new T.Vector3()).length()*.1,.01));
  camera.updateMatrixWorld();return !state.loop&&(state.direction===1?next>=1:next<=0);
 }
 look(dx:number,dy:number){this.yaw-=dx*.005;this.pitch=T.MathUtils.clamp(this.pitch-dy*.005,-Math.PI*.48,Math.PI*.48);}
}

// A visit owns the camera only during an explicit play/seek request. Keep the
// exact requesting pose separately from the route's camera and its fitted view.
export class VisitJourney {
 readonly travel=new ModelVisit();phase:'idle'|'approach'|'travel'|'return'|'done'='idle';
 private command=-1;private seek=-1;private elapsed=0;private origin?:T.PerspectiveCamera;private originTarget=new T.Vector3();private from?:T.PerspectiveCamera;private fromTarget=new T.Vector3();private destination?:T.PerspectiveCamera;private destinationTarget=new T.Vector3();
 release(){if(this.phase==='done')this.phase='idle';}
 get progress(){return this.travel.progress;}
 look(dx:number,dy:number){this.travel.look(dx,dy);}
 update(camera:T.PerspectiveCamera,target:T.Vector3,bounds:T.Box3,state:VisitState,dt:number,points:[number,number,number][]|undefined,fit:(camera:T.PerspectiveCamera,target:T.Vector3)=>void=()=>{}) {
  const seeking=state.seekId>0&&state.seekId!==this.seek;this.seek=state.seekId;
  if(seeking){this.phase='idle';this.travel.update(camera,target,bounds,state,0,points);fit(camera,target);return false;}
  if(!state.playing)return false;
  if(this.command!==state.restart||this.phase==='idle'){
   this.command=state.restart;this.origin=camera.clone();this.originTarget.copy(target);this.from=camera.clone();this.fromTarget.copy(target);
   this.destination=camera.clone();this.destinationTarget.copy(target);this.travel.update(this.destination,this.destinationTarget,bounds,state,0,points);fit(this.destination,this.destinationTarget);
   this.phase='approach';this.elapsed=0;return false;
  }
  if(this.phase==='done')return true;
  if(this.phase==='approach'||this.phase==='return'){
   this.elapsed+=Math.max(0,dt);const t=Math.min(1,this.elapsed/2.4),mix=t*t*t*(t*(t*6-15)+10),to=this.phase==='return'?this.origin!:this.destination!,toTarget=this.phase==='return'?this.originTarget:this.destinationTarget;
   camera.position.lerpVectors(this.from!.position,to.position,mix);camera.quaternion.copy(this.from!.quaternion).slerp(to.quaternion,mix);camera.up.copy(this.from!.up).lerp(to.up,mix).normalize();target.lerpVectors(this.fromTarget,toTarget,mix);camera.zoom=T.MathUtils.lerp(this.from!.zoom,to.zoom,mix);
   const a=this.from!.view,b=to.view;if(a?.enabled||b?.enabled){const av=a?.enabled?a:{fullWidth:camera.aspect,fullHeight:1,offsetX:0,offsetY:0,width:camera.aspect,height:1},bv=b?.enabled?b:{fullWidth:camera.aspect,fullHeight:1,offsetX:0,offsetY:0,width:camera.aspect,height:1};camera.setViewOffset(T.MathUtils.lerp(av.fullWidth,bv.fullWidth,mix),T.MathUtils.lerp(av.fullHeight,bv.fullHeight,mix),T.MathUtils.lerp(av.offsetX,bv.offsetX,mix),T.MathUtils.lerp(av.offsetY,bv.offsetY,mix),T.MathUtils.lerp(av.width,bv.width,mix),T.MathUtils.lerp(av.height,bv.height,mix));}
   camera.updateProjectionMatrix();camera.updateMatrixWorld();
   if(t===1){camera.copy(to);target.copy(toTarget);if(this.phase==='return'){this.phase='done';return true;}this.phase='travel';}
   return false;
  }
  const ended=this.travel.update(camera,target,bounds,state,dt,points);fit(camera,target);
  if(ended){this.from=camera.clone();this.fromTarget.copy(target);this.phase='return';this.elapsed=0;}
  return false;
 }
}
