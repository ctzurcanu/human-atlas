import * as T from 'three';
import type {AlimentaryVisit} from './alimentary-visit';
import {visitCurve,type VisitState} from './model-visit';

const colors=['#819ead','#e47b76','#b994d6','#738fd2','#d9789b','#c96975','#e59b5e','#d1b15a','#98ba65','#61b292','#4ab4c4','#619dd4','#8587d5','#b082cb','#ce7ca7','#ce946c','#a3a887','#849dac','#789aab'];
export type VisitRegion={name:string;start:number;color:string};
export function visitRegions(bounds:T.Box3,path:VisitState['path'],route?:Pick<AlimentaryVisit,'points'|'segments'>):VisitRegion[]{
 if(path==='orbit')return [0,90,180,270].map((n,i)=>({name:`${n}–${n+90}°`,start:i/4,color:colors[i]}));
 if(path!=='alimentary'){const height=bounds.getSize(new T.Vector3()).y*100;return [{name:'Head',start:0},{name:'Torso',start:.14},{name:'Pelvis',start:.5},{name:'Legs',start:.62},{name:'Feet',start:.94}].map((r,i,a)=>({...r,name:`${r.name} · ${Math.round(r.start*height)}–${Math.round((a[i+1]?.start??1)*height)} cm`,color:colors[i]}));}
 if(!route)return [];
 const curve=visitCurve(bounds,path,route.points) as unknown as T.Curve<T.Vector3>&{waypointProgress:number[]};
 return route.segments.map((s,i)=>({name:s.name,start:curve.waypointProgress[s.startIndex]??0,color:colors[i%colors.length]}));
}
export function visitMapSegments(curve:T.Curve<T.Vector3>,regions:VisitRegion[]){
 // Include every slider boundary so a line segment never straddles two colors.
 const samples=[...new Set([...Array.from({length:513},(_,i)=>i/512),...regions.map(r=>r.start)])].sort((a,b)=>a-b);
 const positions:number[]=[],colors:number[]=[],bounds=new T.Box3();
 for(let i=0;i<samples.length-1;i++){
  const middle=(samples[i]+samples[i+1])/2,region=regions.filter(r=>r.start<=middle).at(-1),color=new T.Color(region?.color??'#819ead');
  for(const t of [samples[i],samples[i+1]]){const point=curve.getPointAt(t);bounds.expandByPoint(point);positions.push(...point.toArray());colors.push(...color.toArray());}
 }
 return {positions,colors,bounds};
}
