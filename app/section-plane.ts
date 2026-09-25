export type SectionAxis='axial'|'sagittal'|'coronal'|'oblique';
export type SectionState={enabled:boolean;axis:SectionAxis;position:number;flip:boolean;azimuth?:number;elevation?:number};
export type Point3=[number,number,number];
export type Bounds3={min:Point3;max:Point3};
const clamp=(value:number,min:number,max:number)=>Math.max(min,Math.min(max,value));
export function sectionNormal(section:Pick<SectionState,'axis'|'azimuth'|'elevation'>):Point3{
 if(section.axis==='axial')return [0,-1,0];
 if(section.axis==='sagittal')return [1,0,0];
 if(section.axis==='coronal')return [0,0,-1];
 const azimuth=(section.azimuth??35)*Math.PI/180,elevation=(section.elevation??30)*Math.PI/180;
 const horizontal=Math.cos(elevation);
 return [Math.sin(azimuth)*horizontal,-Math.sin(elevation),-Math.cos(azimuth)*horizontal];
}
export function sectionAngles(section:Pick<SectionState,'axis'|'azimuth'|'elevation'>){
 return section.axis==='axial'?{azimuth:0,elevation:90}:section.axis==='sagittal'?{azimuth:90,elevation:0}:section.axis==='coronal'?{azimuth:0,elevation:0}:{azimuth:section.azimuth??35,elevation:section.elevation??30};
}
export const sectionDot=(a:Point3,b:Point3)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
export function sectionRange(bounds:Bounds3,normal:Point3){
 const center:Point3=[0,1,2].map(i=>(bounds.min[i]+bounds.max[i])/2) as Point3;
 const radius=[0,1,2].reduce((sum,i)=>sum+Math.abs(normal[i])*(bounds.max[i]-bounds.min[i])/2,0);
 const middle=sectionDot(center,normal);
 return {min:middle-radius,max:middle+radius,center,radius};
}
export function sectionPoint(bounds:Bounds3,section:Pick<SectionState,'axis'|'azimuth'|'elevation'|'position'>):Point3{
 const normal=sectionNormal(section),range=sectionRange(bounds,normal),distance=range.min+clamp(section.position,0,1)*(range.max-range.min)-sectionDot(range.center,normal);
 return range.center.map((value,i)=>value+normal[i]*distance) as Point3;
}
export function sectionFraction(bounds:Bounds3,section:Pick<SectionState,'axis'|'azimuth'|'elevation'>,point:Point3){
 const normal=sectionNormal(section),range=sectionRange(bounds,normal);
 return range.radius<1e-8?.5:clamp((sectionDot(point,normal)-range.min)/(range.max-range.min),0,1);
}
export function sectionPlaneEquation(bounds:Bounds3,section:Pick<SectionState,'axis'|'azimuth'|'elevation'|'position'|'flip'>){
 const normal=sectionNormal(section),point=sectionPoint(bounds,section),sign=section.flip?1:-1;
 return {normal:normal.map(value=>value*sign) as Point3,constant:-sign*sectionDot(normal,point),point};
}
