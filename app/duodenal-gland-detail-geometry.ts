/** Original branched submucosal gland illustration, not measured organ anatomy.
 * Values used by the generator are millimetres; output coordinates are metres.
 * The conforming tetrahedral extraction bounds tissue, never fills its lumen.
 */
type Point=[number,number,number];
const branches:{a:Point;b:Point;r:number}[]=[{a:[0,.04,0],b:[0,-.12,0],r:.034}];
for(const side of [-1,1]){
 const junction:Point=[side*.07,-.20,side*.016];branches.push({a:[0,-.12,0],b:junction,r:.035});
 for(const [x,y,z] of [[.14,-.28,.065],[.14,-.28,-.065]])branches.push({a:junction,b:[side*x,y,side*z],r:.043});
}
const mouth=.00317,thickness=.018;
function capsule(p:readonly number[],a:Point,b:Point,r:number){const d=b.map((v,i)=>v-a[i]),t=Math.max(0,Math.min(1,d.reduce((s,v,i)=>s+(p[i]-a[i])*v,0)/d.reduce((s,v)=>s+v*v,0)));return Math.hypot(...d.map((v,i)=>p[i]-a[i]-v*t))-r;}
export function duodenalGlandField(point:readonly number[]){const p=point.map(v=>v*1000),outer=Math.min(...branches.map(b=>capsule(p,b.a,b.b,b.r))),lumen=Math.min(...branches.map(b=>capsule(p,b.a,b.b,b.r-thickness)));return Math.max(outer,-lumen,p[1]-mouth);}
export function duodenalGlandDetailGeometry(){
 const step=.008,min:Point=[-.21389,-.34273,-.11011],dims=[55,49,29],grid=new Float64Array(dims[0]*dims[1]*dims[2]);
 const id=(x:number,y:number,z:number)=>(z*dims[1]+y)*dims[0]+x;
 const coord=(index:number):Point=>{const x=index%dims[0],y=Math.floor(index/dims[0])%dims[1],z=Math.floor(index/(dims[0]*dims[1]));return [min[0]+x*step,min[1]+y*step,min[2]+z*step];};
 for(let i=0;i<grid.length;i++)grid[i]=duodenalGlandField(coord(i).map(v=>v*.001));
 const positions:number[]=[],faces:number[]=[],cache=new Map<string,number>();
 const vertex=(a:number,b:number)=>{const key=Math.min(a,b)+':'+Math.max(a,b),existing=cache.get(key);if(existing!==undefined)return existing;// Keep isovertex distinct from grid endpoints (maximum 13.86 nm shift).
  const A=coord(a),B=coord(b),t=Math.max(.001,Math.min(.999,grid[a]/(grid[a]-grid[b]))),index=positions.length/3;positions.push(...A.map((v,i)=>(v+(B[i]-v)*t)*.001));cache.set(key,index);return index;};
 const tetrahedra=[[0,1,3,7],[0,3,2,7],[0,2,6,7],[0,6,4,7],[0,4,5,7],[0,5,1,7]],edges=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];
 const sub=(a:readonly number[],b:readonly number[])=>a.map((v,i)=>v-b[i]),dot=(a:number[],b:number[])=>a.reduce((s,v,i)=>s+v*b[i],0),cross=(a:number[],b:number[])=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],p=(i:number)=>positions.slice(i*3,i*3+3);
 for(let z=0;z<dims[2]-1;z++)for(let y=0;y<dims[1]-1;y++)for(let x=0;x<dims[0]-1;x++){
  const cube=[id(x,y,z),id(x+1,y,z),id(x,y+1,z),id(x+1,y+1,z),id(x,y,z+1),id(x+1,y,z+1),id(x,y+1,z+1),id(x+1,y+1,z+1)];
  for(const tet of tetrahedra){const ids=tet.map(k=>cube[k]),inside=ids.filter(k=>grid[k]<0),outside=ids.filter(k=>grid[k]>=0);if(!inside.length||!outside.length)continue;
   const polygon=edges.filter(([a,b])=>(grid[ids[a]]<0)!==(grid[ids[b]]<0)).map(([a,b])=>vertex(ids[a],ids[b]));
   const average=(indices:number[])=>[0,1,2].map(axis=>indices.reduce((s,i)=>s+coord(i)[axis],0)/indices.length),n=sub(average(outside),average(inside));
   const center=[0,1,2].map(axis=>polygon.reduce((s,i)=>s+p(i)[axis],0)/polygon.length),u=sub(p(polygon[0]),center),v=cross(n,u);
   polygon.sort((a,b)=>Math.atan2(dot(sub(p(a),center),v),dot(sub(p(a),center),u))-Math.atan2(dot(sub(p(b),center),v),dot(sub(p(b),center),u)));
   for(let j=1;j<polygon.length-1;j++){const tri=[polygon[0],polygon[j],polygon[j+1]],normal=cross(sub(p(tri[1]),p(tri[0])),sub(p(tri[2]),p(tri[0])));if(dot(normal,n)<0)[tri[1],tri[2]]=[tri[2],tri[1]];faces.push(...tri);}
  }
 }
 return {positions:new Float32Array(positions),faces:new Uint16Array(faces),stepMeters:step*.001,illustrativeThicknessMeters:thickness*.001,mouthMeters:mouth*.001,branches:branches.length};
}
