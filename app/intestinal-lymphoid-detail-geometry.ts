/** Original schematic secondary lymphoid follicle in metres. The shared
 * convex boundaries delimit solid germinal-center and mantle tissue, no lumen.
 * Shape, dimensions, six-member layout and central symmetry are illustrative.
 */
export function intestinalLymphoidDetailGeometry(){
 const sides=32,rings=13,count=sides*rings+2,top=count-2,bottom=count-1;
 const radii=[.00048,.00038,.00034],outer=new Float32Array(count*3),core=new Float32Array(count*3),faces:number[]=[];
 for(let k=0;k<rings;k++)for(let j=0;j<sides;j++){
  const latitude=(k+1)/(rings+1)*Math.PI,angle=j/sides*2*Math.PI,id=k*sides+j;
  const direction=[Math.sin(latitude)*Math.cos(angle),Math.cos(latitude),Math.sin(latitude)*Math.sin(angle)];
  outer.set(direction.map((v,d)=>v*radii[d]),id*3);
  // Power-of-two homothety preserves the exact float32 shared topology.
  core.set(Array.from(outer.slice(id*3,id*3+3),v=>v*.5),id*3);
 }
 outer.set([0,radii[1],0],top*3);outer.set([0,-radii[1],0],bottom*3);
 core.set(Array.from(outer.slice(top*3,top*3+3),v=>v*.5),top*3);core.set(Array.from(outer.slice(bottom*3,bottom*3+3),v=>v*.5),bottom*3);
 for(let j=0;j<sides;j++){
  const a=j,b=(j+1)%sides;faces.push(top,b,a);
  for(let k=0;k<rings-1;k++){
   const a=k*sides+j,b=k*sides+(j+1)%sides,A=a+sides,B=b+sides;
   // The ideal latitude quad is planar. Packing can make its four corners
   // slightly nonplanar; choose the convex diagonal from actual stored data.
   const u=[0,1,2].map(d=>outer[b*3+d]-outer[a*3+d]),v=[0,1,2].map(d=>outer[A*3+d]-outer[a*3+d]);
   const n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]],side=n.reduce((s,x,d)=>s+x*(outer[B*3+d]-outer[a*3+d]),0);
   if(side>0)faces.push(a,b,B,a,B,A);else faces.push(a,b,A,b,B,A);
  }
  const A=(rings-1)*sides+j,B=(rings-1)*sides+(j+1)%sides;faces.push(bottom,A,B);
 }
 const centers:readonly (readonly number[])[]=[[-.00115,-.00052,0],[0,-.00052,0],[.00115,-.00052,0],[-.00115,.00052,0],[0,.00052,0],[.00115,.00052,0]];
 return {outer,core,faces:new Uint32Array(faces),radii,centers,count,sides};
}

/** Exact convex polyhedral containment for cut-face picking, rather than an
 * ellipsoid approximation that would extend a cap beyond the rendered wall.
 */
export function insideLymphoidEnvelope(point:readonly number[],positions:Float32Array,faces:Uint32Array){
 for(let i=0;i<faces.length;i+=3){
  const a=faces[i]*3,b=faces[i+1]*3,c=faces[i+2]*3;
  const u=[positions[b]-positions[a],positions[b+1]-positions[a+1],positions[b+2]-positions[a+2]],v=[positions[c]-positions[a],positions[c+1]-positions[a+1],positions[c+2]-positions[a+2]];
  const n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];
  const distance=n.reduce((s,x,d)=>s+x*(point[d]-positions[a+d]),0);
  if(distance>Math.hypot(...n)*1e-12)return false;
 }
 return true;
}
