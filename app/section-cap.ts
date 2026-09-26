import * as T from 'three';

type Point=[number,number];
type Node={point:Point;edges:number[]};
type Edge={a:number;b:number;used:boolean};
type Loop={points:Point[];area:number;parent:number;depth:number};
type Contour={points:Point[];closed:boolean};
export type SectionCapOptions={thinShell?:boolean;hollowWall?:boolean;width?:number;outlineWidth?:number;closureWidth?:number;closureFraction?:number;outermostOnly?:boolean};
const area=(points:Point[])=>points.reduce((sum,p,i)=>{const q=points[(i+1)%points.length];return sum+p[0]*q[1]-q[0]*p[1];},0)/2;
function contains(points:Point[],p:Point){let inside=false;for(let i=0,j=points.length-1;i<points.length;j=i++){const a=points[i],b=points[j];if((a[1]>p[1])!==(b[1]>p[1])&&p[0]<(b[0]-a[0])*(p[1]-a[1])/(b[1]-a[1])+a[0])inside=!inside;}return inside;}

// Intersect a source mesh with a plane and triangulate its closed contours.
// Keeping each source mesh separate prevents neighboring tissues from sharing a cap.
export function sectionCapGeometry(source:T.BufferGeometry,worldPlane:T.Plane,matrix:T.Matrix4,options:SectionCapOptions={}):T.BufferGeometry|null{
 const positions=source.getAttribute('position');const indices=source.getIndex();if(!positions||!indices)return null;
 const localPlane=worldPlane.clone().applyMatrix4(matrix.clone().invert()),normal=worldPlane.normal;
 const origin=normal.clone().multiplyScalar(-worldPlane.constant);
 const u=new T.Vector3().crossVectors(Math.abs(normal.y)>.9?new T.Vector3(0,0,1):new T.Vector3(0,1,0),normal).normalize();
 const v=new T.Vector3().crossVectors(normal,u).normalize();
 const nodes:Node[]=[],edges:Edge[]=[],nodeIds=new Map<string,number>();
 const quant=1e-4;
 const node=(p:Point)=>{const key=Math.round(p[0]/quant)+','+Math.round(p[1]/quant);let id=nodeIds.get(key);if(id===undefined){id=nodes.length;nodeIds.set(key,id);nodes.push({point:p,edges:[]});}return id;};
 const add=(p:Point,q:Point)=>{const a=node(p),b=node(q);if(a===b)return;const id=edges.length;edges.push({a,b,used:false});nodes[a].edges.push(id);nodes[b].edges.push(id);};
 const project=(p:T.Vector3):Point=>{p.applyMatrix4(matrix).sub(origin);return [p.dot(u),p.dot(v)];};
 const crossing=(a:number,b:number,da:number,db:number,out:Point[])=>{
  if(da*db>0||Math.abs(da)<1e-9&&Math.abs(db)<1e-9)return;
  const t=Math.max(0,Math.min(1,da/(da-db))),p=new T.Vector3(
   positions.getX(a)+(positions.getX(b)-positions.getX(a))*t,
   positions.getY(a)+(positions.getY(b)-positions.getY(a))*t,
   positions.getZ(a)+(positions.getZ(b)-positions.getZ(a))*t);
  const projected=project(p);if(!out.some(q=>Math.hypot(q[0]-projected[0],q[1]-projected[1])<quant/2))out.push(projected);
 };
 for(let i=0;i<indices.count;i+=3){
  const a=indices.getX(i),b=indices.getX(i+1),c=indices.getX(i+2);
  const da=localPlane.distanceToPoint(new T.Vector3(positions.getX(a),positions.getY(a),positions.getZ(a)));
  const db=localPlane.distanceToPoint(new T.Vector3(positions.getX(b),positions.getY(b),positions.getZ(b)));
  const dc=localPlane.distanceToPoint(new T.Vector3(positions.getX(c),positions.getY(c),positions.getZ(c)));
  if(Math.min(da,db,dc)>0||Math.max(da,db,dc)<0)continue;
  const hits:Point[]=[];crossing(a,b,da,db,hits);crossing(b,c,db,dc,hits);crossing(c,a,dc,da,hits);
  if(hits.length===2)add(hits[0],hits[1]);
 }
 if(!edges.length)return null;
 const loops:Loop[]=[],ribbons:Contour[]=[],contours:Contour[]=[];
 // An open contour must be walked from an endpoint. Starting in its middle
 // splits it into two unusable paths, which leaves otherwise cap-shaped cuts open.
 const starts=[...nodes.keys()].filter(id=>nodes[id].edges.length===1).concat(edges.map(edge=>edge.a));
 for(const start of starts){
  if(!nodes[start].edges.some(id=>!edges[id].used))continue;
  const path:Point[]=[];let current=start,previous=-1,guard=0;
  while(guard++<edges.length+1){
   path.push(nodes[current].point);
   const choices=nodes[current].edges.filter(id=>!edges[id].used);
   if(!choices.length)break;
   let nextEdge=choices[0];
   if(choices.length>1&&previous>=0){
    const before=nodes[previous].point,here=nodes[current].point;
    let best=-Infinity;for(const id of choices){const edge=edges[id],other=nodes[edge.a===current?edge.b:edge.a].point;
     const dx=here[0]-before[0],dy=here[1]-before[1],ex=other[0]-here[0],ey=other[1]-here[1],score=(dx*ex+dy*ey)/(Math.hypot(dx,dy)*Math.hypot(ex,ey)+1e-12);
     if(score>best){best=score;nextEdge=id;}}
   }
   edges[nextEdge].used=true;const edge=edges[nextEdge],next=edge.a===current?edge.b:edge.a;previous=current;current=next;
   if(current===start)break;
  }
  if(path.length<2)continue;
  contours.push({points:path,closed:current===start});
 }
 // A publisher mesh can split one real boundary into several open chains at
 // vertices shared by more than two cut edges. Reconnect only endpoints that
 // already meet at the cut. This recovers the two sides of bones and muscles
 // without drawing a new edge across an open sheet.
 const stitchDistance=options.thinShell ? .002 : options.closureWidth??.002;
 for(let changed=true;changed;){
  changed=false;let best=stitchDistance,bestPair:[number,number,boolean,boolean]|null=null;
  for(let i=0;i<contours.length;i++){
   if(contours[i].closed)continue;
   for(let j=i+1;j<contours.length;j++){
    if(contours[j].closed)continue;
    for(const reverseA of [false,true])for(const reverseB of [false,true]){
     const a=contours[i].points[reverseA?0:contours[i].points.length-1];
     const b=contours[j].points[reverseB?contours[j].points.length-1:0];
     const distance=Math.hypot(a[0]-b[0],a[1]-b[1]);
     if(distance<best){best=distance;bestPair=[i,j,reverseA,reverseB];}
    }
   }
  }
  if(!bestPair)break;
  const [i,j,reverseA,reverseB]=bestPair;
  const first=reverseA?[...contours[i].points].reverse():contours[i].points;
  const second=reverseB?[...contours[j].points].reverse():contours[j].points;
  contours[i]={points:[...first,...second.slice(best<quant?1:0)],closed:false};
  contours.splice(j,1);changed=true;
 }
 for(const {points:path,closed:tracedClosed} of contours){
  // A long closing edge invents tissue across empty space. Only weld a seam
  // that is small at the scale of its existing contour.
  const lo=[Infinity,Infinity],hi=[-Infinity,-Infinity];
  for(const p of path)for(let axis=0;axis<2;axis++){lo[axis]=Math.min(lo[axis],p[axis]);hi[axis]=Math.max(hi[axis],p[axis]);}
  const span=Math.hypot(hi[0]-lo[0],hi[1]-lo[1]);
  const first=path[0],last=path[path.length-1];
  const gap=Math.hypot(last[0]-first[0],last[1]-first[1]);
  const closed=tracedClosed||gap<=Math.min(stitchDistance,span*(options.closureFraction??.04));
  if(options.thinShell||!closed||path.length<3){ribbons.push({points:path,closed});continue;}
  const signed=area(path);
  const rectangle=(hi[0]-lo[0])*(hi[1]-lo[1]);
  if(Math.abs(signed)<1e-8||!tracedClosed&&Math.abs(signed)<rectangle*.025){ribbons.push({points:path,closed});continue;}
  loops.push({points:path,area:Math.abs(signed),parent:-1,depth:0});
 }
 if(options.outermostOnly){
  // The publisher's derived body surface includes nested cut contours. Only
  // its exterior boundary represents skin; inner contours are not another
  // layer of skin and must not form a second ring around the body cavity.
  const exterior=ribbons.filter((contour,i)=>!ribbons.some((other,j)=>{
   if(j===i||!other.closed||Math.abs(area(other.points))<=Math.abs(area(contour.points))+1e-7)return false;
   const samples=contour.closed?[contour.points[0]]:[contour.points[Math.floor(contour.points.length*.25)],contour.points[Math.floor(contour.points.length*.5)],contour.points[Math.floor(contour.points.length*.75)]];
   return samples.every(point=>contains(other.points,point));
  }));
  ribbons.splice(0,ribbons.length,...exterior);
 }
 if(!loops.length&&!ribbons.length)return null;
 for(let i=0;i<loops.length;i++){
  let closest=-1,size=Infinity;
  for(let j=0;j<loops.length;j++)if(i!==j&&loops[j].area>loops[i].area&&loops[j].area<size&&contains(loops[j].points,loops[i].points[0])){closest=j;size=loops[j].area;}
  loops[i].parent=closest;
 }
 for(const loop of loops){let parent=loop.parent;while(parent>=0){loop.depth++;parent=loops[parent].parent;}}
 const vertices:number[]=[],normals:number[]=[],uvs:number[]=[],triangleIndices:number[]=[];
 const append=(p:Point)=>{const world=origin.clone().addScaledVector(u,p[0]).addScaledVector(v,p[1]).addScaledVector(normal,-.00015);vertices.push(world.x,world.y,world.z);normals.push(normal.x,normal.y,normal.z);uvs.push(p[0],p[1]);};
 for(let i=0;i<loops.length;i++){
  const outer=loops[i];if(outer.depth%2)continue;
  const holes=loops.filter(loop=>loop.parent===i&&loop.depth%2).map(loop=>loop.points);
  if(options.hollowWall&&!holes.length){ribbons.push({points:outer.points,closed:true});continue;}
  const points=[...outer.points,...holes.flat()],faces=T.ShapeUtils.triangulateShape(outer.points.map(p=>new T.Vector2(...p)),holes.map(hole=>hole.map(p=>new T.Vector2(...p))));
  if(!faces.length){ribbons.push({points:outer.points,closed:true},...holes.map(points=>({points,closed:true})));continue;}
  const expected=outer.area-holes.reduce((sum,hole)=>sum+Math.abs(area(hole)),0);
  const triangulated=faces.reduce((sum,[a,b,c])=>sum+Math.abs((points[b][0]-points[a][0])*(points[c][1]-points[a][1])-(points[b][1]-points[a][1])*(points[c][0]-points[a][0]))*.5,0);
  if(expected<=0||Math.abs(triangulated-expected)>Math.max(.000002,expected*.08)){
   ribbons.push({points:outer.points,closed:true},...holes.map(points=>({points,closed:true})));
   continue;
  }
  const offset=vertices.length/3;
  for(const p of points)append(p);
  for(const face of faces)triangleIndices.push(offset+face[0],offset+face[1],offset+face[2]);
  if(options.outlineWidth){
   ribbons.push({points:outer.points,closed:true},...holes.map(points=>({points,closed:true})));
  }
 }
 const fillIndexCount=triangleIndices.length;
 for(const {points,closed} of ribbons){
  if(points.length<2)continue;
  const lo=[Infinity,Infinity],hi=[-Infinity,-Infinity];for(const p of points)for(let axis=0;axis<2;axis++){lo[axis]=Math.min(lo[axis],p[axis]);hi[axis]=Math.max(hi[axis],p[axis]);}
  const width=Math.min(options.outlineWidth??options.width??.0015,Math.max(.0002,Math.hypot(hi[0]-lo[0],hi[1]-lo[1])*.1));
  const offset=vertices.length/3;
  for(let i=0;i<points.length;i++){
   const before=points[i===0?(closed?points.length-1:0):i-1],after=points[i===points.length-1?(closed?0:points.length-1):i+1];
   const dx=after[0]-before[0],dy=after[1]-before[1],length=Math.hypot(dx,dy)||1;
   const nx=-dy/length*width*.5,ny=dx/length*width*.5;
   append([points[i][0]+nx,points[i][1]+ny]);append([points[i][0]-nx,points[i][1]-ny]);
  }
  const segments=closed?points.length:points.length-1;
  for(let i=0;i<segments;i++){const a=offset+i*2,b=offset+((i+1)%points.length)*2;triangleIndices.push(a,a+1,b,b,a+1,b+1);}
 }
 if(!triangleIndices.length)return null;
 const cap=new T.BufferGeometry();cap.setAttribute('position',new T.Float32BufferAttribute(vertices,3));cap.setAttribute('normal',new T.Float32BufferAttribute(normals,3));cap.setAttribute('cutUv',new T.Float32BufferAttribute(uvs,2));cap.setAttribute('uv',new T.Float32BufferAttribute(uvs.map(value=>value*50),2));cap.setIndex(triangleIndices);
 cap.userData={fillIndexCount,outlineIndexCount:triangleIndices.length-fillIndexCount,closedContours:loops.length,openContours:ribbons.filter(contour=>!contour.closed).length};
 if(options.outlineWidth){if(fillIndexCount)cap.addGroup(0,fillIndexCount,0);if(triangleIndices.length>fillIndexCount)cap.addGroup(fillIndexCount,triangleIndices.length-fillIndexCount,1);}
 cap.computeBoundingSphere();return cap;
}
