import type * as T from 'three';

export type MeshTopology={triangles:number;boundaryEdges:number;nonManifoldEdges:number;windingConflicts:number;duplicateTriangles:number;degenerateTriangles:number;reversedNormals:number;hasUv:boolean;closed:boolean;boundaryRatio:number};

// Weld only for diagnosis. The renderer retains the publisher's original
// vertices, normals, UV seams, and triangles.
export function analyzeMeshTopology(geometry:T.BufferGeometry):MeshTopology{
 const position=geometry.getAttribute('position'),normal=geometry.getAttribute('normal'),index=geometry.getIndex();
 const empty={triangles:0,boundaryEdges:0,nonManifoldEdges:0,windingConflicts:0,duplicateTriangles:0,degenerateTriangles:0,reversedNormals:0,hasUv:!!geometry.getAttribute('uv'),closed:false,boundaryRatio:1};
 if(!position||!index)return empty;
 const box=geometry.boundingBox??(geometry.computeBoundingBox(),geometry.boundingBox!);
 const extent=Math.hypot(box.max.x-box.min.x,box.max.y-box.min.y,box.max.z-box.min.z);
 const unit=Math.max(1e-7,extent*1e-5),ids=new Int32Array(position.count),welded=new Map<string,number>();
 for(let i=0;i<position.count;i++){
  const key=`${Math.round(position.getX(i)/unit)},${Math.round(position.getY(i)/unit)},${Math.round(position.getZ(i)/unit)}`;
  let id=welded.get(key);if(id===undefined){id=welded.size;welded.set(key,id);}ids[i]=id;
 }
 const edges=new Map<string,{count:number;balance:number}>(),faces=new Set<string>();
 let duplicateTriangles=0,degenerateTriangles=0,reversedNormals=0;
 for(let i=0;i<index.count;i+=3){
  const raw=[index.getX(i),index.getX(i+1),index.getX(i+2)],v=raw.map(j=>ids[j]);
  if(new Set(v).size<3){degenerateTriangles++;continue;}
  const ax=position.getX(raw[0]),ay=position.getY(raw[0]),az=position.getZ(raw[0]);
  const bx=position.getX(raw[1])-ax,by=position.getY(raw[1])-ay,bz=position.getZ(raw[1])-az;
  const cx=position.getX(raw[2])-ax,cy=position.getY(raw[2])-ay,cz=position.getZ(raw[2])-az;
  const nx=by*cz-bz*cy,ny=bz*cx-bx*cz,nz=bx*cy-by*cx;
  if(nx*nx+ny*ny+nz*nz<unit*unit*unit*unit){degenerateTriangles++;continue;}
  if(normal){const dot=nx*(normal.getX(raw[0])+normal.getX(raw[1])+normal.getX(raw[2]))+ny*(normal.getY(raw[0])+normal.getY(raw[1])+normal.getY(raw[2]))+nz*(normal.getZ(raw[0])+normal.getZ(raw[1])+normal.getZ(raw[2]));if(dot<0)reversedNormals++;}
  const face=[...v].sort((a,b)=>a-b).join(',');if(faces.has(face))duplicateTriangles++;else faces.add(face);
  for(let j=0;j<3;j++){const a=v[j],b=v[(j+1)%3],key=a<b?`${a},${b}`:`${b},${a}`,edge=edges.get(key)??{count:0,balance:0};edge.count++;edge.balance+=a<b?1:-1;edges.set(key,edge);}
 }
 let boundaryEdges=0,nonManifoldEdges=0,windingConflicts=0;
 for(const edge of edges.values()){if(edge.count===1)boundaryEdges++;else if(edge.count>2)nonManifoldEdges++;else if(edge.balance!==0)windingConflicts++;}
 return {triangles:index.count/3,boundaryEdges,nonManifoldEdges,windingConflicts,duplicateTriangles,degenerateTriangles,reversedNormals,hasUv:!!geometry.getAttribute('uv'),closed:boundaryEdges===0&&nonManifoldEdges===0&&windingConflicts===0,boundaryRatio:edges.size?boundaryEdges/edges.size:1};
}
