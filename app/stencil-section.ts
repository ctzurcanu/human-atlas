import * as T from 'three';

/** A private, closed counting shell. The displayed geometry is never changed. */
export function stencilShell(source:T.BufferGeometry,inset=0):T.BufferGeometry{
 const position=source.getAttribute('position'),index=source.getIndex();
 if(!position||!index)throw new Error('Section counting needs indexed positions.');
 const weld=new Int32Array(position.count),representatives:number[]=[],vertices=new Map<string,number>();
 for(let i=0;i<position.count;i++){
  const key=`${Math.round(position.getX(i)*1e5)},${Math.round(position.getY(i)*1e5)},${Math.round(position.getZ(i)*1e5)}`;
  let id=vertices.get(key);if(id===undefined){id=representatives.length;vertices.set(key,id);representatives.push(i);}weld[i]=id;
 }
 const directed=new Map<string,number>();
 const directedKey=(a:number,b:number)=>`${a},${b}`;
 for(let i=0;i<index.count;i+=3){const ids=[weld[index.getX(i)],weld[index.getX(i+1)],weld[index.getX(i+2)]];
  if(ids[0]===ids[1]||ids[1]===ids[2]||ids[2]===ids[0])continue;
  for(let j=0;j<3;j++){const key=directedKey(ids[j],ids[(j+1)%3]);directed.set(key,(directed.get(key)??0)+1);}
 }
 const outgoing=new Map<number,number[]>();
 for(const [key,count] of directed){const [a,b]=key.split(',').map(Number),surplus=count-(directed.get(directedKey(b,a))??0);
  if(surplus<=0)continue;const ends=outgoing.get(a)??[];for(let i=0;i<surplus;i++)ends.push(b);outgoing.set(a,ends);
 }
 const extra:number[]=[];
 for(const start of outgoing.keys()){
  const starts=outgoing.get(start)!;
  while(starts.length){const path=[start];let next=starts.pop(),guard=0;
   while(next!==undefined&&next!==start&&guard++<directed.size){path.push(next);next=outgoing.get(next)?.pop();}
   if(next!==start||path.length<3)continue;
   const pivot=representatives[path[0]];
   for(let j=1;j<path.length-1;j++)extra.push(pivot,representatives[path[j+1]],representatives[path[j]]);
  }
 }
 const indices=new (position.count>65535?Uint32Array:Uint16Array)(index.count+extra.length);
 for(let i=0;i<index.count;i++)indices[i]=index.getX(i);
 indices.set(extra,index.count);
 const geometry=new T.BufferGeometry();geometry.setAttribute('position',position);geometry.setIndex(new T.BufferAttribute(indices,1));
 geometry.boundingBox=source.boundingBox?.clone()??null;geometry.boundingSphere=source.boundingSphere?.clone()??null;
 geometry.userData.originalTriangles=index.count/3;geometry.userData.repairTriangles=extra.length/3;
 if(inset>0)geometry.setAttribute('capInset',new T.Float32BufferAttribute(safeInset(position,indices,weld,representatives.length,index.count,inset),3));
 return geometry;
}

function safeInset(position:T.BufferAttribute|T.InterleavedBufferAttribute,indices:Uint16Array|Uint32Array,weld:Int32Array,weldCount:number,originalIndexCount:number,requested:number){
 const normal=new Float32Array(weldCount*3),a=new T.Vector3(),b=new T.Vector3(),c=new T.Vector3(),ab=new T.Vector3(),ac=new T.Vector3(),face=new T.Vector3();
 let volume=0;
 for(let i=0;i<indices.length;i+=3){const ia=indices[i],ib=indices[i+1],ic=indices[i+2];
  a.fromBufferAttribute(position,ia);b.fromBufferAttribute(position,ib);c.fromBufferAttribute(position,ic);
  face.crossVectors(ab.subVectors(b,a),ac.subVectors(c,a));volume+=a.x*(b.y*c.z-b.z*c.y)+a.y*(b.z*c.x-b.x*c.z)+a.z*(b.x*c.y-b.y*c.x);
  for(const j of [ia,ib,ic]){const k=weld[j]*3;normal[k]+=face.x;normal[k+1]+=face.y;normal[k+2]+=face.z;}
 }
 const sign=volume<0?-1:1;
 for(let i=0;i<weldCount;i++){const k=i*3,length=Math.hypot(normal[k],normal[k+1],normal[k+2])||1;normal[k]=normal[k]/length*sign;normal[k+1]=normal[k+1]/length*sign;normal[k+2]=normal[k+2]/length*sign;}
 const limits=new Float32Array(weldCount).fill(requested),p=[new T.Vector3(),new T.Vector3(),new T.Vector3()],q=[new T.Vector3(),new T.Vector3(),new T.Vector3()],before=new T.Vector3(),after=new T.Vector3(),u=new T.Vector3(),v=new T.Vector3();
 const du=new T.Vector3(),dv=new T.Vector3(),linear=new T.Vector3(),quadratic=new T.Vector3(),temp=new T.Vector3();
 // Bound each vertex by the first area reversal of any incident triangle.
 // A tube's local radius therefore limits its wall without thinning the
 // wider part of that same structure.
 for(let i=0;i<originalIndexCount;i+=3){
  const i0=indices[i],i1=indices[i+1],i2=indices[i+2],w0=weld[i0],w1=weld[i1],w2=weld[i2];
  p[0].fromBufferAttribute(position,i0);p[1].fromBufferAttribute(position,i1);p[2].fromBufferAttribute(position,i2);
  u.subVectors(p[1],p[0]);v.subVectors(p[2],p[0]);
  du.set(normal[w1*3]-normal[w0*3],normal[w1*3+1]-normal[w0*3+1],normal[w1*3+2]-normal[w0*3+2]);
  dv.set(normal[w2*3]-normal[w0*3],normal[w2*3+1]-normal[w0*3+1],normal[w2*3+2]-normal[w0*3+2]);
  before.crossVectors(u,v);const c0=before.lengthSq();if(c0<1e-24)continue;
  linear.crossVectors(du,v).add(temp.crossVectors(u,dv)).negate();quadratic.crossVectors(du,dv);
  const c1=before.dot(linear),c2=before.dot(quadratic);
  let root=Infinity;
  if(Math.abs(c2)<1e-20){if(c1<0)root=-c0/c1;}
  else{const discriminant=c1*c1-4*c2*c0;if(discriminant>=0){const sqrt=Math.sqrt(discriminant),r1=(-c1-sqrt)/(2*c2),r2=(-c1+sqrt)/(2*c2);if(r1>0)root=Math.min(root,r1);if(r2>0)root=Math.min(root,r2);}}
  if(root<Infinity){const safe=root*.8;limits[w0]=Math.min(limits[w0],safe);limits[w1]=Math.min(limits[w1],safe);limits[w2]=Math.min(limits[w2],safe);}
 }
 // Back off locally when an inset would flip an original triangle. This also
 // keeps small vessel or bone tips from acquiring an inverted, oversized cap.
 for(let pass=0;pass<7;pass++){let changed=false;
  for(let i=0;i<originalIndexCount;i+=3){const ids=[indices[i],indices[i+1],indices[i+2]],w=ids.map(id=>weld[id]);
   for(let j=0;j<3;j++){p[j].fromBufferAttribute(position,ids[j]);const k=w[j]*3,s=limits[w[j]];q[j].set(p[j].x-normal[k]*s,p[j].y-normal[k+1]*s,p[j].z-normal[k+2]*s);}
   before.crossVectors(u.subVectors(p[1],p[0]),v.subVectors(p[2],p[0]));
   if(before.lengthSq()<1e-20)continue;
   after.crossVectors(u.subVectors(q[1],q[0]),v.subVectors(q[2],q[0]));
   if(before.dot(after)>=before.lengthSq()*.25)continue;
   for(const id of w)limits[id]*=.5;changed=true;
  }
  if(!changed)break;
 }
 const displacement=new Float32Array(position.count*3);
 for(let i=0;i<position.count;i++){const k=weld[i]*3,target=i*3,limit=limits[weld[i]];displacement[target]=normal[k]*limit;displacement[target+1]=normal[k+1]*limit;displacement[target+2]=normal[k+2]*limit;}
 return displacement;
}
