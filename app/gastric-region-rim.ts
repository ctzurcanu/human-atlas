import {BufferAttribute,BufferGeometry} from 'three';

/** Bridge paired gastric wall boundaries; never fill a lumen opening. */
export function gastricRegionRim(source:BufferGeometry):BufferGeometry|null{
 const position=source.getAttribute('position'),indices=source.index;if(!indices)return null;
 const points:number[][]=[],canonical=new Map<string,number>(),vertices:number[]=[];
 for(let i=0;i<position.count;i++){
  const p=[position.getX(i),position.getY(i),position.getZ(i)],key=p.map(v=>Math.round(v*1e7)).join(',');
  let id=canonical.get(key);if(id===undefined){id=points.length;canonical.set(key,id);points.push(p);}vertices.push(id);
 }
 const edges=new Map<string,{a:number;b:number;count:number}>();
 for(let i=0;i<indices.count;i+=3){const face=[0,1,2].map(j=>vertices[indices.getX(i+j)]);
  for(let j=0;j<3;j++){const a=face[j],b=face[(j+1)%3],key=a<b?`${a}:${b}`:`${b}:${a}`,edge=edges.get(key);if(edge)edge.count++;else edges.set(key,{a,b,count:1});}
 }
 const boundary=[...edges.values()].filter(edge=>edge.count===1),next=new Map(boundary.map(edge=>[edge.a,edge.b]));
 if(next.size!==boundary.length)return null;
 const remaining=new Set(next.keys()),loops:number[][]=[];
 while(remaining.size){const start=remaining.values().next().value!;let current=start;const loop:number[]=[];
  do{if(!remaining.delete(current))return null;loop.push(current);const following=next.get(current);if(following===undefined)return null;current=following;}while(current!==start);
  loops.push(loop);
 }
 if(!loops.length||loops.length%2)return null;
 const center=(loop:number[])=>[0,1,2].map(k=>loop.reduce((sum,id)=>sum+points[id][k],0)/loop.length);
 const distance=(a:number[],b:number[])=>a.reduce((sum,v,k)=>sum+(v-b[k])**2,0);
 const pending=[...loops],pairs:[number[],number[]][]=[],out:number[]=[];
 while(pending.length){const a=pending.shift()!,ca=center(a);let nearest=0,best=Infinity;
  pending.forEach((b,i)=>{const cost=distance(ca,center(b));if(cost<best){nearest=i;best=cost;}});
  pairs.push([a,pending.splice(nearest,1)[0]]);
 }
 for(const [a,b] of pairs){
  const n=a.length,m=b.length;
  if(n===m){
   // Preserve the existing equal-edge pyloric rim and its diagonal choices.
   let offset=0,best=Infinity;
   for(let shift=0;shift<n;shift++){let cost=0;for(let i=0;i<n;i++)cost+=distance(points[a[i]],points[b[(shift-i+n)%n]]);if(cost<best){best=cost;offset=shift;}}
   if(Math.sqrt(best/n)>.006)return null;
   const paired=a.map((_,i)=>b[(offset-i+n)%n]);
   for(let i=0;i<n;i++){const j=(i+1)%n;for(const id of [a[j],a[i],paired[i],a[j],paired[i],paired[j]])out.push(...points[id]);}
  }else{
   // A cyclic zipper accommodates different source sampling on the two layers.
   // It uses only existing positions and joins opposite boundary windings.
   let best=Infinity,chosen:number[]=[],steps:Uint8Array|undefined;
   for(let offset=0;offset<m;offset++){
    const bb=b.map((_,j)=>b[(offset-j+m)%m]),width=m+1,dp=new Float64Array((n+1)*width).fill(Infinity),prev=new Uint8Array(dp.length);dp[0]=0;
    for(let i=0;i<=n;i++)for(let j=0;j<=m;j++){
     const k=i*width+j;
     if(i<n){const next=(i+1)*width+j,cost=dp[k]+distance(points[a[(i+1)%n]],points[bb[j%m]]);if(cost<dp[next]){dp[next]=cost;prev[next]=1;}}
     if(j<m){const next=k+1,cost=dp[k]+distance(points[a[i%n]],points[bb[(j+1)%m]]);if(cost<dp[next]){dp[next]=cost;prev[next]=2;}}
    }
    if(dp[n*width+m]<best){best=dp[n*width+m];chosen=bb;steps=prev;}
   }
   // Distant rings are different physiological ports, not two wall layers.
   if(!steps||Math.sqrt(best/(n+m))>.006)return null;
   let i=n,j=m;
   while(i||j){if(steps[i*(m+1)+j]===1){for(const id of [a[i%n],a[(i-1)%n],chosen[j%m]])out.push(...points[id]);i--;}
    else{for(const id of [a[i%n],chosen[(j-1)%m],chosen[j%m]])out.push(...points[id]);j--;}}
  }
 }
 const rim=new BufferGeometry();rim.setAttribute('position',new BufferAttribute(new Float32Array(out),3));rim.computeVertexNormals();rim.computeBoundingSphere();return rim;
}
