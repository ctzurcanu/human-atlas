// Join only existing, bit-exact nine-interface mouths. No source deformation.
export function joinBowelWallPorts(organs) {
 const mouths=[];
 const key=(o,i)=>Array.from({length:9},(_,k)=>Array.from(o.bank.subarray((k*o.n+i)*3,(k*o.n+i)*3+3)).join(',')).join('|');
 for(const o of organs){
  const adjacency=new Map();
  for(let j=0;j<o.rims.length;j+=2){const a=o.rims[j],b=o.rims[j+1];for(const [x,y] of [[a,b],[b,a]])adjacency.set(x,[...(adjacency.get(x)??[]),y]);}
  if([...adjacency.values()].some(a=>a.length!==2))throw Error('Bowel mouth is not a simple boundary');
  const unseen=new Set(adjacency.keys());
  while(unseen.size){const queue=[unseen.values().next().value],vertices=[];unseen.delete(queue[0]);while(queue.length){const i=queue.pop();vertices.push(i);for(const j of adjacency.get(i))if(unseen.delete(j))queue.push(j);}
   const keys=new Map(vertices.map(i=>[key(o,i),i]));if(keys.size!==vertices.length)throw Error('Duplicate paired mouth control');
   mouths.push({o,vertices,keys,adjacency});
  }
 }
 const joins=[],used=new Set();
 for(let a=0;a<mouths.length;a++)for(let b=a+1;b<mouths.length;b++){
  const A=mouths[a],B=mouths[b];if(A.o===B.o||A.keys.size!==B.keys.size||![...A.keys.keys()].every(k=>B.keys.has(k)))continue;
  if(used.has(a)||used.has(b))throw Error('Ambiguous exact bowel junction');
  const pairs=A.vertices.map(i=>[i,B.keys.get(key(A.o,i))]),map=new Map(pairs);
  // Matching controls must also have the same boundary-edge adjacency.
  for(const [i,j] of pairs)if(!A.adjacency.get(i).every(v=>B.adjacency.get(j).includes(map.get(v))))throw Error('Mouth topology differs');
  used.add(a);used.add(b);joins.push({A,B,pairs});
 }
 const normalSums=new Map();
 for(const o of new Set(joins.flatMap(j=>[j.A.o,j.B.o]))){
  const sums=new Float64Array(o.bank.length);
  for(let k=0;k<9;k++)for(let f=0;f<o.faces.length;f+=3){const ids=[o.faces[f],o.faces[f+1],o.faces[f+2]].map(i=>(k*o.n+i)*3),[a,b,c]=ids,p=o.bank,u=[0,1,2].map(d=>p[b+d]-p[a+d]),v=[0,1,2].map(d=>p[c+d]-p[a+d]),n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];for(const i of ids)for(let d=0;d<3;d++)sums[i+d]+=n[d];}
  normalSums.set(o,sums);
 }
 const evidence=[];
 for(const {A,B,pairs} of joins){
  for(let k=0;k<9;k++)for(const [a,b] of pairs){const ia=(k*A.o.n+a)*3,ib=(k*B.o.n+b)*3,n=[0,1,2].map(d=>normalSums.get(A.o)[ia+d]+normalSums.get(B.o)[ib+d]),len=Math.hypot(...n);if(!(len>0))throw Error('Undefined junction normal');for(let d=0;d<3;d++)A.o.normals[ia+d]=B.o.normals[ib+d]=Math.round(n[d]/len*32767);}
  const counts=[];
  for(const [mouth,other] of [[A,B],[B,A]]){
   const edgeKeys=new Set(mouth.vertices.flatMap(i=>mouth.adjacency.get(i).map(j=>[i,j].sort((a,b)=>a-b).join(':'))));let count=0;
   for(const cell of mouth.o.cells){
    const groups=[];
    for(const g of cell.part.tissueDrawGroups){if(!cell.portStarts.has(g.start)){groups.push(g);continue;}
     for(let at=g.start;at<g.start+g.count;at+=6){const six=Array.from(cell.indices.subarray(at,at+6)),base=[...new Set(six.map(i=>i%mouth.o.n))].sort((a,b)=>a-b),matched=base.length===2&&edgeKeys.has(base.join(':'));
      const neighbor=matched?other.o.cells.find(c=>c.k===cell.k&&!c.feature)?.part.id:g.neighbor;
      if(matched&&!neighbor)throw Error('Missing corresponding junction tissue');
      const last=groups.at(-1);if(last&&last.start+last.count===at&&last.neighbor===neighbor&&last.port)last.count+=6;else groups.push({start:at,count:6,...(neighbor?{neighbor}:{}),port:true});if(matched)count+=6;
     }
    }
    cell.portStarts=new Set(groups.filter(g=>g.port).map(g=>g.start));
    cell.part.tissueDrawGroups=groups.map(({port,...g})=>g);
   }
   counts.push(count);
  }
  if(counts[0]!==pairs.length*8*6||counts[1]!==counts[0])throw Error('Incomplete junction tissue rims '+JSON.stringify({organs:[A.o.id,B.o.id],counts,expected:pairs.length*8*6}));
  evidence.push({organs:[A.o.id,B.o.id],controls:pairs.length,interfaces:9,sharedNormals:true,annularIndicesPerSide:counts[0],pairs});
 }
 return evidence;
}
