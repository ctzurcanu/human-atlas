import {BufferAttribute, BufferGeometry, Box3, Sphere, Vector3} from 'three';
import {createSharedTissueWall, type SharedWallInput} from './shared-tissue-wall';

/** Source face ownership partitions tissue, never the lumen. All regions use
 * the same position and normal attributes, including at every shared port. */
export function createRegionalTissueWall(input: SharedWallInput, owners: Uint8Array, regionIds: readonly string[]) {
  if (owners.length !== input.faces.length/3 || owners.some(i=>i>=regionIds.length)
    || new Set(regionIds).size !== regionIds.length) throw Error('Invalid tissue regions');
  const pool=createSharedTissueWall(input),count=input.interfaces[0].length/3,L=input.layerIds.length;
  const edges=new Map<string,{a:number;b:number;region:number}[]>();
  const regionalFaces=regionIds.map(()=>[] as number[]);
  for(let j=0;j<owners.length;j++) {
    const region=owners[j],f=Array.from(input.faces.subarray(j*3,j*3+3));regionalFaces[region].push(...f);
    for(let k=0;k<3;k++){const a=f[k],b=f[(k+1)%3],key=a<b?`${a}:${b}`:`${b}:${a}`;
      const rows=edges.get(key)??[];rows.push({a,b,region});edges.set(key,rows);}
  }
  const ports=regionIds.map(()=>new Map<number,number[]>());
  for(const rows of edges.values()) {
    if(rows.length>2 || rows.length===2 && (rows[0].a!==rows[1].b || rows[0].b!==rows[1].a)) throw Error('Nonmanifold source interface');
    if(rows.length===2 && rows[0].region===rows[1].region)continue;
    for(const row of rows){const neighbor=rows.find(r=>r!==row)?.region??-1,arr=ports[row.region].get(neighbor)??[];arr.push(row.a,row.b);ports[row.region].set(neighbor,arr);}
  }
  const full:BufferGeometry[]=[],displayed:BufferGeometry[]=[],rimGroups:{neighbor:number;start:number;count:number}[][]=[];
  const Index=pool.position.count<=65536?Uint16Array:Uint32Array;
  for(let r=0;r<regionIds.length;r++)for(let k=0;k<L;k++) {
    const faces=regionalFaces[r],idx:number[]=[];
    for(const v of faces)idx.push(k*count+v);
    for(let j=0;j<faces.length;j+=3)for(const v of [faces[j],faces[j+2],faces[j+1]])idx.push((k+1)*count+v);
    const groups:{neighbor:number;start:number;count:number}[]=[];
    for(const [neighbor,port] of ports[r]){const start=idx.length;
      for(let j=0;j<port.length;j+=2){const a=k*count+port[j],b=k*count+port[j+1],A=a+count,B=b+count;idx.push(b,a,A,b,A,B);}
      groups.push({neighbor,start,count:idx.length-start});}
    const box=new Box3(),p=new Vector3();for(const id of idx)box.expandByPoint(p.fromBufferAttribute(pool.position,id));
    const make=()=>{const g=new BufferGeometry();g.setAttribute('position',pool.position);g.setAttribute('normal',pool.normal);g.boundingBox=box.clone();g.boundingSphere=box.getBoundingSphere(new Sphere());g.userData.pairedTissueWall=true;g.userData.tissueLayer=input.layerIds[k];g.userData.region=regionIds[r];return g;};
    const g=make();g.setIndex(new BufferAttribute(new Index(idx),1));g.addGroup(0,faces.length,0);g.addGroup(faces.length,faces.length,1);for(const group of groups)g.addGroup(group.start,group.count,2);full.push(g);
    const view=make();view.setIndex(g.index);displayed.push(view);rimGroups.push(groups);
  }
  // Whole-organ indices were used only to calculate one shared normal field.
  pool.dispose();let signature='';
  return {full,displayed,position:pool.position,normal:pool.normal,
    regionFaceCounts:regionalFaces.map(f=>f.length/3),
    ports:ports.map(m=>Array.from(m,([neighbor,edges])=>({neighbor,edges:edges.length/2}))),
    update(visible:readonly boolean[]){if(visible.length!==full.length)throw Error('Invalid regional visibility');const key=visible.map(Number).join('');if(key===signature)return false;signature=key;
      for(let r=0;r<regionIds.length;r++)for(let k=0;k<L;k++){const i=r*L+k,g=displayed[i],size=regionalFaces[r].length;g.clearGroups();g.setDrawRange(0,visible[i]?Infinity:0);if(!visible[i])continue;
        if(k===0||!visible[i-1])g.addGroup(0,size,0);
        if(k===L-1||!visible[i+1])g.addGroup(size,size,1);
        // Keep closed counting shells; hide only the displayed shared annuli.
        for(const group of rimGroups[i])if(group.neighbor<0||!visible[group.neighbor*L+k])g.addGroup(group.start,group.count,2);
      }return true;},
    dispose(){for(const g of [...full,...displayed])g.dispose();}
  };
}
