import {BufferAttribute, BufferGeometry} from 'three';
import {createSharedTissueWall, type SharedWallInput} from './shared-tissue-wall';

export type TissueSubset = {id:string;layer:number};

/** Arbitrary disjoint tissue subsets per coat. Every contact is drawn only
 * when its actual neighbor is hidden. Counting shells remain closed. */
export function createPartitionedTissueWall(input:SharedWallInput, ownership:readonly Uint16Array[], subsets:readonly TissueSubset[]) {
  const L=input.layerIds.length,F=input.faces.length/3,n=input.interfaces[0].length/3;
  if(ownership.length!==L||new Set(subsets.map(s=>s.id)).size!==subsets.length
    ||subsets.some(s=>!Number.isInteger(s.layer)||s.layer<0||s.layer>=L)
    ||ownership.some((row,k)=>row.length!==F||row.some(i=>!subsets[i]||subsets[i].layer!==k)))throw Error('Invalid partitioned tissue ownership');
  const pool=createSharedTissueWall(input);
  const edgeRows=new Map<string,{face:number;a:number;b:number}[]>();
  for(let fi=0;fi<F;fi++)for(let j=0;j<3;j++) {
    const a=input.faces[fi*3+j],b=input.faces[fi*3+(j+1)%3],key=a<b?`${a}:${b}`:`${b}:${a}`;
    const rows=edgeRows.get(key)??[];rows.push({face:fi,a,b});edgeRows.set(key,rows);
  }
  for(const rows of edgeRows.values())if(rows.length>2||rows.length===2&&(rows[0].a!==rows[1].b||rows[0].b!==rows[1].a))throw Error('Nonmanifold partition carrier');
  type Contact={neighbor:number;material:number;start:number;count:number};
  const full:BufferGeometry[]=[],displayed:BufferGeometry[]=[],contacts:Contact[][]=[];
  const Index=pool.position.count<=65536?Uint16Array:Uint32Array;
  const groups=subsets.map(()=>new Map<string,{neighbor:number;material:number;indices:number[]}>());
  function add(owner:number,neighbor:number,material:number,indices:number[]) {
    const key=`${material}:${neighbor}`,group=groups[owner].get(key)??{neighbor,material,indices:[]};
    group.indices.push(...indices);groups[owner].set(key,group);
  }
  for(let k=0;k<L;k++) {
    for(let fi=0;fi<F;fi++) {
      const owner=ownership[k][fi],a=input.faces[fi*3],b=input.faces[fi*3+1],c=input.faces[fi*3+2];
      add(owner,k===0?-1:ownership[k-1][fi],0,[k*n+a,k*n+b,k*n+c]);
      add(owner,k===L-1?-1:ownership[k+1][fi],1,[(k+1)*n+a,(k+1)*n+c,(k+1)*n+b]);
    }
    for(const rows of edgeRows.values())for(const row of rows) {
      const owner=ownership[k][row.face],other=rows.find(r=>r!==row),neighbor=other?ownership[k][other.face]:-1;
      if(neighbor===owner)continue;
      const a=k*n+row.a,b=k*n+row.b,A=a+n,B=b+n;add(owner,neighbor,2,[b,a,A,b,A,B]);
    }
  }
  for(let i=0;i<subsets.length;i++) {
    const indices:number[]=[],cc:Contact[]=[];
    for(const group of groups[i].values()) {cc.push({neighbor:group.neighbor,material:group.material,start:indices.length,count:group.indices.length});for(const id of group.indices)indices.push(id);}
    if(!indices.length)throw Error('Empty tissue subset');
    const base=pool.full[subsets[i].layer];
    function make(){const g=new BufferGeometry();g.setAttribute('position',pool.position);g.setAttribute('normal',pool.normal);g.boundingBox=base.boundingBox!.clone();g.boundingSphere=base.boundingSphere!.clone();g.userData.pairedTissueWall=true;g.userData.tissueLayer=input.layerIds[subsets[i].layer];g.userData.tissueSubset=subsets[i].id;return g;}
    const g=make();g.setIndex(new BufferAttribute(new Index(indices),1));for(const c of cc)g.addGroup(c.start,c.count,c.material);full.push(g);
    const view=make();view.setIndex(g.index);displayed.push(view);contacts.push(cc);
  }
  pool.dispose();let signature='';
  return {full,displayed,position:pool.position,normal:pool.normal,contacts,
    bytes:{position:pool.position.array.byteLength,normal:pool.normal.array.byteLength,indices:full.reduce((sum,g)=>sum+g.index!.array.byteLength,0)},
    update(visible:readonly boolean[]) {
      if(visible.length!==subsets.length)throw Error('Invalid subset visibility');
      const key=visible.map(Number).join('');if(key===signature)return false;signature=key;
      displayed.forEach((g,i)=>{g.clearGroups();g.setDrawRange(0,visible[i]?Infinity:0);if(!visible[i])return;for(const c of contacts[i])if(c.neighbor<0||!visible[c.neighbor])g.addGroup(c.start,c.count,c.material);});
      return true;
    },
    dispose(){for(const g of [...full,...displayed])g.dispose();}
  };
}
