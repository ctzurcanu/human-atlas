import {SYSTEMS,type Atlas,type Part} from './anatomy';
import {REGION_ORDER,buildAnatomyNodes,systemName,type AnatomyEntry,type AnatomyNode} from './anatomy-hierarchy';
import {DEPTH_LAYERS,depthLayerFor} from './depth-layers';
import type {ResolvedGuestNode} from './guest-hierarchy';

export type ExplodeHierarchy='systems'|'regions'|'depth'|'guest';
type Bounds={cx:number;cy:number;cz:number;width:number;height:number};
type Node={parts:Part[];children:Node[];bounds:Bounds;name?:string};
type Measured={node:Node;width:number;height:number;children:{measured:Measured;x:number;y:number}[]|null};
export type ExplosionGroup={id:string;name:string;depth:number;x:number;y:number;z:number;height:number;count:number};
export type ExplosionStage={positions:Float32Array;width:number;height:number;groups:ExplosionGroup[];clusterIds:Int32Array;clusters:number[][]};
export type HierarchicalExplosionLayout={stages:ExplosionStage[];steps:number};

const partCenter=(part:Part,axis:number)=>(part.bounds[0][axis]+part.bounds[1][axis])/2;
function nodeBounds(parts:Part[]):Bounds{
 const low=[Infinity,Infinity,Infinity],high=[-Infinity,-Infinity,-Infinity];
 for(const part of parts)for(let axis=0;axis<3;axis++){low[axis]=Math.min(low[axis],part.bounds[0][axis]);high[axis]=Math.max(high[axis],part.bounds[1][axis]);}
 return {cx:(low[0]+high[0])/2,cy:(low[1]+high[1])/2,cz:(low[2]+high[2])/2,width:high[0]-low[0],height:high[1]-low[1]};
}
const node=(parts:Part[],children:Node[]=[],name?:string):Node=>({parts,children,name,bounds:nodeBounds(parts)});
function entryNode(entry:AnatomyEntry):Node{
 return node(entry.parts,entry.parts.map(part=>node([part])));
}
function anatomyNode(item:AnatomyNode):Node{
 return item.kind==='group'?node(item.parts,item.nodes.map(anatomyNode),item.name)
  :item.kind==='bilateral'?node(item.parts,item.entries.map(entryNode),item.name):entryNode(item.entry);
}
function hierarchyTree(atlas:Atlas,entries:AnatomyEntry[],visibleIds:Set<string>,mode:ExplodeHierarchy,guestNodes?:ResolvedGuestNode[]):Node|null{
 const visibleEntries=entries.flatMap(entry=>{const parts=entry.parts.filter(part=>visibleIds.has(part.id));return parts.length?[{...entry,parts}]:[]});
 if(!visibleEntries.length)return null;
 if(mode==='guest'&&guestNodes){
  const remaining=new Set(visibleIds);
  const convert=(item:ResolvedGuestNode):Node|null=>{
   const children=item.children.map(convert).filter((child):child is Node=>!!child);
   const own=item.directParts.filter(part=>remaining.delete(part.id));
   if(own.length)children.push(...own.map(part=>node([part])));
   const parts=[...new Map(children.flatMap(child=>child.parts).map(part=>[part.id,part])).values()];
   return parts.length?node(parts,children,item.name):null;
  };
  const categories=guestNodes.map(convert).filter((item):item is Node=>!!item);
  const other=atlas.parts.filter(part=>remaining.has(part.id));
  if(other.length)categories.push(node(other,other.map(part=>node([part])),'Other anatomy'));
  return node(visibleEntries.flatMap(entry=>entry.parts),categories);
 }
 if(mode==='systems'){
  const categories=SYSTEMS.map(system=>{
   const members=visibleEntries.filter(entry=>entry.system===system.id);
   return members.length?node(members.flatMap(entry=>entry.parts),buildAnatomyNodes(members,entry=>atlas.scope==='cell'?[]:[entry.region,...entry.location]).map(anatomyNode),system.name):null;
  }).filter((item):item is Node=>!!item);
  return node(visibleEntries.flatMap(entry=>entry.parts),categories);
 }
 if(mode==='regions'){
  const names=atlas.scope==='cell'?['Cell boundary','Nucleus','Cytoplasm']:REGION_ORDER;
  const categories=names.map(name=>{
   const members=visibleEntries.filter(entry=>entry.region===name);
   return members.length?node(members.flatMap(entry=>entry.parts),buildAnatomyNodes(members,entry=>[...entry.location,systemName(entry.system)]).map(anatomyNode),name):null;
  }).filter((item):item is Node=>!!item);
  return node(visibleEntries.flatMap(entry=>entry.parts),categories);
 }
 const categories=DEPTH_LAYERS.map(layer=>{
  const members=visibleEntries.flatMap(entry=>{const parts=entry.parts.filter(part=>depthLayerFor(part)===layer.id);return parts.length?[{...entry,parts}]:[]});
  return members.length?node(members.flatMap(entry=>entry.parts),buildAnatomyNodes(members,entry=>[entry.region,...entry.location.filter(name=>name.toLowerCase()!==layer.name.toLowerCase())]).map(anatomyNode),layer.name):null;
 }).filter((item):item is Node=>!!item);
 return node(visibleEntries.flatMap(entry=>entry.parts),categories);
}

function pack(children:Measured[],aspect:number,gap:number){
 const area=children.reduce((total,child)=>total+(child.width+gap)*(child.height+gap),0);
 const maxWidth=Math.max(...children.map(child=>child.width+gap));
 const targetWidth=Math.max(maxWidth,Math.sqrt(area*Math.max(.45,Math.min(1.8,aspect)))*1.14);
 let x=0,y=0,rowHeight=0,usedWidth=0;
 const cells=children.map(child=>{
  const width=child.width+gap,height=child.height+gap;
  if(x>0&&x+width>targetWidth){x=0;y+=rowHeight;rowHeight=0;}
  const cell={measured:child,x:x+width/2,y:-y-height/2};
  x+=width;rowHeight=Math.max(rowHeight,height);usedWidth=Math.max(usedWidth,x);
  return cell;
 });
 const usedHeight=y+rowHeight;
 for(const cell of cells){cell.x-=usedWidth/2;cell.y+=usedHeight/2;}
 return {cells,width:usedWidth,height:usedHeight};
}

/** Each slider interval separates one more tree level. Unopened descendants retain their assembled anatomy. */
export function createHierarchicalExplosionLayout(atlas:Atlas,entries:AnatomyEntry[],visibleIds:Set<string>,mode:ExplodeHierarchy,aspect=1,guestNodes?:ResolvedGuestNode[]):HierarchicalExplosionLayout{
 const root=hierarchyTree(atlas,entries,visibleIds,mode,guestNodes),index=new Map(atlas.parts.map((part,i)=>[part.id,i]));
 const initial=new Float32Array(atlas.parts.length*3);
 atlas.parts.forEach((part,i)=>{initial[i*3]=partCenter(part,0);initial[i*3+1]=partCenter(part,1);initial[i*3+2]=partCenter(part,2);});
 if(!root)return {stages:[{positions:initial,width:0,height:0,groups:[],clusterIds:new Int32Array(atlas.parts.length).fill(-1),clusters:[]}],steps:0};
 const gap=Math.max(.01,root.bounds.height*.018);
 const gapAt=(depth:number)=>gap*Math.pow(.7,depth);
 const maxDepth=(current:Node,depth:number):number=>current.children.length>1?Math.max(depth+1,...current.children.map(child=>maxDepth(child,depth+1))):current.children.length?maxDepth(current.children[0],depth+1):0;
 const steps=Math.max(1,maxDepth(root,0));
 const extent=(positions:Float32Array)=>{
  let minX=Infinity,maxX=-Infinity,minY=Infinity,maxY=-Infinity;
  for(const part of root.parts){const i=index.get(part.id)!;minX=Math.min(minX,positions[i*3]-(part.bounds[1][0]-part.bounds[0][0])/2);maxX=Math.max(maxX,positions[i*3]+(part.bounds[1][0]-part.bounds[0][0])/2);minY=Math.min(minY,positions[i*3+1]-(part.bounds[1][1]-part.bounds[0][1])/2);maxY=Math.max(maxY,positions[i*3+1]+(part.bounds[1][1]-part.bounds[0][1])/2);}
  return {width:maxX-minX,height:maxY-minY};
 };
 const stageFor=(level:number):ExplosionStage=>{
  const measure=(current:Node,depth:number):Measured=>{
   if(depth>=level||!current.children.length)return {node:current,width:Math.max(gapAt(depth),current.bounds.width),height:Math.max(gapAt(depth),current.bounds.height),children:null};
   const children=current.children.map(child=>measure(child,depth+1));
   const packed=pack(children,aspect,gapAt(depth));
   return {node:current,width:packed.width,height:packed.height,children:packed.cells};
  };
  const positions=initial.slice(),groups:ExplosionGroup[]=[],measured=measure(root,0),clusterIds=new Int32Array(atlas.parts.length).fill(-1),clusters:number[][]=[];
  const recordClusters=(item:Measured)=>{if(item.children){item.children.forEach(child=>recordClusters(child.measured));return;}const members=item.node.parts.map(part=>index.get(part.id)!);const cluster=clusters.length;clusters.push(members);for(const member of members)clusterIds[member]=cluster;};
  recordClusters(measured);
  const recordCollapsed=(current:Node,x:number,y:number,z:number,depth:number,id:string)=>{
   if(current.name&&current.children.length)groups.push({id,name:current.name,depth,x,y,z,height:current.bounds.height,count:current.parts.length});
   current.children.forEach((child,index)=>recordCollapsed(child,x+child.bounds.cx-current.bounds.cx,y+child.bounds.cy-current.bounds.cy,z+child.bounds.cz-current.bounds.cz,depth+1,`${id}/${index}`));
  };
  const place=(item:Measured,x:number,y:number,z:number,depth:number,id:string)=>{
   if(item.node.name&&item.node.children.length)groups.push({id,name:item.node.name,depth,x,y,z,height:item.height,count:item.node.parts.length});
   if(item.children){for(let index=0;index<item.children.length;index++){
    const child=item.children[index];
    place(child.measured,x+child.x,y+child.y,0,depth+1,`${id}/${index}`);
   }return;}
   const origin=item.node.bounds;
   for(const part of item.node.parts){const i=index.get(part.id)!;positions[i*3]=x+partCenter(part,0)-origin.cx;positions[i*3+1]=y+partCenter(part,1)-origin.cy;positions[i*3+2]=z+partCenter(part,2)-origin.cz;}
   item.node.children.forEach((child,index)=>recordCollapsed(child,x+child.bounds.cx-origin.cx,y+child.bounds.cy-origin.cy,z+child.bounds.cz-origin.cz,depth+1,`${id}/${index}`));
  };
  place(measured,0,.85,0,0,'root');
  return {positions,...extent(positions),groups,clusterIds,clusters};
 };
 const initialGroups:ExplosionGroup[]=[];
 const recordInitial=(current:Node,depth:number,id:string)=>{
  if(current.name&&current.children.length)initialGroups.push({id,name:current.name,depth,x:current.bounds.cx,y:current.bounds.cy,z:current.bounds.cz,height:current.bounds.height,count:current.parts.length});
  current.children.forEach((child,index)=>recordInitial(child,depth+1,`${id}/${index}`));
 };
 recordInitial(root,0,'root');
 const initialClusterIds=new Int32Array(atlas.parts.length).fill(-1),initialMembers=root.parts.map(part=>index.get(part.id)!);for(const member of initialMembers)initialClusterIds[member]=0;
 const stages:ExplosionStage[]=[{positions:initial,...extent(initial),groups:initialGroups,clusterIds:initialClusterIds,clusters:[initialMembers]}];
 for(let level=1;level<=steps;level++)stages.push(stageFor(level));
 return {stages,steps};
}
