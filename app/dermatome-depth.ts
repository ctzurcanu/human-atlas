import * as T from 'three';
import {GLTFLoader} from 'three/examples/jsm/loaders/GLTFLoader.js';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';
import {DEPTH_LAYERS,type DepthLayerId} from './depth-layers';
import {depthLayerOpacity,depthPosition} from './depth-control';
import type {SceneState} from './anatomy';

export type NativeDepthManifest={nativeBodyAlbedo:string;layers:{id:DepthLayerId;url:string;geometryBytes:number}[]};
type Loaded={group:T.Group;materials:Set<T.Material>;geometries:Set<T.BufferGeometry>;textures:Set<T.Texture>};
const release=(record:Loaded)=>{record.group.removeFromParent();record.geometries.forEach(g=>g.dispose());record.materials.forEach(m=>m.dispose());record.textures.forEach(t=>t.dispose());};
/** Same source coordinates as the skin; never registers this frame to atlas tissue. */
export function nativeDermatomeDepth(holder:T.Group,center:T.Vector3,manifest:NativeDepthManifest,changed:()=>void,failed:(message:string)=>void,clippingPlanes:T.Plane[]=[]){
 const loaded=new Map<string,Loaded>(),desired=new Set<string>(),loader=new GLTFLoader();let disposed=false,busy=false,timer:ReturnType<typeof setTimeout>|undefined,current:SceneState;
 const order=new Map(DEPTH_LAYERS.map((layer,index)=>[layer.id,index]));
 const layers=[...manifest.layers].sort((a,b)=>order.get(a.id)!-order.get(b.id)!);
 const opacity=(id:DepthLayerId)=>current.depthHidden?.includes(id)?0:depthLayerOpacity(id,current);
 const apply=()=>{for(const [id,record] of loaded){const alpha=opacity(id as DepthLayerId);record.group.visible=alpha>0;for(const original of record.materials){const material=original as T.MeshStandardMaterial;const transparent=alpha<.999;if(material.transparent!==transparent){material.transparent=transparent;material.needsUpdate=true;}material.opacity=alpha;material.depthWrite=!transparent;}}changed();};
 const batch=(source:T.Group):Loaded=>{
  source.updateMatrixWorld(true);const group=new T.Group(),materials=new Set<T.Material>(),geometries=new Set<T.BufferGeometry>(),textures=new Set<T.Texture>();
  const buckets=new Map<string,{material:T.Material;geometry:T.BufferGeometry[]}>(),retained:T.Mesh[]=[];
  source.traverse(object=>{if(!(object instanceof T.Mesh))return;
   const ms=Array.isArray(object.material)?object.material:[object.material];for(const m of ms){m.clippingPlanes=clippingPlanes;materials.add(m);for(const value of Object.values(m))if(value instanceof T.Texture)textures.add(value);}
   if(ms.length!==1||Object.keys(object.geometry.morphAttributes).length){retained.push(object);return;}
   const geometry=object.geometry.clone().applyMatrix4(object.matrixWorld),key=ms[0].uuid+':'+Object.entries(geometry.attributes as Record<string,T.BufferAttribute>).map(([k,a])=>`${k}:${a.itemSize}:${a.normalized}`).sort().join(',')+':'+!!geometry.index;
   const bucket=buckets.get(key)??{material:ms[0],geometry:[] as T.BufferGeometry[]};bucket.geometry.push(geometry);buckets.set(key,bucket);
  });
  for(const bucket of buckets.values()){
   const merged=bucket.geometry.length===1?bucket.geometry[0]:mergeGeometries(bucket.geometry,false);
   if(merged){if(bucket.geometry.length>1)bucket.geometry.forEach(g=>g.dispose());geometries.add(merged);group.add(new T.Mesh(merged,bucket.material));}
   else for(const geometry of bucket.geometry){geometries.add(geometry);group.add(new T.Mesh(geometry,bucket.material));}
  }
  for(const mesh of retained){group.attach(mesh);geometries.add(mesh.geometry);}
  source.traverse(o=>{if(o instanceof T.Mesh&&!geometries.has(o.geometry))o.geometry.dispose();});
  group.position.sub(center);return {group,materials,geometries,textures};
 };
 const pump=async()=>{if(disposed||busy)return;const next=layers.find(layer=>desired.has(layer.id)&&!loaded.has(layer.id));if(!next)return;busy=true;
  try{const gltf=await loader.loadAsync(next.url);const record=batch(gltf.scene);
   if(disposed||!desired.has(next.id)){release(record);}else{loaded.set(next.id,record);holder.add(record.group);apply();}
  }catch(error){if(!disposed){desired.delete(next.id);failed(`Could not load the native ${next.id} depth layer: ${String(error)}`);}}
  finally{busy=false;if(!disposed)void pump();}
 };
 return {
  update(state:SceneState){current=state;const progress=depthPosition(state)*DEPTH_LAYERS.length,skinAlpha=state.depthHidden?.includes('skin')?0:depthLayerOpacity('skin',state);
   const next=new Set<string>();let bytes=0;
   if(skinAlpha<.999)for(const layer of layers){if(order.get(layer.id)!<Math.floor(progress)||state.depthHidden?.includes(layer.id))continue;if(next.size>=6||next.size>1&&bytes+layer.geometryBytes>40*1024*1024)break;next.add(layer.id);bytes+=layer.geometryBytes;}
   const different=next.size!==desired.size||[...next].some(id=>!desired.has(id));
   if(different){desired.clear();next.forEach(id=>desired.add(id));for(const [id,record] of loaded)if(!desired.has(id)){release(record);loaded.delete(id);}if(timer)clearTimeout(timer);timer=setTimeout(()=>void pump(),150);}
   apply();
  },
  stats(){return {layers:[...loaded.keys()],geometryBytes:layers.filter(layer=>loaded.has(layer.id)).reduce((sum,layer)=>sum+layer.geometryBytes,0),drawMeshes:[...loaded.values()].reduce((sum,record)=>sum+record.group.children.length,0)};},
  dispose(){disposed=true;if(timer)clearTimeout(timer);loaded.forEach(release);loaded.clear();},
 };
}
