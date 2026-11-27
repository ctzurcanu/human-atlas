import * as T from 'three';
import type {Part} from './anatomy';
import {organGeometryPoint} from './organ-frame';
type Point=[number,number,number];
type SurfaceSelection={conceptId:string;representation?:string;surface?:{partId:string;faces:number[]};anchor?:{partId:string;point:Point};curves?:{partId:string;points:Point[];closed:boolean}[]};
export const surfaceFeatureCarrier=(feature:SurfaceSelection|undefined)=>feature?.surface?.partId??feature?.anchor?.partId??feature?.curves?.[0]?.partId;

/** Selection ink on existing facets. These display-only overlays are never
 * anatomical meshes, section counters, depth layers or explosion members. */
export function createSurfaceFeatureOverlay(parent:T.Group,planes:T.Plane[]){
 const fillMaterial=new T.MeshBasicMaterial({color:0x77d1db,opacity:.3,transparent:true,side:T.DoubleSide,depthWrite:false,polygonOffset:true,polygonOffsetFactor:-3,polygonOffsetUnits:-3,clippingPlanes:planes});
 const edgeMaterial=new T.LineBasicMaterial({color:0x9ae8ee,depthWrite:false,transparent:true,opacity:.95,clippingPlanes:planes});
 const pointMaterial=new T.PointsMaterial({color:0x9ae8ee,size:5,sizeAttenuation:false,depthWrite:false,clippingPlanes:planes});
 const fill=new T.Mesh(new T.BufferGeometry(),fillMaterial),outline=new T.LineSegments(new T.BufferGeometry(),edgeMaterial);
 const marker=new T.Points(new T.BufferGeometry(),pointMaterial);
 for(const object of [fill,outline,marker]){object.matrixAutoUpdate=false;object.frustumCulled=false;object.visible=false;object.renderOrder=3;parent.add(object);}
 let key='',source:T.BufferGeometry|undefined;
 const clear=()=>{fill.geometry.deleteAttribute('position');fill.geometry.dispose();outline.geometry.dispose();marker.geometry.dispose();};
 const update=(feature:SurfaceSelection|undefined,host:T.Mesh|undefined,carrier?:Part)=>{
  fill.visible=outline.visible=marker.visible=false;if(!feature||!host)return;
  fill.visible=!!feature.surface;outline.visible=!!feature.surface||!!feature.curves?.length;
  marker.visible=feature.representation==='point-landmark'&&!!feature.anchor;
  fillMaterial.opacity=feature.conceptId.startsWith('TA98HEART:') ? .6 : .3;
  fillMaterial.depthTest=edgeMaterial.depthTest=!feature.conceptId.includes(':label-site:');
  if(key!==feature.conceptId||source!==host.geometry){
   clear();key=feature.conceptId;source=host.geometry;
   const p=source.getAttribute('position'),index=source.getIndex()!,chosen=new Uint32Array((feature.surface?.faces.length??0)*3);
   const edges=new Map<string,{a:number;b:number;uses:number}>();
   const vertexKey=(i:number)=>`${p.getX(i)},${p.getY(i)},${p.getZ(i)}`;
   feature.surface?.faces.forEach((face,j)=>{
    const ids=[index.getX(face*3),index.getX(face*3+1),index.getX(face*3+2)];chosen.set(ids,j*3);
    for(let i=0;i<3;i++){const a=ids[i],b=ids[(i+1)%3],key=[vertexKey(a),vertexKey(b)].sort().join('|'),edge=edges.get(key);if(edge)edge.uses++;else edges.set(key,{a,b,uses:1});}
   });
   fill.geometry=new T.BufferGeometry().setAttribute('position',p).setIndex(new T.BufferAttribute(chosen,1));
   const lines:number[]=[];for(const edge of edges.values())if(edge.uses===1)for(const i of [edge.a,edge.b])lines.push(p.getX(i),p.getY(i),p.getZ(i));
   const local=(point:Point)=>carrier?organGeometryPoint(carrier,new T.Vector3(...point)).toArray():point;
   for(const curve of feature.curves??[])for(let i=0;i<curve.points.length-(curve.closed?0:1);i++)for(const point of [curve.points[i],curve.points[(i+1)%curve.points.length]])lines.push(...local(point));
   outline.geometry=new T.BufferGeometry().setAttribute('position',new T.Float32BufferAttribute(lines,3));
   marker.geometry=new T.BufferGeometry().setAttribute('position',new T.Float32BufferAttribute(feature.anchor?local(feature.anchor.point):[],3));
  }
  for(const object of [fill,outline,marker]){object.matrix.copy(host.matrixWorld);object.matrixWorldNeedsUpdate=true;}
 };
 const dispose=()=>{clear();fill.removeFromParent();outline.removeFromParent();marker.removeFromParent();fillMaterial.dispose();edgeMaterial.dispose();pointMaterial.dispose();};
 return {update,dispose};
}
