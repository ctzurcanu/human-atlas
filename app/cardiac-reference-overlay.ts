import * as T from 'three';
import type {Atlas,Part} from './anatomy';
import {cardiacLabelRegistry,type CardiacLabelSite} from './cardiac-label-sites';
import {organAnatomicalPoint,organGeometryPoint} from './organ-frame';
import {createSurfaceFeatureOverlay} from './surface-feature-overlay';

export function cardiacReferenceLocalPoint(carrier:Part,point:number[],coordinateFrame?:'atlas-world-metres',epicardium?:Part){
 const p=new T.Vector3().fromArray(point);
 return coordinateFrame==='atlas-world-metres'?organGeometryPoint(carrier,p):epicardium?organGeometryPoint(carrier,organAnatomicalPoint(epicardium,p)):undefined;
}

/** Typed reference display never changes anatomy, clipping, depth or the camera.
 * Indexed regions share the native position buffer. Empty spaces use wireframes. */
function createSingleCardiacReferenceOverlay(parent:T.Group,planes:T.Plane[],atlas:Atlas,beginVolumes:(renderer:T.WebGLRenderer)=>void){
 const parts=new Map(atlas.parts.map(p=>[p.id,p]));
 const epicardium=atlas.parts.find(p=>/^TA98HEART:(male|female):epicardium$/.test(p.id));
 let key='',objects:{object:T.Object3D;owner:string;geometry:T.BufferGeometry;material:T.Material}[]=[],surfaces:{overlay:ReturnType<typeof createSurfaceFeatureOverlay>;feature:{conceptId:string;surface:{partId:string;faces:number[]}}}[]=[];
 const clear=()=>{for(const row of objects){row.object.removeFromParent();row.geometry.dispose();row.material.dispose();}objects=[];for(const row of surfaces)row.overlay.dispose();surfaces=[];};
 const update=(feature:CardiacLabelSite|undefined,hostFor:(id:string)=>T.Mesh|undefined)=>{
  const host=feature&&hostFor(feature.anchor.partId),geometry=feature?.geometry;
  const next=feature&&geometry&&host?feature.conceptId+':'+String(host.geometry.id):'';
  if(next!==key){clear();key=next;
   if(feature&&geometry&&host&&(feature.coordinateFrame==='atlas-world-metres'||epicardium)){
    const owner=feature.anchor.partId,carrier=parts.get(owner)!;
    const local=(p:number[])=>cardiacReferenceLocalPoint(carrier,p,feature.coordinateFrame,epicardium)!;
    const add=(object:T.Object3D,g:T.BufferGeometry,m:T.Material)=>{object.matrixAutoUpdate=false;object.frustumCulled=false;object.renderOrder=4;parent.add(object);objects.push({object,owner,geometry:g,material:m});};
    if(geometry.type==='surface-region'){
     for(const s of geometry.surfaces)surfaces.push({overlay:createSurfaceFeatureOverlay(parent,planes),feature:{conceptId:feature.conceptId,surface:s}});
    }else if(geometry.type==='boundary-curve'||geometry.type==='opening-rim'){
     const points:number[]=[];
     for(const path of geometry.paths)for(let i=0;i<path.points.length-(path.closed?0:1);i++)for(const p of [path.points[i],path.points[(i+1)%path.points.length]])points.push(...local(p).toArray());
     const g=new T.BufferGeometry().setAttribute('position',new T.Float32BufferAttribute(points,3));
     const m=new T.LineBasicMaterial({color:0x8cdae3,transparent:true,opacity:1,depthTest:false,depthWrite:false,clippingPlanes:planes});add(new T.LineSegments(g,m),g,m);
    }else if(geometry.type==='point-set'){
     const points=geometry.points.flatMap(p=>local(p).toArray());
     const g=new T.BufferGeometry().setAttribute('position',new T.Float32BufferAttribute(points,3));
     const m=new T.PointsMaterial({color:0x9ae8ee,size:5,sizeAttenuation:false,depthTest:false,depthWrite:false,clippingPlanes:planes});add(new T.Points(g,m),g,m);
    }else if((geometry.type==='finite-mesh'||geometry.type==='space-region')&&'vertices'in geometry){
     const g=new T.BufferGeometry().setAttribute('position',new T.Float32BufferAttribute(geometry.vertices.flatMap(p=>local(p).toArray()),3)).setIndex(geometry.faces.flat());g.computeVertexNormals();
     const space=geometry.type==='space-region',color=space?0x8cdae3:feature.term.startsWith('A12.1.06.')?0xe4c66b:feature.term.startsWith('A12.1.08.')?0xd6c1b4:0xc3867f;
     // This is selection ink for a typed anatomical reference. Keep a tiny
     // intramural volume visible through its carrier without hiding tissue.
     const m=new T.MeshStandardMaterial({color,roughness:.65,side:space?T.DoubleSide:T.FrontSide,wireframe:space,transparent:true,opacity:space?.65:.94,depthTest:!space,depthWrite:!space,clippingPlanes:planes});
     const mesh=new T.Mesh(g,m);if(!space)mesh.onBeforeRender=beginVolumes;add(mesh,g,m);
    }
   }
  }
  for(const row of surfaces)row.overlay.update(row.feature,hostFor(row.feature.surface.partId));
  for(const row of objects){const owner=hostFor(row.owner);row.object.visible=!!owner;if(owner){row.object.matrix.copy(owner.matrixWorld);row.object.matrixWorldNeedsUpdate=true;}}
 };
 return {update,dispose:clear};
}

/** Collective names display their typed descendants; they never manufacture a
 * point or one generic solid in place of the constituent anatomy. */
export function createCardiacReferenceOverlay(parent:T.Group,planes:T.Plane[],atlas:Atlas){
 const sites=new Map(cardiacLabelRegistry(atlas)?.termGroups.map(g=>[g.term,g])??[]);
 const active=new Map<string,ReturnType<typeof createSingleCardiacReferenceOverlay>>();
 const volumeFrames=new WeakMap<T.WebGLRenderer,number>();
 const beginVolumes=(renderer:T.WebGLRenderer)=>{const frame=renderer.info.render.frame;if(volumeFrames.get(renderer)!==frame){renderer.clearDepth();volumeFrames.set(renderer,frame);}};
 const update=(feature:CardiacLabelSite|undefined,hostFor:(id:string)=>T.Mesh|undefined)=>{
  const visible=new Map<string,CardiacLabelSite>(),visited=new Set<string>();
  const collect=(g:CardiacLabelSite)=>{
   if(visited.has(g.term))return;visited.add(g.term);
   if(g.geometry?.type==='collection'){
    const terms=g.geometry.terms.filter(t=>!(t==='A12.1.08.006'&&g.geometry?.type==='collection'&&g.geometry.terms.includes('A12.1.08.008')&&g.geometry.terms.includes('A12.1.08.009')));
    for(const term of terms){const child=sites.get(term);if(child)collect(child);}
   }else if(g.geometry)visible.set(g.conceptId,g);
  };
  if(feature)collect(feature);
  for(const [id,overlay]of active)if(!visible.has(id)){overlay.dispose();active.delete(id);}
  for(const [id,g]of visible){let overlay=active.get(id);if(!overlay){overlay=createSingleCardiacReferenceOverlay(parent,planes,atlas,beginVolumes);active.set(id,overlay);}overlay.update(g,hostFor);}
 };
 return {update,dispose:()=>{for(const overlay of active.values())overlay.dispose();active.clear();}};
}
