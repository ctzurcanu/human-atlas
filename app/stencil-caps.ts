import {tissueWallMaterial} from './tissue-wall-material';
import {organBounds} from './organ-frame';
import * as T from 'three';
import type {Part} from './anatomy';
import {sectionTissue,type SectionTissue} from './section-tissue';
import {stencilShell,disposeStencilShell} from './stencil-section';
import {sectionCapProfile} from './section-cap-profile';
import {applyHepaticCut,HEPATIC_CUT_COLOR} from './hepatic-appearance';

export type StencilCapCandidate={part:Part;index:number;cutIndex:number;geometry:T.BufferGeometry;matrix:T.Matrix4;mode?:'wall'|'solid';wallWidth?:number;closedWall?:boolean;opacity?:number};
type Counter={outer:T.Mesh[];inner:T.Mesh[];geometry:T.BufferGeometry;source:T.BufferGeometry;assembly:boolean;innerPick?:T.BufferGeometry;role:SectionTissue;mode:'wall'|'solid'};
type Style={fill:number;rim?:number;inset?:number;tile?:string};
const styles:Partial<Record<SectionTissue,Style>>={
 'gastric-wall':{fill:0xA86A60},mucosa:{fill:0xC47782},submucosa:{fill:0xC99F8F},
 muscular:{fill:0xB0605A,rim:0x8A443C,inset:.0007,tile:'muscle'},
 cardiac:{fill:0xAE3239,rim:0x6E1219,inset:.0009,tile:'muscle'},
 skeletal:{fill:0xCBA286,rim:0xEFE9DC,inset:.0016,tile:'/models/cancellous-bone-seamless.png'},
 cartilage:{fill:0x86AEBD},
 tendon:{fill:0xC3B79B},
 adipose:{fill:0xD8BE91,tile:'/models/adipose-cut-v2.png'},
 lung:{fill:0xD99DA3,rim:0xAE777B,inset:.0006,tile:'/models/lung-cut-seamless.png'},
 liver:{fill:HEPATIC_CUT_COLOR},
 spleen:{fill:0xCBA85E,rim:0x98753F,inset:.0006,tile:'/models/spleen-cut-seamless.png'},
 cns:{fill:0xC4BBB2,rim:0x958B82,inset:.0008,tile:'/models/cns-cut-seamless.png'},
 organ:{fill:0xB98580,rim:0x8A615D,inset:.0006,tile:'/models/organ-cut-seamless.png'},
 nervous:{fill:0xE7C55C,rim:0x9E7619,inset:.00024},
 lymphatic:{fill:0x7FB894},
 connective:{fill:0xB7B4AB},
 attachments:{fill:0x8E5049},
 integumentary:{fill:0xBBAA9F},
 fascia:{fill:0xBEA99C},
 serosa:{fill:0xBC8B82},
 arterial:{fill:0xC8393F},
 venous:{fill:0x466FBA},
};
const roleOrder:SectionTissue[]=['integumentary','fascia','serosa','submucosa','mucosa','gastric-wall','adipose','tendon','muscular','cardiac','skeletal','cartilage','lung','liver','organ','spleen','lymphatic','cns','nervous','arterial','venous','connective','attachments'];

/** GPU cut faces: repaired source shells count stencil coverage at the live plane. */
export function createStencilCaps(scene:T.Scene|T.Group,planes:[T.Plane,T.Plane],textureFor:(url:string)=>T.Texture){
 const counters=new Map<number,Counter>(),countMaterials:T.MeshBasicMaterial[]=[],quadMaterials:(T.MeshBasicMaterial|T.MeshStandardMaterial)[]=[],quads:T.Mesh[]=[];
 const planeGeometry=new T.PlaneGeometry(1,1),active=new Set<string>();
 let visibleCandidates:StencilCapCandidate[]=[];
 const pickMaterial=new T.MeshBasicMaterial({side:T.DoubleSide}),pickMesh:T.Mesh<T.BufferGeometry,T.MeshBasicMaterial>=new T.Mesh(planeGeometry,pickMaterial);
 const insideRay=new T.Raycaster(),pickBox=new T.Box3(),pickPoint=new T.Vector3(),pickOrigin=new T.Vector3();
 const axes=[0,1].map(()=>({u:new T.Vector3(1,0,0),v:new T.Vector3(0,1,0)}));
 const countMaterial=(cutIndex:number,side:T.Side,inside:boolean,reverse=false)=>{
  // Winding at this cut comes from its one-sided closed source. Clipping the
  // counters by another cut opens that volume and cancels real tissue coverage.
  const material=new T.MeshBasicMaterial({side,colorWrite:false,depthWrite:false,depthTest:false,clippingPlanes:[planes[cutIndex]]});
  material.stencilWrite=true;material.stencilFunc=T.AlwaysStencilFunc;
  const op=(side===T.BackSide)!==reverse?T.IncrementWrapStencilOp:T.DecrementWrapStencilOp;
  material.stencilFail=op;material.stencilZFail=op;material.stencilZPass=op;
  if(inside){material.onBeforeCompile=shader=>{shader.vertexShader='attribute vec3 capInset;\n'+shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\ntransformed-=capInset;');};material.customProgramCacheKey=()=>`stencil-section-inside-${cutIndex}-${side}`;}
  countMaterials.push(material);return material;
 };
 const counts=([0,1] as const).map(cut=>({back:countMaterial(cut,T.BackSide,false),front:countMaterial(cut,T.FrontSide,false),innerBack:countMaterial(cut,T.BackSide,true),innerFront:countMaterial(cut,T.FrontSide,true),wallBack:countMaterial(cut,T.BackSide,true,true),wallFront:countMaterial(cut,T.FrontSide,true,true)}));
 const quadMaterial=(role:SectionTissue,cutIndex:number,rim=false)=>{
  const style=styles[role]??{fill:0x8B9099},parameters={color:rim?style.rim??style.fill:style.tile?0xFFFFFF:style.fill,side:T.DoubleSide,depthTest:true,depthWrite:true,polygonOffset:true,polygonOffsetFactor:-1,polygonOffsetUnits:-1,clippingPlanes:[planes[1-cutIndex]]};
  const material=role==='gastric-wall'?new T.MeshStandardMaterial({...parameters,roughness:.86}):new T.MeshBasicMaterial(parameters);
  material.stencilWrite=true;material.stencilRef=0;material.stencilFunc=T.NotEqualStencilFunc;
  material.stencilFail=T.ReplaceStencilOp;material.stencilZFail=T.ReplaceStencilOp;material.stencilZPass=T.ReplaceStencilOp;
  if(!rim){
   material.onBeforeCompile=shader=>{
    shader.uniforms.capAxisU={value:axes[cutIndex].u};shader.uniforms.capAxisV={value:axes[cutIndex].v};shader.uniforms.capNormal={value:planes[cutIndex].normal};
    if(style.tile==='muscle'){
     shader.uniforms.capAxial={value:textureFor('/models/muscle-axial-seamless.png')};
     shader.uniforms.capCoronal={value:textureFor('/models/muscle-coronal-seamless.png')};
     shader.uniforms.capSagittal={value:textureFor('/models/muscle-sagittal-seamless.png')};
    }else if(style.tile)shader.uniforms.capTile={value:textureFor(style.tile)};
    shader.vertexShader='varying vec3 capWorld;\n'+shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\ncapWorld=(modelMatrix*vec4(transformed,1.0)).xyz;');
    const uniforms=style.tile==='muscle'?'uniform sampler2D capAxial,capCoronal,capSagittal;':style.tile?'uniform sampler2D capTile;':'';
    const sample=style.tile==='muscle'?`vec3 w=pow(abs(normalize(capNormal)),vec3(4.0));w/=max(w.x+w.y+w.z,1e-5);diffuseColor.rgb=texture2D(capSagittal,cutUv*25.0).rgb*w.x+texture2D(capAxial,cutUv*25.0).rgb*w.y+texture2D(capCoronal,cutUv*25.0).rgb*w.z;`:style.tile?'diffuseColor.rgb=texture2D(capTile,cutUv*24.0).rgb;':['gastric-wall','mucosa','submucosa'].includes(role)?'diffuseColor.rgb*=.96+.04*sin(cutUv.x*1700.0)*sin(cutUv.y*1300.0);':'';
    shader.fragmentShader=`varying vec3 capWorld;uniform vec3 capAxisU,capAxisV,capNormal;${uniforms}\n`+shader.fragmentShader.replace('#include <color_fragment>',`#include <color_fragment>\nvec2 cutUv=vec2(dot(capWorld,capAxisU),dot(capWorld,capAxisV));${sample}`);
   };
   material.customProgramCacheKey=()=>`stencil-cap-${role}-${cutIndex}`;
   if(role==='liver')applyHepaticCut(material);
  }
  quadMaterials.push(material);return material;
 };
 const quadMap=new Map<string,{rim?:T.Mesh;fill:T.Mesh}>();
 // A fill clipped by the other plane cannot reset discarded stencil fragments
 // through its Replace ops. Clear after the final fill of each counting batch,
 // including depth-occluded pixels, before another tissue or cut is counted.
 const orderFor=(role:SectionTissue,cutIndex:number,mode:'wall'|'solid')=>-1000+roleOrder.indexOf(role)*16+cutIndex*8+(mode==='wall'?4:0);
 const quadFor=(role:SectionTissue,cutIndex:number,mode:'wall'|'solid')=>{
  const key=`${cutIndex}:${role}:${mode}`,existing=quadMap.get(key);if(existing)return existing;
  const style=styles[role]??{fill:0x8B9099},base=orderFor(role,cutIndex,mode);
  const make=(rim:boolean,order:number)=>{const mesh=new T.Mesh(planeGeometry,quadMaterial(role,cutIndex,rim));mesh.visible=false;mesh.frustumCulled=false;mesh.renderOrder=order;mesh.userData.cutIndex=cutIndex;if(!rim)mesh.onAfterRender=renderer=>renderer.clearStencil();scene.add(mesh);quads.push(mesh);return mesh;};
  const result={rim:mode==='solid'&&style.inset?make(true,base+1):undefined,fill:make(false,base+(mode==='wall'?2:style.inset?3:1))};quadMap.set(key,result);return result;
 };
 const makeCounter=(candidate:StencilCapCandidate)=>{
  const {index,geometry,part}=candidate,role=sectionTissue(part),style=styles[role]??{fill:0x8B9099},mode=candidate.mode??'solid',shell=stencilShell(geometry,mode==='wall'?candidate.wallWidth??.002:candidate.part.pairedTissue?0:style.inset??0,!!candidate.closedWall,mode==='wall'&&!part.sectionAssembly&&!!sectionCapProfile(part).hollowWall);
  const entry:Counter={geometry:shell,source:geometry,assembly:!!part.sectionAssembly,role,mode,outer:[],inner:[]};
  if(mode==='wall'&&!shell.userData.openPortWall){
   const positions=shell.getAttribute('position'),insets=shell.getAttribute('capInset'),values=new Float32Array(positions.count*3);
   for(let i=0;i<positions.count;i++)for(let axis=0;axis<3;axis++)values[i*3+axis]=positions.getComponent(i,axis)-insets.getComponent(i,axis);
   entry.innerPick=new T.BufferGeometry();entry.innerPick.setAttribute('position',new T.Float32BufferAttribute(values,3));entry.innerPick.setIndex(shell.getIndex());entry.innerPick.userData.borrowedStencilIndex=true;
   entry.innerPick.boundingBox=shell.boundingBox?.clone()??null;entry.innerPick.boundingSphere=shell.boundingSphere?.clone()??null;
  }
  for(const cutIndex of [0,1]){
   const base=part.pairedTissue?-100000+index*8+cutIndex*4:orderFor(role,cutIndex,mode),mat=counts[cutIndex];
   const make=(material:T.MeshBasicMaterial,order:number)=>{const mesh=new T.Mesh(shell,material);mesh.visible=false;mesh.frustumCulled=false;mesh.matrixAutoUpdate=false;mesh.renderOrder=order;scene.add(mesh);return mesh;};
   entry.outer.push(make(mat.back,base),make(mat.front,base));
   if(mode==='wall'&&!shell.userData.openPortWall)entry.inner.push(make(mat.wallBack,base+1),make(mat.wallFront,base+1));
   else if(mode==='solid'&&style.inset&&!part.pairedTissue)entry.inner.push(make(mat.innerBack,base+2),make(mat.innerFront,base+2));
  }
  counters.set(index,entry);return entry;
 };
 const transform=(mesh:T.Mesh,matrix:T.Matrix4)=>{mesh.matrix.copy(matrix);mesh.matrixWorldNeedsUpdate=true;};
 const update=(candidates:StencilCapCandidate[])=>{
  const present=new Set(candidates.map(candidate=>candidate.index));
  release([...counters].filter(([index,entry])=>entry.assembly&&!present.has(index)).map(([index])=>index));
  for(const candidate of candidates){const entry=counters.get(candidate.index);if(entry&&entry.source!==candidate.geometry)release([candidate.index]);}
  visibleCandidates=candidates;
  active.clear();for(const entry of counters.values())for(const mesh of [...entry.outer,...entry.inner])mesh.visible=false;
  for(const mesh of quads)mesh.visible=false;
  const byCut:[StencilCapCandidate[],StencilCapCandidate[]]=[[],[]];
  for(const candidate of candidates){const entry=counters.get(candidate.index)??makeCounter(candidate),cut=candidate.cutIndex;
   for(const mesh of entry.outer.slice(cut*2,cut*2+2)){mesh.visible=true;transform(mesh,candidate.matrix);}
   for(const mesh of entry.inner.slice(cut*2,cut*2+2)){mesh.visible=true;transform(mesh,candidate.matrix);}
   active.add(`${cut}:${candidate.index}`);byCut[cut].push(candidate);
  }
  for(const cutIndex of [0,1]){
   const list=byCut[cutIndex];if(!list.length)continue;
   const plane=planes[cutIndex],normal=plane.normal,rotation=new T.Quaternion().setFromUnitVectors(new T.Vector3(0,0,1),normal.clone().negate());
   const {u,v}=axes[cutIndex];u.set(1,0,0).applyQuaternion(rotation);v.set(0,1,0).applyQuaternion(rotation);
   const origin=normal.clone().multiplyScalar(-plane.constant),low=[Infinity,Infinity],high=[-Infinity,-Infinity],corner=new T.Vector3();
   for(const {part,matrix} of list){for(let x=0;x<2;x++)for(let y=0;y<2;y++)for(let z=0;z<2;z++){
    const bounds=organBounds(part,matrix);corner.set(x?bounds.max.x:bounds.min.x,y?bounds.max.y:bounds.min.y,z?bounds.max.z:bounds.min.z);
    corner.addScaledVector(normal,-plane.distanceToPoint(corner));
    const a=corner.dot(u),b=corner.dot(v);low[0]=Math.min(low[0],a);high[0]=Math.max(high[0],a);low[1]=Math.min(low[1],b);high[1]=Math.max(high[1],b);
   }}
   const center=origin.addScaledVector(u,(low[0]+high[0])*.5).addScaledVector(v,(low[1]+high[1])*.5);
   const width=(high[0]-low[0])*1.06+.0001,height=(high[1]-low[1])*1.06+.0001;
   for(const candidate of list.filter(item=>item.part.pairedTissue)){
    const key=`tissue:${cutIndex}:${candidate.index}`;let pair=quadMap.get(key);
    if(!pair){const appearance=candidate.part.tissueAppearance;if(!appearance)throw new Error('Paired tissue requires a generated appearance');const tissue=tissueWallMaterial(appearance,candidate.matrix.clone().invert()),material=tissue.material;material.side=T.DoubleSide;material.clippingPlanes=[planes[1-cutIndex]];material.stencilWrite=true;material.stencilFunc=T.NotEqualStencilFunc;material.stencilRef=0;material.stencilFail=T.ReplaceStencilOp;material.stencilZFail=T.ReplaceStencilOp;material.stencilZPass=T.ReplaceStencilOp;material.polygonOffset=true;material.polygonOffsetFactor=-1;material.polygonOffsetUnits=-1;quadMaterials.push(material);const mesh=new T.Mesh(planeGeometry,material);mesh.userData.tissueFrame=tissue.tissueFrame;mesh.userData.cutIndex=cutIndex;mesh.frustumCulled=false;mesh.renderOrder=-100000+candidate.index*8+cutIndex*4+1;mesh.onAfterRender=renderer=>renderer.clearStencil();scene.add(mesh);quads.push(mesh);pair={fill:mesh};quadMap.set(key,pair);}
    const capMaterial=pair.fill.material as T.MeshStandardMaterial,opacity=Math.max(0,Math.min(1,candidate.opacity??1));
    // Count/fill/reset must stay in the opaque render queue even when the
    // anatomical context is translucent. Custom alpha blending preserves that
    // exact sequence; transparent sorting would mix unrelated tissue stencils.
    capMaterial.opacity=opacity;capMaterial.depthWrite=opacity>=.999;capMaterial.blending=opacity>=.999?T.NoBlending:T.CustomBlending;capMaterial.blendSrc=T.SrcAlphaFactor;capMaterial.blendDst=T.OneMinusSrcAlphaFactor;capMaterial.blendEquation=T.AddEquation;
    pair.fill.userData.tissueFrame.value.copy(candidate.matrix).invert();pair.fill.visible=true;pair.fill.position.copy(center);pair.fill.quaternion.copy(rotation);pair.fill.scale.set(width,height,1);
   }
   for(const key of new Set(list.filter(candidate=>!candidate.part.pairedTissue).map(candidate=>`${sectionTissue(candidate.part)}:${candidate.mode??'solid'}`))){
    const [role,mode]=key.split(':') as [SectionTissue,'wall'|'solid'];
    const pair=quadFor(role,cutIndex,mode);for(const mesh of [pair.rim,pair.fill])if(mesh){mesh.visible=true;mesh.position.copy(center);mesh.quaternion.copy(rotation);mesh.scale.set(width,height,1);}
   }
  }
 };
 const pick=(ray:T.Ray,otherVisible:(index:number)=>boolean)=>{
  let result:{index:number;distance:number;order:number}|undefined;
  for(const candidate of visibleCandidates){
   if(!otherVisible(candidate.index))continue;
   const plane=planes[candidate.cutIndex];
   if(!ray.intersectPlane(plane,pickPoint)||planes[1-candidate.cutIndex].distanceToPoint(pickPoint)<-.00001)continue;
   const distance=ray.origin.distanceTo(pickPoint);
   if(result&&distance>result.distance+.00001)continue;
   pickBox.copy(organBounds(candidate.part,candidate.matrix));
   if(!pickBox.expandByScalar(.00002).containsPoint(pickPoint))continue;
   const entry=counters.get(candidate.index);if(!entry)continue;
   pickMesh.matrix.copy(candidate.matrix);pickMesh.matrixAutoUpdate=false;pickMesh.updateMatrixWorld(true);
   pickOrigin.copy(pickPoint).addScaledVector(plane.normal,.00002);
   insideRay.set(pickOrigin,plane.normal);insideRay.near=0;insideRay.far=Infinity;
   const inside=(geometry:T.BufferGeometry)=>{pickMesh.geometry=geometry;const intersections=insideRay.intersectObject(pickMesh,false);let crossings=0,last=-Infinity;for(const hit of intersections){if(hit.distance-last<.00002)continue;crossings++;last=hit.distance;}return crossings%2===1;};
   if(!inside(entry.geometry)||entry.innerPick&&inside(entry.innerPick))continue;
   const order=orderFor(entry.role,candidate.cutIndex,entry.mode);
   if(!result||distance<result.distance-.00001||order>result.order)result={index:candidate.index,distance,order};
  }
  return result;
 };
 const dispose=()=>{for(const entry of counters.values()){for(const mesh of [...entry.outer,...entry.inner])scene.remove(mesh);disposeStencilShell(entry.geometry);if(entry.innerPick)disposeStencilShell(entry.innerPick);}for(const mesh of quads)scene.remove(mesh);for(const material of [...countMaterials,...quadMaterials,pickMaterial])material.dispose();planeGeometry.dispose();};
 const release=(indices:number[])=>{
  const removed=new Set(indices);visibleCandidates=visibleCandidates.filter(candidate=>!removed.has(candidate.index));
  for(const index of indices){const entry=counters.get(index);if(!entry)continue;for(const mesh of [...entry.outer,...entry.inner])scene.remove(mesh);disposeStencilShell(entry.geometry);if(entry.innerPick)disposeStencilShell(entry.innerPick);counters.delete(index);for(const cut of [0,1])active.delete(`${cut}:${index}`);}
 };
 return {update,pick,dispose,release,planeError:()=>Math.max(0,...quads.filter(mesh=>mesh.visible).map(mesh=>Math.abs(planes[mesh.userData.cutIndex].distanceToPoint(mesh.position)))),activeCount:()=>active.size};
}
