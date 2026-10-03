import * as T from 'three';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';
import {assetUrl} from './asset-url';
import {loadModelBuffer} from './model-download';
import {fetchAssetJson} from './asset-cache';
import {partLayerOpacity} from './depth-control';
import {SYSTEMS,type Atlas,type SceneState} from './anatomy';

interface VrPart {id:string;sourceVertexCount:number;sourceIndexCount:number;positions:number;normals:number;indices:number;vertexCount:number;indexCount:number}
interface VrModel {version:number;sourceTriangles:number;triangles:number;parts:VrPart[];chunk:{url:string;bytes:number;gzip:boolean}}
interface Uniforms {rotationState:T.DataTexture;rotationPivots:T.DataTexture;rotation:{value:T.Vector4};pivot:{value:T.Vector3};planes:T.Plane[]}

/** A single lightweight draw per eye; hidden parts never reach the GPU. */
export async function loadVrAnatomy(atlas:Atlas,url:string,uniforms:Uniforms,signal:AbortSignal){
 const model=await fetchAssetJson<VrModel>(assetUrl(url),{signal},'The lightweight VR model is unavailable. Reload the viewer or choose another model.');
 if(model.version!==1||model.parts.length!==atlas.parts.length||model.parts.some((part,i)=>part.id!==atlas.parts[i].id||part.sourceVertexCount!==atlas.parts[i].vertexCount||part.sourceIndexCount!==atlas.parts[i].indexCount))throw new Error('The lightweight VR model does not match this anatomy. Reload the viewer or choose another model.');
 const buffer=await loadModelBuffer(assetUrl(model.chunk.url),model.chunk.bytes,model.chunk.gzip,signal);
 if(signal.aborted)throw new DOMException('Aborted','AbortError');
 const width=T.MathUtils.ceilPowerOfTwo(atlas.parts.length),states=new Float32Array(width*4),texture=new T.DataTexture(states,width,1,T.RGBAFormat,T.FloatType);
 const pickMaterial=new T.MeshBasicMaterial({side:T.DoubleSide});
 const pickers=model.parts.map((part,i)=>{
  const source=atlas.parts[i],spec=atlas.materials?.[source.material??''],geometry=new T.BufferGeometry();
  geometry.setAttribute('position',new T.BufferAttribute(new Float32Array(buffer,part.positions,part.vertexCount*3),3));
  geometry.setAttribute('normal',new T.BufferAttribute(new Int16Array(buffer,part.normals,part.vertexCount*3),3,true));
  geometry.setIndex(new T.BufferAttribute(new Uint32Array(buffer,part.indices,part.indexCount),1));
  geometry.setAttribute('partIndex',new T.BufferAttribute(new Float32Array(part.vertexCount).fill(i),1));
  const color=spec?.color?new T.Color().setRGB(spec.color[0]/255,spec.color[1]/255,spec.color[2]/255,T.SRGBColorSpace):new T.Color(source.tissue==='tendon'?'#d6c8b0':source.tissue==='cartilage'?'#afc0cb':SYSTEMS.find(system=>system.id===source.system)?.color??'#aebbb8');
  if(/^spleen(?:\s|$)/i.test(source.name))color.setRGB(123/255,65/255,84/255,T.SRGBColorSpace);
  const colors=new Float32Array(part.vertexCount*3);
  for(let v=0;v<part.vertexCount;v++)color.toArray(colors,v*3);
  geometry.setAttribute('color',new T.BufferAttribute(colors,3));geometry.computeBoundingBox();geometry.computeBoundingSphere();
  const picker=new T.Mesh(geometry,pickMaterial);picker.matrixAutoUpdate=false;return picker;
 });
 const geometry=mergeGeometries(pickers.map(picker=>picker.geometry));
 if(!geometry)throw new Error('The lightweight VR geometry could not be assembled.');
 const originalIndices=Uint32Array.from(geometry.getIndex()!.array),activeIndices=new Uint32Array(originalIndices.length),shown=new Uint8Array(model.parts.length);
 const indices=new T.BufferAttribute(activeIndices,1).setUsage(T.DynamicDrawUsage);geometry.setIndex(indices);geometry.setDrawRange(0,0);
 const material=new T.MeshLambertMaterial({vertexColors:true,side:T.DoubleSide,alphaHash:true,clippingPlanes:uniforms.planes});
 material.onBeforeCompile=shader=>{
  shader.uniforms.vrState={value:texture};shader.uniforms.stateWidth={value:width};shader.uniforms.rotationState={value:uniforms.rotationState};shader.uniforms.rotationPivots={value:uniforms.rotationPivots};shader.uniforms.selectionRotation=uniforms.rotation;shader.uniforms.selectionPivot=uniforms.pivot;
  shader.vertexShader='attribute float partIndex; uniform sampler2D vrState; uniform sampler2D rotationState; uniform sampler2D rotationPivots; uniform float stateWidth; uniform vec4 selectionRotation; uniform vec3 selectionPivot; varying float vrOpacity; vec3 rotateSelected(vec3 v){return v+2.0*cross(selectionRotation.xyz,cross(selectionRotation.xyz,v)+selectionRotation.w*v);}\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nif(texture2D(rotationState,vec2((partIndex+0.5)/stateWidth,0.5)).r>0.5)objectNormal=rotateSelected(objectNormal);');
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvec2 uvState=vec2((partIndex+0.5)/stateWidth,0.5);vec4 state=texture2D(vrState,uvState);transformed+=state.xyz;vrOpacity=state.w;if(texture2D(rotationState,uvState).r>0.5){vec3 pivot=texture2D(rotationPivots,uvState).xyz;transformed=pivot+rotateSelected(transformed-pivot);}');
  shader.fragmentShader='varying float vrOpacity;\n'+shader.fragmentShader;
  // Hashed opacity keeps context see-through in one depth-writing pass.
  shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\ndiffuseColor.a=vrOpacity;');
 };
 material.customProgramCacheKey=()=> 'quest-anatomy-v1';
 const mesh=new T.Mesh(geometry,material);mesh.frustumCulled=false;mesh.visible=false;mesh.matrixAutoUpdate=false;
 const sync=(state:SceneState,offsets:Float32Array,sourcePickers:(T.Mesh|undefined)[],context:number,solid=false)=>{
  const selected=new Set(state.selected);let changed=false;
  for(let i=0;i<model.parts.length;i++){
   const part=atlas.parts[i],sourceAlpha=(selected.has(part.id)?1:context)*partLayerOpacity(part,state)*(atlas.materials?.[part.material??'']?.opacity??1);
   // Quest's single-pass alpha hash looks like pixel noise in a headset.
   // Keep its visible parts solid; the ray peel then exposes the next layer.
   const alpha=solid?(sourceAlpha>=.5?1:0):sourceAlpha,visible=offsets[i*4+3]>.5&&alpha>=.001;
   if(shown[i]!==Number(visible)){shown[i]=Number(visible);changed=true;}
   states.set([offsets[i*4],offsets[i*4+1],offsets[i*4+2],visible?alpha:0],i*4);
   const source=sourcePickers[i];if(source){pickers[i].matrix.copy(source.matrixWorld);pickers[i].matrixWorld.copy(source.matrixWorld);}
  }
  texture.needsUpdate=true;
  if(changed){let sourceOffset=0,count=0;for(let i=0;i<model.parts.length;i++){const length=model.parts[i].indexCount;if(shown[i]){activeIndices.set(originalIndices.subarray(sourceOffset,sourceOffset+length),count);count+=length;}sourceOffset+=length;}
   geometry.setDrawRange(0,count);indices.clearUpdateRanges();indices.addUpdateRange(0,count);indices.needsUpdate=true;
  }
 };
 const dispose=()=>{geometry.dispose();material.dispose();texture.dispose();pickMaterial.dispose();for(const picker of pickers)picker.geometry.dispose();};
 return {mesh,pickers,sync,dispose,triangles:model.triangles,sourceTriangles:model.sourceTriangles};
}
export type VrAnatomy=Awaited<ReturnType<typeof loadVrAnatomy>>;
