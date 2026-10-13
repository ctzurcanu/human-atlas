import {Color,MeshBasicMaterial,MeshStandardMaterial} from 'three';
import type {Part} from './anatomy';
import {isNativeWholeLiver} from './hepatic-relations';

// Authored gross tissue colors. Utah WebPath LIVER002/LIVER083 show a smooth
// brown normal exterior and cut surface. These are illustrative appearance
// fields, not a microscopic TA98 representation or sampled reference pixels.
export const HEPATIC_EXTERIOR_COLOR=0x6b3d2b;
export const HEPATIC_CUT_COLOR=0x824533;
export function hepaticSurfaceFinish(part:Part){
 return isNativeWholeLiver(part)?{hepaticSurface:true as const}:undefined;
}
export function applyHepaticSurface(material:MeshStandardMaterial){
 const previous=material.onBeforeCompile,key=material.customProgramCacheKey();
 material.color.setHex(0xffffff);material.metalness=0;material.roughness=.68;
 material.onBeforeCompile=(shader,renderer)=>{
  previous.call(material,shader,renderer);
  shader.uniforms.hepaticExteriorColor={value:new Color(HEPATIC_EXTERIOR_COLOR)};
  shader.fragmentShader='uniform vec3 hepaticExteriorColor;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>',`#include <color_fragment>
   float hepaticSourceGrain=clamp(dot(diffuseColor.rgb,vec3(.2126,.7152,.0722)),0.0,1.0);
   diffuseColor.rgb=hepaticExteriorColor*(.86+.14*hepaticSourceGrain);`);
 };
 material.customProgramCacheKey=()=>key+':hepatic-gross-surface-v1';
}

// Stable three-dimensional, low-contrast cut grain. One eight-corner value
// field, zero image samples. Derivative filtering removes unresolved detail.
// No depicted hexagonal lobules, vessels or histological structures are
// inferred from this surface finish. True channels still come from geometry.
export const hepaticCutGrainGLSL=`
 float hepaticHash(vec3 p){return fract(sin(dot(mod(p,97.0),vec3(127.1,311.7,74.7)))*43758.5453);}
 float hepaticNoise(vec3 p){
  vec3 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);
  return mix(mix(mix(hepaticHash(i),hepaticHash(i+vec3(1,0,0)),f.x),
                 mix(hepaticHash(i+vec3(0,1,0)),hepaticHash(i+vec3(1,1,0)),f.x),f.y),
             mix(mix(hepaticHash(i+vec3(0,0,1)),hepaticHash(i+vec3(1,0,1)),f.x),
                 mix(hepaticHash(i+vec3(0,1,1)),hepaticHash(i+vec3(1,1,1)),f.x),f.y),f.z);
 }
 float hepaticGrain(vec3 p){
  vec3 q=p/.0005,footprint=fwidth(q);
  float detail=1.0-smoothstep(.35,1.2,max(footprint.x,max(footprint.y,footprint.z)));
  return .96+.08*mix(.5,hepaticNoise(q),detail);
 }`;
export function applyHepaticCut(material:MeshBasicMaterial|MeshStandardMaterial){
 const previous=material.onBeforeCompile,key=material.customProgramCacheKey();
 material.onBeforeCompile=(shader,renderer)=>{
  previous.call(material,shader,renderer);
  shader.uniforms.hepaticCutColor={value:new Color(HEPATIC_CUT_COLOR)};
  shader.vertexShader='varying vec3 hepaticCutPoint;\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>',
   '#include <begin_vertex>\nhepaticCutPoint=(modelMatrix*vec4(transformed,1.0)).xyz;');
  shader.fragmentShader='uniform vec3 hepaticCutColor; varying vec3 hepaticCutPoint;\n'+hepaticCutGrainGLSL+shader.fragmentShader;
  // Insert after existing cap color operations, preserving clipping, stencil,
  // transparency and the original silhouette. Only appearance is replaced.
  shader.fragmentShader=shader.fragmentShader.replace('#include <alphamap_fragment>',
   'diffuseColor.rgb=hepaticCutColor*hepaticGrain(hepaticCutPoint);\n#include <alphamap_fragment>');
 };
 material.customProgramCacheKey=()=>key+':hepatic-gross-cut-v1';
}
