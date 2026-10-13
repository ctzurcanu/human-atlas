import {Color,MeshStandardMaterial} from 'three';

/** Native ZA oral-skin pigment only; follows the surface through all presentation transforms. */
export function applyMaleLipAppearance(material:MeshStandardMaterial){
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>`${key.call(material)}:native-male-vermilion-v1`;
 material.onBeforeCompile=(shader,renderer)=>{
  compile.call(material,shader,renderer);
  shader.uniforms.nativeLipColor={value:new Color('#a15b54')};
  shader.vertexShader='varying vec3 nativeLipPosition;\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nnativeLipPosition=position;');
  shader.fragmentShader='uniform vec3 nativeLipColor;varying vec3 nativeLipPosition;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',`#include <map_fragment>
   float lipX=abs(nativeLipPosition.x)/.027;
   float lipArch=sqrt(max(0.0,1.0-lipX*lipX));
   float lipTop=1.638+.006*lipArch-.0015*exp(-pow(nativeLipPosition.x/.004,2.0));
   float lipBottom=1.638-.010*lipArch;
   float lipPigment=(1.0-smoothstep(.97,1.02,lipX))*smoothstep(lipBottom-.0005,lipBottom+.0005,nativeLipPosition.y)*(1.0-smoothstep(lipTop-.0005,lipTop+.0005,nativeLipPosition.y))*smoothstep(.085,.092,nativeLipPosition.z);
   diffuseColor.rgb=mix(diffuseColor.rgb,nativeLipColor,lipPigment*.85);
  `);
 };
}
