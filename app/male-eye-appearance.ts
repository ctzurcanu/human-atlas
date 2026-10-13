import {MeshStandardMaterial} from 'three';

/** Cosmetic optical darkness on the existing ZA eye surface; adds no pupil solid. */
export function applyMaleEyeAppearance(material:MeshStandardMaterial){
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>`${key.call(material)}:native-male-pupil-v1`;
 material.onBeforeCompile=(shader,renderer)=>{
  compile.call(material,shader,renderer);
  shader.vertexShader='varying vec3 nativeEyePosition;\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nnativeEyePosition=position;');
  shader.fragmentShader='varying vec3 nativeEyePosition;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <opaque_fragment>',`
   float nativePupilRadius=length(vec2(abs(nativeEyePosition.x)-.033358805,nativeEyePosition.y-1.703712772));
   float nativePupil=(1.0-smoothstep(.0021,.0023,nativePupilRadius))*step(.073,nativeEyePosition.z);
   outgoingLight=mix(outgoingLight,vec3(.001),nativePupil);
   #include <opaque_fragment>
  `);
 };
}
