import {MeshStandardMaterial} from 'three';
/** Low-resolution whole-globe appearance; no separate histological coat claim. */
export function applyOcularSurface(material:MeshStandardMaterial){
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>`${key.call(material)}:ocular-surface-v1`;
 material.onBeforeCompile=(shader,renderer)=>{
  compile.call(material,shader,renderer);
  shader.vertexShader='varying vec3 ocularPoint;\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\n ocularPoint=position;');
  shader.fragmentShader='varying vec3 ocularPoint;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',`#include <map_fragment>
   if(ocularPoint.z>.006){
    float radius=length(ocularPoint.xy);float iris=1.0-smoothstep(.0055,.0058,radius);float pupil=1.0-smoothstep(.0021,.0023,radius);
    float fiber=.8+.2*sin(atan(ocularPoint.y,ocularPoint.x)*52.0+radius*3100.0);
    vec3 color=mix(vec3(.82,.80,.74),vec3(.065,.026,.009)*fiber,iris);color=mix(color,vec3(.002),pupil);diffuseColor.rgb=color;
   }`);
 };
}
