import {MeshStandardMaterial,Vector2,Vector3} from 'three';
export interface SkinEyeAppearance {centers:[number,number][];halfWidth:number;halfHeight:number;contour?:'eyelid';irisRadius?:number;pupilRadius?:number}
export function skinEyeHoleAt(eyes:SkinEyeAppearance,uv:{x:number;y:number}):boolean{
 return eyes.centers.some(([x,y],i)=>{const u=(uv.x-x)/eyes.halfWidth,vertical=uv.y-y-(eyes.contour==='eyelid'?.00035*u*(i===0?-1:1):0);const arch=Math.pow(Math.max(0,1-u*u),eyes.contour==='eyelid'?.75:.5);return Math.abs(u)<1&&Math.abs(vertical)<eyes.halfHeight*arch*(eyes.contour==='eyelid'&&vertical<0?.78:1);});
}
/** Binary eyelid openings: discard exposes real underlying anatomy without sorting. */
export function applySkinEyeAppearance(material:MeshStandardMaterial,eyes:SkinEyeAppearance){
 if(!material.map)return;
 if(!eyes.centers.length||eyes.centers.length>2||[eyes.halfWidth,eyes.halfHeight].some(v=>!Number.isFinite(v)||v<=0))throw Error('Invalid skin eye openings');
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>`${key.call(material)}:skin-eye-openings:${JSON.stringify(eyes)}`;
 material.onBeforeCompile=(shader,renderer)=>{
  compile.call(material,shader,renderer);
  shader.uniforms.skinEyeCenters={value:eyes.centers.map(p=>new Vector2(...p))};shader.uniforms.skinEyeSize={value:new Vector2(eyes.halfWidth,eyes.halfHeight)};
  shader.fragmentShader=`uniform vec2 skinEyeCenters[${eyes.centers.length}]; uniform vec2 skinEyeSize;\n`+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',`#include <map_fragment>
   #ifdef USE_MAP
    for(int eye=0;eye<${eyes.centers.length};eye++){
     vec2 delta=vMapUv-skinEyeCenters[eye];float x=delta.x/skinEyeSize.x;
     float eyeVertical=delta.y${eyes.contour==='eyelid'?'- .00035*x*(eye==0?-1.0:1.0)':''};
     float aperture=skinEyeSize.y*pow(max(0.0,1.0-x*x),${eyes.contour==='eyelid'?'0.75':'0.5'})${eyes.contour==='eyelid'?'*(eyeVertical<0.0?.78:1.0)':''};
     if(abs(x)<1.0&&abs(eyeVertical)<aperture)discard;
     // Remove the legacy painted sclera outside the narrower physical opening.
     float oldX=delta.x/.027;float oldHeight=.015*sqrt(max(0.0,1.0-oldX*oldX));
     if(abs(oldX)<1.0&&abs(delta.y)<oldHeight){
      vec3 original=texture2D(map,vMapUv).rgb;vec3 gain=diffuseColor.rgb/max(original,vec3(.01));
      diffuseColor.rgb=texture2D(map,vec2(vMapUv.x,skinEyeCenters[eye].y-.030)).rgb*gain;
     }
    }
   #endif`);
 };
}
export function applySkinToneGain(material:MeshStandardMaterial,gain:[number,number,number],neckGain?:[number,number,number]){
 if(gain.some(v=>!Number.isFinite(v)||v<=0))throw Error('Invalid skin tone gain');
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>`${key.call(material)}:skin-tone:${gain.join(',')}:${neckGain?.join(',')??''}`;
 material.onBeforeCompile=(shader,renderer)=>{compile.call(material,shader,renderer);shader.uniforms.skinToneGain={value:new Vector3(...gain)};
  if(neckGain){shader.uniforms.skinToneNeckGain={value:new Vector3(...neckGain)};shader.vertexShader='varying float skinToneHeight;\n'+shader.vertexShader;shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\n skinToneHeight=position.y;');shader.fragmentShader='uniform vec3 skinToneNeckGain; varying float skinToneHeight;\n'+shader.fragmentShader;}
  shader.fragmentShader='uniform vec3 skinToneGain;\n'+shader.fragmentShader;shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>','#include <map_fragment>\n diffuseColor.rgb*='+ (neckGain?'mix(skinToneNeckGain,skinToneGain,smoothstep(1.55,1.60,skinToneHeight))':'skinToneGain')+';');};
}
