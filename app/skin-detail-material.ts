import {MeshStandardMaterial,Texture,Vector3} from 'three';

export interface SkinDetail {map:string;repeat:number;strength:number;mean:[number,number,number]}

/** A small repeatable albedo detail map preserves the atlas's anatomical marks.
 * No additional vertex attributes, geometry or CPU animation are required. */
export function applySkinDetail(material:MeshStandardMaterial,texture:Texture,detail:SkinDetail){
 if(!material.map)return;
 if(!(detail.repeat>0)||!Number.isFinite(detail.repeat)||!Number.isFinite(detail.strength)||detail.strength<0||detail.strength>1||detail.mean.some(c=>!Number.isFinite(c)||c<=0))throw Error('Invalid skin detail');
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>`${key.call(material)}:skin-detail:${JSON.stringify(detail)}`;
 material.onBeforeCompile=(shader,renderer)=>{
  compile.call(material,shader,renderer);
  shader.uniforms.skinDetailMap={value:texture};shader.uniforms.skinDetailMean={value:new Vector3(...detail.mean)};
  shader.uniforms.skinDetailRepeat={value:detail.repeat};shader.uniforms.skinDetailStrength={value:detail.strength};
  shader.fragmentShader='uniform sampler2D skinDetailMap; uniform vec3 skinDetailMean; uniform float skinDetailRepeat; uniform float skinDetailStrength;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',`#include <map_fragment>
   #ifdef USE_MAP
    vec3 skinMicroColor=texture2D(skinDetailMap,vMapUv*skinDetailRepeat).rgb;
    vec3 skinMicroRatio=clamp(skinMicroColor/skinDetailMean,vec3(.88),vec3(1.12));
    diffuseColor.rgb*=mix(vec3(1.0),skinMicroRatio,skinDetailStrength);
   #endif`);
 };
}
