import {MeshStandardMaterial,Texture} from 'three';

/** Generated UV height relief. No vertex displacement, added topology or UV
 * duplication. The luminal mask prevents rugae on the serosal exterior. */
export function applyGastricRelief(material:MeshStandardMaterial,height:Texture,scale=.0012){
 if(!(scale>0&&Number.isFinite(scale)))throw new Error('Invalid gastric relief scale');
 material.bumpMap=height;material.bumpScale=scale;
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;
 material.customProgramCacheKey=()=>key.call(material)+':gastric-uv-relief-v3';
 material.onBeforeCompile=(shader,renderer)=>{
  compile(shader,renderer);
  shader.vertexShader='attribute float gastricLumen;varying float vGastricLumen;\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\nvGastricLumen=gastricLumen;');
  shader.fragmentShader='varying float vGastricLumen;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>',`#include <color_fragment>
   #ifdef USE_BUMPMAP
   vec2 gastricReliefUv=vBumpMapUv;
   if(vGastricLumen>.5){
    vec3 gastricDx=dFdx(-vViewPosition),gastricDy=dFdy(-vViewPosition);
    vec2 gastricUx=dFdx(vBumpMapUv),gastricUy=dFdy(vBumpMapUv);
    float gastricDet=gastricUx.x*gastricUy.y-gastricUx.y*gastricUy.x;
    vec2 gastricUv=vBumpMapUv;
    if(abs(gastricDet)>1e-10){
     vec3 gastricU=(gastricDx*gastricUy.y-gastricDy*gastricUx.y)/gastricDet;
     vec3 gastricV=(-gastricDx*gastricUy.x+gastricDy*gastricUx.x)/gastricDet;
     vec3 gastricN=normalize(cross(gastricDx,gastricDy))*(gl_FrontFacing?1.0:-1.0);
     vec3 gastricEye=normalize(vViewPosition);
     float gastricFacing=max(abs(dot(gastricEye,gastricN)),.25);
     vec2 gastricRay=vec2(dot(gastricEye,normalize(gastricU))/max(length(gastricU),1e-6),dot(gastricEye,normalize(gastricV))/max(length(gastricV),1e-6))*bumpScale/gastricFacing;
     gastricRay=clamp(gastricRay,vec2(-.035),vec2(.035));
     float gastricDepth=0.0;
     for(int step=0;step<8;step++){
      float gastricSurface=1.0-texture2D(bumpMap,gastricUv).r;
      if(gastricDepth>=gastricSurface)break;
      gastricUv-=gastricRay*.125;gastricDepth+=.125;
     }
     if(gastricDepth>0.0){
      float gastricAfter=gastricDepth-1.0+texture2D(bumpMap,gastricUv).r;
      vec2 gastricPreviousUv=gastricUv+gastricRay*.125;
      float gastricBefore=gastricDepth-.125-1.0+texture2D(bumpMap,gastricPreviousUv).r;
      float gastricSpan=gastricAfter-gastricBefore;
      if(abs(gastricSpan)>1e-5)gastricUv=mix(gastricUv,gastricPreviousUv,clamp(gastricAfter/gastricSpan,0.0,1.0));
     }
    }
    gastricReliefUv=gastricUv;
    float gastricHeight=texture2D(bumpMap,gastricUv).r;
    diffuseColor.rgb=mix(vec3(.46,.12,.085),vec3(.70,.32,.20),gastricHeight);
   }
   #endif`);
  shader.fragmentShader=shader.fragmentShader.replace('#include <normal_fragment_maps>',`vec3 gastricUnperturbedNormal=normal;
   #include <normal_fragment_maps>
   #ifdef USE_BUMPMAP
   if(vGastricLumen>.5){
    float gastricBase=bumpScale*texture2D(bumpMap,gastricReliefUv).r;
    vec2 gastricGradient=vec2(
     bumpScale*texture2D(bumpMap,gastricReliefUv+dFdx(vBumpMapUv)).r-gastricBase,
     bumpScale*texture2D(bumpMap,gastricReliefUv+dFdy(vBumpMapUv)).r-gastricBase);
    normal=perturbNormalArb(-vViewPosition,gastricUnperturbedNormal,gastricGradient,faceDirection);
   }else normal=gastricUnperturbedNormal;
   #endif`);
 };
}
