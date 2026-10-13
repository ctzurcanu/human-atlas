import {MeshStandardMaterial,ShaderChunk} from 'three';

/** Source meshes can contain a zero shading normal even on a valid triangle.
 * Use that triangle's actual normal rather than normalize(0) in the shader.
 * Stored anatomy and valid smooth normals remain unchanged.
 */
export function applySafeSurfaceNormals(material:MeshStandardMaterial){
 const base=material.onBeforeCompile,key=material.customProgramCacheKey();
 material.customProgramCacheKey=()=>key+':finite-surface-normal-v1';
 material.onBeforeCompile=(shader,renderer)=>{
  base(shader,renderer);
  shader.fragmentShader=shader.fragmentShader.replace('#include <normal_fragment_begin>',
   ShaderChunk.normal_fragment_begin.replace('vec3 normal = normalize( vNormal );',`
    vec3 sourceFaceNormal=cross(dFdx(vViewPosition),dFdy(vViewPosition));
    vec3 normal=dot(vNormal,vNormal)>1e-12?normalize(vNormal):
      (dot(sourceFaceNormal,sourceFaceNormal)>1e-30?normalize(sourceFaceNormal)*faceDirection:vec3(0.0,0.0,1.0));
   `));
 };
}
