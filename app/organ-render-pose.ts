import {Vector3,Vector4} from 'three';
import type {Part} from './anatomy';
import {organOrientation} from './organ-frame';

export type OrganRenderPose=Pick<Part,'geometryOrigin'|'geometryInstance'>;
type PoseShader={uniforms:Record<string,{value:unknown}>;vertexShader:string};

/** A material may be shared by parts with different current placements. */
export function organRenderPoseKey(part?:OrganRenderPose){
 return part?.geometryOrigin||part?.geometryInstance?`:pose:${JSON.stringify([part.geometryOrigin??[0,0,0],part.geometryInstance??null])}`:'';
}

/** Use the same authoritative part pose as picking, sections and geometry audits.
 * Legacy material placement fields describe an earlier export and must not move anatomy.
 */
export function applyOrganRenderPose(shader:PoseShader,part?:OrganRenderPose){
 shader.uniforms.geometryOrigin={value:new Vector3().fromArray(part?.geometryOrigin??[0,0,0])};
 if(!part?.geometryInstance)return;
 organOrientation(part); // Validate the common CPU/GPU placement contract.
 shader.uniforms.geometryInstanceRotation={value:new Vector4(...part.geometryInstance.rotation)};
 shader.uniforms.geometryInstanceScale={value:part.geometryInstance.scale};
 shader.vertexShader='uniform vec4 geometryInstanceRotation;uniform float geometryInstanceScale;vec3 instanceVector(vec3 v){return v+2.0*cross(geometryInstanceRotation.xyz,cross(geometryInstanceRotation.xyz,v)+geometryInstanceRotation.w*v);}\n'+shader.vertexShader;
 shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nobjectNormal=instanceVector(objectNormal);');
 shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\ntransformed=instanceVector(transformed*geometryInstanceScale);');
}
