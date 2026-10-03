/** Shared illustrative texture projected in native metre coordinates.
 * Three samples per visible fragment; no UV buffer, geometry or histology added.
 */
export const triplanarSurfaceTexture={
 vertex:'varying vec3 surfaceMapPosition; varying vec3 surfaceMapNormal;\n',
 position:'\nsurfaceMapPosition=position;surfaceMapNormal=normal;',
 fragment:'varying vec3 surfaceMapPosition; varying vec3 surfaceMapNormal;\n',
 map:`
 #ifdef USE_MAP
 vec3 surfaceWeights=abs(normalize(surfaceMapNormal));
 surfaceWeights/=max(dot(surfaceWeights,vec3(1.0)),0.00001);
 vec3 surfaceCoords=surfaceMapPosition*60.0;
 vec4 sampledDiffuseColor=texture2D(map,surfaceCoords.yz)*surfaceWeights.x
  +texture2D(map,surfaceCoords.xz)*surfaceWeights.y
  +texture2D(map,surfaceCoords.xy)*surfaceWeights.z;
 diffuseColor*=sampledDiffuseColor;
 #endif
 `,
};
