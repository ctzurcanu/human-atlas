export type EmbryoSurface='embryo-skin'|'placenta-maternal'|'placenta-fetal'|'membrane'|'umbilical-covering'|'vessel';
// Procedural surface finish, independent of source UVs. It changes shading
// only; it does not assert a modeled wall, villus or microscopic structure.
export function embryoSurfaceShader(surface:EmbryoSurface){
 const skin=surface==='embryo-skin',maternal=surface==='placenta-maternal',placenta=maternal||surface==='placenta-fetal';
 const frequency=skin?1450:placenta?360:surface==='vessel'?900:1100;
 const amplitude=skin?.000018:maternal?.00014:placenta?.000035:.000012;
 return {pars:`
 varying vec3 embryoSurfacePosition;
 float embryoHash(vec3 p){p=fract(p*.1031);p+=dot(p,p.yzx+33.33);return fract((p.x+p.y)*p.z);}
 float embryoNoise(vec3 p){vec3 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);
 return mix(mix(mix(embryoHash(i),embryoHash(i+vec3(1,0,0)),f.x),mix(embryoHash(i+vec3(0,1,0)),embryoHash(i+vec3(1,1,0)),f.x),f.y),mix(mix(embryoHash(i+vec3(0,0,1)),embryoHash(i+vec3(1,0,1)),f.x),mix(embryoHash(i+vec3(0,1,1)),embryoHash(i+vec3(1,1,1)),f.x),f.y),f.z);}
 float embryoFinish(vec3 p){return .68*embryoNoise(p)+.32*embryoNoise(p*2.13);}
 `,color:`
 float embryoGrain=embryoFinish(embryoSurfacePosition*${frequency.toFixed(1)});
 float embryoMottle=embryoNoise(embryoSurfacePosition*${(placenta?80:125).toFixed(1)});
 diffuseColor.rgb*=1.0+${placenta?'.16':'.045'}*(embryoMottle-.5)+${maternal?'.11':'.035'}*(embryoGrain-.5);
 `,normal:`
 float embryoHeight=embryoFinish(embryoSurfacePosition*${frequency.toFixed(1)})*${amplitude.toFixed(7)};
 vec3 embryoQ0=dFdx(-vViewPosition),embryoQ1=dFdy(-vViewPosition);
 vec3 embryoS=cross(embryoQ1,normal),embryoT=cross(normal,embryoQ0);float embryoDet=dot(embryoQ0,embryoS);
 normal=normalize(abs(embryoDet)*normal-sign(embryoDet)*(dFdx(embryoHeight)*embryoS+dFdy(embryoHeight)*embryoT));
 `};
}
