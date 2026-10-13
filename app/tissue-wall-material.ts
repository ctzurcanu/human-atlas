import {MeshStandardMaterial, Matrix4, SRGBColorSpace} from 'three';

export type TissuePattern = 'solid'|'cellular'|'longitudinal'|'circular'|'oblique'|'fibrous'|'mucosal'|'intestinal'|'gallbladder-folds';
export interface TissueAppearance {color: [number,number,number]; pattern: TissuePattern; scale: number; relief: number; courseLength?:number;foldFrame?:{origin:[number,number,number];axes:[[number,number,number],[number,number,number],[number,number,number]]}}

/** Two procedural fields (color and relief), evaluated in unchanged organ
 * coordinates. A cap can use the same material with its tissue frame supplied.
 * No UV buffer, sampled image, or copied reference pixels are needed.
 */
export function tissueWallMaterial(appearance: TissueAppearance, frame?: Matrix4, inner=false) {
  if (!(appearance.scale>0) || !Number.isFinite(appearance.scale) || appearance.relief<0
    || !Number.isFinite(appearance.relief) || appearance.color.some(c=>c<0 || c>1 || !Number.isFinite(c)))
    throw new Error('Invalid tissue appearance');
  if(appearance.courseLength!==undefined&&(!(appearance.courseLength>0)||!Number.isFinite(appearance.courseLength)))throw new Error('Invalid tissue course length');
  if(appearance.pattern==='gallbladder-folds'&&(!appearance.foldFrame||appearance.foldFrame.origin.length!==3||appearance.foldFrame.axes.length!==3||appearance.foldFrame.axes.some(a=>a.length!==3)||[...appearance.foldFrame.origin,...appearance.foldFrame.axes.flat()].some(x=>!Number.isFinite(x))))throw new Error('Invalid gallbladder fold frame');
  const material = new MeshStandardMaterial({roughness:.82});
  material.color.setRGB(...appearance.color,SRGBColorSpace);
  const tissueFrame = {value: frame ?? new Matrix4()};
  const tissueWorld = {value: frame ? 1 : 0};
  material.onBeforeCompile = shader => {
    shader.uniforms.tissueFrame = tissueFrame; shader.uniforms.tissueWorld = tissueWorld;
    shader.vertexShader = 'uniform mat4 tissueFrame; uniform float tissueWorld; varying vec3 tissuePoint;\n'+shader.vertexShader;
    shader.vertexShader = shader.vertexShader.replace('#include <begin_vertex>',
      '#include <begin_vertex>\ntissuePoint=tissueWorld>0.5?(tissueFrame*modelMatrix*vec4(position,1.0)).xyz:position;');
    if(appearance.courseLength!==undefined){
      shader.vertexShader='attribute vec3 tissueChart; varying vec3 tissueCourse;\n'+shader.vertexShader;
      shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\ntissueCourse=tissueChart;');
    }
    if(inner)shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>',
      '#include <beginnormal_vertex>\nobjectNormal=-objectNormal;');
    const field = appearance.pattern==='gallbladder-folds'?'gallbladderFoldField(p)':appearance.pattern==='solid'?'0.5':appearance.pattern==='longitudinal'
      ? '.5+.5*sin(tissueAngle(p.xz)*30.0+sin(p.y)*.15)'
      : appearance.pattern==='circular'?'.5+.5*sin(p.y*6.2831853+.12*sin(p.x))'
      : appearance.pattern==='oblique'?'.5+.5*sin(p.y*4.5+tissueAngle(p.xz)*24.0)'
      : appearance.pattern==='fibrous'?'.5+.30*sin(p.y*4.0+p.x*1.3+.2*sin(p.z*3.0))+.20*sin(p.y*7.7-p.z*1.5)'
      : appearance.pattern==='mucosal'?'.5+.5*sin(p.x*6.2831853)*sin(p.y*6.2831853)*sin(p.z*6.2831853)'
      : '.5+.25*sin(p.x*3.7+p.y*2.9)+.25*sin(p.z*4.1-p.y*3.3)';
    shader.fragmentShader = `varying vec3 tissuePoint;
      float tissueAngle(vec2 p){return dot(p,p)>1e-12?atan(p.y,p.x):0.0;}
      float tissueField(vec3 p){return ${field};}\n`+shader.fragmentShader;
    if(appearance.pattern==='gallbladder-folds'){
      const frame=appearance.foldFrame!,vector=(v:number[])=>`vec3(${v.map(n=>n.toPrecision(10)).join(',')})`,origin=vector(frame.origin.map(n=>n/appearance.scale));
      shader.fragmentShader=`
       float gallbladderFoldField(vec3 point){
        vec3 q=point-${origin};
        vec3 p=vec3(dot(q,${vector(frame.axes[0])}),dot(q,${vector(frame.axes[1])}),dot(q,${vector(frame.axes[2])}));
        float phase=6.2831853*p.x/3.4+1.1*sin(6.2831853*p.y/9.0)+.6*sin(6.2831853*p.z/8.0);
        float branch=6.2831853*p.y/5.2+.9*sin(6.2831853*p.x/8.0)+.5*sin(6.2831853*p.z/6.0);
        return max(pow(.5+.5*cos(phase),5.0),.9*pow(.5+.5*cos(branch),6.0)*(.5+.5*sin(6.2831853*p.x/13.0)));
       }\n`+shader.fragmentShader;
    }
    shader.fragmentShader = shader.fragmentShader.replace('#include <color_fragment>',
      `#include <color_fragment>\nfloat tissueDetail=tissueField(tissuePoint/${appearance.scale.toPrecision(8)});diffuseColor.rgb*=${appearance.pattern==='gallbladder-folds'?'.72+.38*tissueDetail':'.92+.16*tissueDetail'};`);
    if(appearance.courseLength!==undefined){
      // Unit-circle coordinates interpolate across the circumferential seam.
      // World-frame section caps have no interpolated course and keep their
      // existing local cut field; only carrier surfaces use this chart.
      const courseField=appearance.pattern==='longitudinal'?'.5+.5*sin(angle*30.0+.12*sin(course*.2))'
        :appearance.pattern==='circular'?'.5+.5*sin(course*6.2831853+.12*sin(angle*3.0))'
        :appearance.pattern==='fibrous'?'.5+.30*sin(course*4.0+angle*13.0+.2*sin(angle*3.0))+.20*sin(course*7.7-angle*5.0)'
        :appearance.pattern==='intestinal'?'intestinalSurfaceField(course,angle)'
        :appearance.pattern==='mucosal'?'.5+.25*sin(course*6.2831853+sin(angle*7.0))+.25*sin(angle*35.0+.3*sin(course*4.1))'
        :appearance.pattern==='cellular'?'.5+.25*sin(course*3.7+angle*29.0)+.25*sin(angle*41.0-course*3.3)'
        :appearance.pattern==='oblique'?'.5+.5*sin(course*4.5+angle*24.0)':'0.5';
      shader.fragmentShader='uniform float tissueWorld; varying vec3 tissueCourse;\n'+shader.fragmentShader;
      if(appearance.pattern==='intestinal')shader.fragmentShader=`
        // Continuous at the angular seam: integer circumference cycles and
        // periodic warping. This is shallow surface lighting, not villus mesh.
        float intestinalSurfaceField(float course,float angle){
         vec2 coordinates=vec2(course+.17*sin(angle*17.0),angle*(192.0/6.28318530718)+.17*sin(course*.9));
         vec2 id=floor(coordinates);id.y=mod(id.y,192.0);
         vec2 jitter=fract(sin(vec2(dot(id,vec2(127.1,311.7)),dot(id,vec2(269.5,183.3))))*43758.5453);
         vec2 cell=fract(coordinates)-.5;
         float window=max(0.0,1.0-4.0*cell.x*cell.x)*max(0.0,1.0-4.0*cell.y*cell.y);
         cell+=(jitter-.5)*.25*window;
         float rounded=max(0.0,1.0-4.0*dot(cell,cell));
         return rounded*rounded*(.75+.25*jitter.x);
        }\n`+shader.fragmentShader;
      shader.fragmentShader=shader.fragmentShader.replace('diffuseColor.rgb*=.92+.16*tissueDetail;',
        `if(tissueWorld<.5){float course=tissueCourse.x*${(appearance.courseLength/appearance.scale).toPrecision(8)};float angle=dot(tissueCourse.yz,tissueCourse.yz)>1e-12?atan(tissueCourse.z,tissueCourse.y):0.0;tissueDetail=${courseField};${appearance.pattern==='intestinal'?'vec2 circle=tissueCourse.yz;float radius2=max(dot(circle,circle),1e-8);float angularDx=(circle.x*dFdx(circle.y)-circle.y*dFdx(circle.x))/radius2;float angularDy=(circle.x*dFdy(circle.y)-circle.y*dFdy(circle.x))/radius2;float footprint=max(fwidth(course),max(abs(angularDx),abs(angularDy))*(192.0/6.28318530718));tissueDetail=mix(.26,tissueDetail,1.0-smoothstep(.4,1.5,footprint));':''}}diffuseColor.rgb*=.92+.16*tissueDetail;`);
    }
    shader.fragmentShader = shader.fragmentShader.replace('#include <normal_fragment_maps>',
      `#include <normal_fragment_maps>
       float height=tissueDetail*${appearance.relief.toPrecision(8)};
       vec3 dx=dFdx(-vViewPosition),dy=dFdy(-vViewPosition);
       vec3 rx=cross(dy,normal),ry=cross(normal,dx);
       float det=dot(dx,rx);
       if(abs(det)>1e-18)normal=normalize(abs(det)*normal-sign(det)*(dFdx(height)*rx+dFdy(height)*ry));`);
  };
  material.customProgramCacheKey = () => 'tissue-wall-v1:'+JSON.stringify(appearance)+':'+Number(inner);
  return {material, tissueFrame};
}
