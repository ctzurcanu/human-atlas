import {MeshStandardMaterial, Matrix4, Vector2, Vector3, SRGBColorSpace} from 'three';

export type TissuePattern = 'solid'|'cellular'|'longitudinal'|'circular'|'oblique'|'fibrous'|'mucosal'|'intestinal'|'gallbladder-folds'|'thyroid-follicles';
export interface ThyroidFollicleAppearance {
  /** Physical metres in the unchanged tissue coordinate frame. */
  diameterMin?:number; diameterMax?:number;
  /** Visual shell fraction of the radius, not a measured epithelial thickness. */
  epithelialFraction?:number;
}
export interface TissueAppearance {color: [number,number,number]; pattern: TissuePattern; scale: number; relief: number; courseLength?:number;foldFrame?:{origin:[number,number,number];axes:[[number,number,number],[number,number,number],[number,number,number]]};thyroidFollicles?:ThyroidFollicleAppearance;fibres?:{directions:[[number,number,number],[number,number,number]];weights?:[number,number];control?:'tongue-tv'}}

function thyroidFollicleParameters(appearance:TissueAppearance){
  const {diameterMin=.00005,diameterMax=.0005,epithelialFraction=.06}=appearance.thyroidFollicles??{};
  if(![diameterMin,diameterMax,epithelialFraction].every(Number.isFinite)
    ||diameterMin<.00005||diameterMax>.0005||diameterMin>diameterMax
    ||epithelialFraction<.015||epithelialFraction>.2||appearance.courseLength!==undefined)
    throw new Error('Invalid thyroid follicle appearance');
  // Spherical follicles 50–500 µm, colloid and cuboidal epithelium:
  // https://histologyguide.com/slideview/MH-151-thyroid/13-slide-1.html
  // This field supplies visual tissue organization, not individual cell meshes.
  // Spheres stay inside disjoint FCC containers. Smaller follicles can move
  // farther from their sites; additional small ones fill the octahedral voids.
  // The placement is visual, not a reconstruction of an individual's histology.
  const power=.3,padding=Math.min(diameterMax*.004,(diameterMax-diameterMin)*.05),stepLength=diameterMax/Math.SQRT2;
  const meanCube=(minimum:number,maximum:number)=>{
    const r0=minimum/2,dr=(maximum-minimum)/2;
    return r0**3+3*r0*r0*dr/(1+power)+3*r0*dr*dr/(1+2*power)+dr**3/(1+3*power);
  };
  const smallMaximum=2*(stepLength-diameterMax/2-2*padding);
  const meanRadiusCube=meanCube(diameterMin,diameterMax-2*padding)+(smallMaximum>=diameterMin?meanCube(diameterMin,smallMaximum):0);
  const coverage=(4*Math.PI/3)*meanRadiusCube/(2*stepLength**3);
  const coreAverage=coverage*(1-epithelialFraction)**3,ringAverage=coverage-coreAverage;
  return {diameterMin,diameterMax,epithelialFraction,coreAverage,ringAverage,padding};
}

function thyroidFollicleShader(appearance:TissueAppearance){
  const p=thyroidFollicleParameters(appearance),n=(v:number)=>v.toPrecision(10);
  return `
    vec3 thyroidFollicleHash(vec3 p){
      // Polynomial hashing avoids trigonometric precision noise at large origins.
      vec3 q=fract(p*vec3(.1031,.1030,.0973));q+=dot(q,q.yxz+33.33);
      return fract((q.xxy+q.yzz)*q.zyx);
    }
    vec3 thyroidFollicleField(vec3 point){
      const float minimum=${n(p.diameterMin/appearance.scale)};
      const float maximum=${n(p.diameterMax/appearance.scale)};
      const float shell=${n(p.epithelialFraction)};
      const float padding=${n(p.padding/appearance.scale)};
      const vec3 average=vec3(${n(p.coreAverage)},${n(p.ringAverage)},0.0);
      float footprint=max(length(dFdx(point)),length(dFdy(point)));
      float detail=1.0-smoothstep(maximum*.25,maximum*.9,footprint);
      // No sphere neighbourhood is evaluated once follicles are subpixel.
      if(detail<=.001)return average;
      float stepLength=maximum*.7071067811865475;
      // An orthonormal rotation avoids aligned rows on common section planes.
      // It preserves physical diameters and does not move the tissue frame.
      mat3 packingFrame=mat3(
        vec3(.6424709519829305,.742614819633407,.1890880893102192),
        vec3(-.6537507461112396,.40241799334757167,.6408351742758127),
        vec3(.39980124789343074,-.5353344639484331,.7440268636896784));
      vec3 q=(packingFrame*point)/stepLength,base=floor(q);
      float nearest=1e6,nearestRadius=maximum*.5;
      for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++)for(int z=-1;z<=1;z++){
        vec3 id=base+vec3(float(x),float(y),float(z));
        bool large=mod(id.x+id.y+id.z,2.0)<.5;
        float container=large?maximum*.5:stepLength-maximum*.5-padding;
        float maximumDiameter=2.0*(container-padding);
        if(maximumDiameter>=minimum){
          vec3 random=thyroidFollicleHash(id);
          vec3 sizeRandom=thyroidFollicleHash(id+vec3(19.13,3.73,11.17));
          float diameter=mix(minimum,maximumDiameter,pow(sizeRandom.x,.3));
          float radius=diameter*.5;
          vec3 direction=random*2.0-1.0;
          direction/=max(length(direction),.0001);
          float travel=max(0.0,container-radius-padding)*pow(sizeRandom.y,.3333333333333333);
          vec3 centre=id+direction*(travel/stepLength);
          float distance=length((q-centre)*stepLength)/radius;
          if(distance<nearest){nearest=distance;nearestRadius=radius;}
        }
      }
      // Small follicles must disappear earlier than the largest ones. Limit the
      // edge filter rather than expanding a tiny lining across its whole core.
      float localDetail=detail*(1.0-smoothstep(nearestRadius*.5,nearestRadius*1.8,footprint));
      float edge=min(.25,max(.002,footprint/nearestRadius*.65));
      float tissue=1.0-smoothstep(1.0-edge,1.0+edge,nearest);
      float core=1.0-smoothstep(1.0-shell-edge,1.0-shell+edge,nearest);
      core=min(core,tissue);
      float ring=max(0.0,tissue-core);
      // Relief is visual surface lighting; the true section geometry is unchanged.
      float relief=(.22*core+.5*ring)*localDetail;
      return vec3(mix(average.xy,vec2(core,ring),localDetail),relief);
    }
    vec3 thyroidFollicleColor(vec3 color,vec3 field){
      vec3 stroma=color*vec3(.88,.80,.74);
      vec3 epithelium=color*vec3(.72,.60,.52);
      vec3 colloid=mix(color*.95,vec3(.62,.54,.42),.3);
      return stroma*(1.0-field.x-field.y)+colloid*field.x+epithelium*field.y;
    }\n`;
}

/** Two procedural fields (color and relief), evaluated in unchanged organ
 * coordinates. A cap can use the same material with its tissue frame supplied.
 * No UV buffer, sampled image, or copied reference pixels are needed.
 */
export function tissueWallMaterial(appearance: TissueAppearance, frame?: Matrix4, inner=false, localOrigin:readonly number[]=[0,0,0]) {
  if(localOrigin.length!==3||localOrigin.some(value=>!Number.isFinite(value)))throw new Error('Invalid tissue coordinate origin');
  if (!(appearance.scale>0) || !Number.isFinite(appearance.scale) || appearance.relief<0
    || !Number.isFinite(appearance.relief) || appearance.color.some(c=>c<0 || c>1 || !Number.isFinite(c)))
    throw new Error('Invalid tissue appearance');
  if(appearance.courseLength!==undefined&&(!(appearance.courseLength>0)||!Number.isFinite(appearance.courseLength)))throw new Error('Invalid tissue course length');
  if(appearance.pattern==='gallbladder-folds'&&(!appearance.foldFrame||appearance.foldFrame.origin.length!==3||appearance.foldFrame.axes.length!==3||appearance.foldFrame.axes.some(a=>a.length!==3)||[...appearance.foldFrame.origin,...appearance.foldFrame.axes.flat()].some(x=>!Number.isFinite(x))))throw new Error('Invalid gallbladder fold frame');
  if(appearance.pattern==='thyroid-follicles')thyroidFollicleParameters(appearance);
  else if(appearance.thyroidFollicles!==undefined)throw new Error('Thyroid follicle parameters require thyroid-follicles pattern');
  const fibres=appearance.fibres;
  if(fibres&&(appearance.pattern!=='fibrous'||appearance.courseLength!==undefined||fibres.directions.length!==2||fibres.directions.some(d=>d.length!==3||d.some(n=>!Number.isFinite(n))||Math.hypot(...d)<1e-10)||fibres.weights&&(fibres.weights.length!==2||fibres.weights.some(n=>!Number.isFinite(n)||n<0)||Math.abs(fibres.weights[0]+fibres.weights[1]-1)>1e-10)))throw new Error('Invalid directional fibre field');
  const fibreWeights={value:new Vector2(...(fibres?.weights??[.5,.5]))};
  const material = new MeshStandardMaterial({roughness:.82});
  material.color.setRGB(...appearance.color,SRGBColorSpace);
  const tissueFrame = {value: frame ?? new Matrix4()};
  const tissueWorld = {value: frame ? 1 : 0};
  const tissueLocalOrigin = {value: new Vector3().fromArray(localOrigin)};
  material.onBeforeCompile = shader => {
    shader.uniforms.tissueFrame = tissueFrame; shader.uniforms.tissueWorld = tissueWorld; shader.uniforms.tissueLocalOrigin=tissueLocalOrigin;
    if(fibres)shader.uniforms.tissueFibreWeights=fibreWeights;
    shader.vertexShader = 'uniform mat4 tissueFrame; uniform float tissueWorld; uniform vec3 tissueLocalOrigin; varying vec3 tissuePoint;\n'+shader.vertexShader;
    shader.vertexShader = shader.vertexShader.replace('#include <begin_vertex>',
      '#include <begin_vertex>\ntissuePoint=(tissueWorld>0.5?(tissueFrame*modelMatrix*vec4(position,1.0)).xyz:position)+tissueLocalOrigin;');
    if(appearance.courseLength!==undefined){
      shader.vertexShader='attribute vec3 tissueChart; varying vec3 tissueCourse;\n'+shader.vertexShader;
      shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\ntissueCourse=tissueChart;');
    }
    if(inner)shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>',
      '#include <beginnormal_vertex>\nobjectNormal=-objectNormal;');
    const vector=(d:number[])=>`vec3(${d.map(n=>n.toPrecision(10)).join(',')})`;
    const field = fibres?`tissueFibreWeights.x*tissueDirectionalFibre(p,normalize(${vector(fibres.directions[0])}))+tissueFibreWeights.y*tissueDirectionalFibre(p,normalize(${vector(fibres.directions[1])}))`:appearance.pattern==='thyroid-follicles'?'thyroidFollicleField(p).z':appearance.pattern==='gallbladder-folds'?'gallbladderFoldField(p)':appearance.pattern==='solid'?'0.5':appearance.pattern==='longitudinal'
      ? '.5+.5*sin(tissueAngle(p.xz)*30.0+sin(p.y)*.15)'
      : appearance.pattern==='circular'?'.5+.5*sin(p.y*6.2831853+.12*sin(p.x))'
      : appearance.pattern==='oblique'?'.5+.5*sin(p.y*4.5+tissueAngle(p.xz)*24.0)'
      : appearance.pattern==='fibrous'?'.5+.30*sin(p.y*4.0+p.x*1.3+.2*sin(p.z*3.0))+.20*sin(p.y*7.7-p.z*1.5)'
      : appearance.pattern==='mucosal'?'.5+.5*sin(p.x*6.2831853)*sin(p.y*6.2831853)*sin(p.z*6.2831853)'
      : '.5+.25*sin(p.x*3.7+p.y*2.9)+.25*sin(p.z*4.1-p.y*3.3)';
    shader.fragmentShader = `varying vec3 tissuePoint;
      float tissueAngle(vec2 p){return dot(p,p)>1e-12?atan(p.y,p.x):0.0;}
      float tissueField(vec3 p){return ${field};}\n`+shader.fragmentShader;
    if(fibres)shader.fragmentShader=`
      uniform vec2 tissueFibreWeights;
      float tissueDirectionalFibre(vec3 p,vec3 direction){
        vec3 u=normalize(cross(direction,abs(direction.y)<.8?vec3(0.0,1.0,0.0):vec3(0.0,0.0,1.0)));
        vec3 v=cross(direction,u);
        vec2 across=vec2(dot(p,u),dot(p,v));
        float along=dot(p,direction);
        float detail=.5+.25*sin(across.x*6.2831853+.12*sin(along*.17))+.25*sin(across.y*6.2831853+.15*sin(along*.13));
        float footprint=max(fwidth(across.x),fwidth(across.y));
        return mix(detail,.5,smoothstep(.25,.8,footprint));
      }\n`+shader.fragmentShader;
    if(appearance.pattern==='thyroid-follicles')shader.fragmentShader=thyroidFollicleShader(appearance)+shader.fragmentShader;
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
      appearance.pattern==='thyroid-follicles'
        ?`#include <color_fragment>\nvec3 thyroidDetail=thyroidFollicleField(tissuePoint/${appearance.scale.toPrecision(8)});float tissueDetail=thyroidDetail.z;diffuseColor.rgb=thyroidFollicleColor(diffuseColor.rgb,thyroidDetail);`
        :`#include <color_fragment>\nfloat tissueDetail=tissueField(tissuePoint/${appearance.scale.toPrecision(8)});diffuseColor.rgb*=${appearance.pattern==='gallbladder-folds'?'.72+.38*tissueDetail':'.92+.16*tissueDetail'};`);
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
  return {material, tissueFrame, fibreWeights};
}

/** Opposite owners of a shared tissue chart carry opposite normal signs.
 * A shared material may also draw an ordinary mesh without that attribute;
 * WebGL supplies zero there, which must preserve its original normal. */
export function applyPairedNormalSigns(material:MeshStandardMaterial) {
  const base=material.onBeforeCompile,key=material.customProgramCacheKey.bind(material);
  material.onBeforeCompile=(shader,renderer)=>{
    base(shader,renderer);
    shader.vertexShader='attribute float normalSign;\n'+shader.vertexShader;
    shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\nobjectNormal*=normalSign<0.0?-1.0:1.0;');
  };
  material.customProgramCacheKey=()=>key()+':paired-normal-signs-v2';
}
