import * as T from 'three';
export interface HairField {url:string;sha256:string;count:number;positions:number;normals:number;uvs:number;length:number;diameter:number;spacing:number;bytes:number;gzipBytes:number}
export function metricHairGeometry(field:HairField,buffer:ArrayBuffer,index:number):T.InstancedBufferGeometry{
 if(buffer.byteLength!==field.bytes||field.count<=0||field.length<=0||field.diameter<=0)throw Error('Invalid metric hair field');
 const prototype=new T.CylinderGeometry(field.diameter/2,field.diameter/2,field.length,6,1,false);prototype.translate(0,field.length/2,0);
 const g=new T.InstancedBufferGeometry();for(const name of ['position','normal','uv'])g.setAttribute(name,prototype.getAttribute(name));g.setIndex(prototype.getIndex());
 g.setAttribute('hairRoot',new T.InstancedBufferAttribute(new Float32Array(buffer,field.positions,field.count*3),3));
 g.setAttribute('hairRootNormal',new T.InstancedBufferAttribute(new Int16Array(buffer,field.normals,field.count*3),3,true));
 g.setAttribute('hairRootUV',new T.InstancedBufferAttribute(new Float32Array(buffer,field.uvs,field.count*2),2));
 g.setAttribute('partIndex',new T.BufferAttribute(new Float32Array(g.getAttribute('position').count).fill(index),1));g.instanceCount=field.count;
 return g;
}
/** Place a shared, UV-authored closed shaft before the viewer's pose transform. */
export function applyMetricHairMaterial(material:T.MeshStandardMaterial){
 const compile=material.onBeforeCompile,key=material.customProgramCacheKey;material.customProgramCacheKey=()=>`${key.call(material)}:metric-hair-v1`;
 material.onBeforeCompile=(shader,renderer)=>{
  compile.call(material,shader,renderer);
  shader.vertexShader=`attribute vec3 hairRoot;attribute vec3 hairRootNormal;attribute vec2 hairRootUV;varying vec2 hairStrandUV;varying vec2 scalpRootUV;
   vec3 hairVector(vec3 v){vec3 n=normalize(hairRootNormal);vec3 helper=abs(n.y)<.9?vec3(0.,1.,0.):vec3(1.,0.,0.);vec3 t=normalize(cross(helper,n));return t*v.x+n*v.y+cross(n,t)*v.z;}\n`+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>','#include <begin_vertex>\n transformed=hairRoot+hairVector(position);hairStrandUV=uv;scalpRootUV=hairRootUV;');
  shader.vertexShader=shader.vertexShader.replace('#include <beginnormal_vertex>','#include <beginnormal_vertex>\n objectNormal=hairVector(normal);');
  shader.fragmentShader='varying vec2 hairStrandUV;varying vec2 scalpRootUV;\n'+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <map_fragment>',`#include <map_fragment>
   float rootVariation=.9+.1*sin(dot(scalpRootUV,vec2(71.,53.)));float shaftFiber=.97+.03*cos(hairStrandUV.x*37.699);diffuseColor.rgb*=rootVariation*shaftFiber;`);
 };
}
export class MetricHairPicking {
 private cells=new Map<string,number[]>();private roots:Float32Array;private normals:Int16Array;private cell=.002;private base=new T.Vector3();private tip=new T.Vector3();private localRay=new T.Ray();private localPoint=new T.Vector3();private inverse=new T.Matrix4();
 constructor(private field:HairField,buffer:ArrayBuffer){this.roots=new Float32Array(buffer,field.positions,field.count*3);this.normals=new Int16Array(buffer,field.normals,field.count*3);}
 hit(ray:T.Ray,point:T.Vector3,matrix:T.Matrix4):boolean{
  if(!this.cells.size)for(let i=0;i<this.field.count;i++){const j=i*3,key=[0,1,2].map(k=>Math.floor(this.roots[j+k]/this.cell)).join(':');const ids=this.cells.get(key)??[];ids.push(i);this.cells.set(key,ids);}
  this.inverse.copy(matrix).invert();this.localRay.copy(ray).applyMatrix4(this.inverse);this.localPoint.copy(point).applyMatrix4(this.inverse);const c=this.localPoint.toArray().map(x=>Math.floor(x/this.cell)),radius2=(this.field.diameter/2)**2;
  for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++)for(const i of this.cells.get(`${c[0]+x}:${c[1]+y}:${c[2]+z}`)??[]){const j=i*3;this.base.fromArray(this.roots,j);this.tip.set(this.normals[j],this.normals[j+1],this.normals[j+2]).normalize().multiplyScalar(this.field.length).add(this.base);if(this.localRay.distanceSqToSegment(this.base,this.tip)<=radius2)return true;}
  return false;
 }
}
