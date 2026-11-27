import * as T from 'three';
import type {Atlas,Part} from './anatomy';
import {liverControlRegistry} from './liver-controls';
import {liverCompletionRegistry} from './liver-completion-controls';

type Volume={id:string;partId?:string;positions:number;indices:number;vertexCount:number;indexCount:number;origin:[number,number,number];geometryScale?:number;geometrySha256:string};
type Marker={schema:string;bankUrl:string;bankSha256:string;renderOnly:boolean;volumes:Volume[];exclusions:Volume[]};
type Registry={byId:Map<string,Volume>;exclusions:Volume[]};
const cache=new WeakMap<Atlas,Registry|null>();

/** Audited closed volumes classify cut faces only. Their computational mouth
 * disks must never be added to the visible open hilar surface. */
export function liverSectionRegistry(atlas:Atlas):Registry|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;},liver=liverControlRegistry(atlas);
 const m=(atlas as Atlas&{liverSectionIntegration?:Marker}).liverSectionIntegration;
 if(!liver||!m||m.schema!=='liver-section-volumes/v1'||m.renderOnly!==true||m.bankUrl!==liver.marker.bankUrl
  ||m.bankSha256!==liver.marker.bankSha256||!Array.isArray(m.volumes)||![4,7].includes(m.volumes.length)
  ||!Array.isArray(m.exclusions)||m.exclusions.length!==5)return reject();
 const expected=['right','left','quadrate','caudate',...(liver.marker.schema.endsWith('/v2')?['fibrous-appendix','papillary-process','caudate-process']:[])];
 const whole=[...m.volumes,...m.exclusions];
 const integer=(n:number)=>Number.isSafeInteger(n)&&n>=0;
 for(const v of whole){
  if(!v||!integer(v.positions)||v.positions%4||!integer(v.indices)||v.indices%4
   ||!integer(v.vertexCount)||v.vertexCount<4||!integer(v.indexCount)||v.indexCount<12||v.indexCount%3
   ||!Array.isArray(v.origin)||v.origin.length!==3||!v.origin.every(Number.isFinite)
   ||!/^[a-f0-9]{64}$/.test(v.geometrySha256))return reject();
 }
 if(new Set(m.volumes.map(v=>v.id)).size!==expected.length||m.volumes.length!==expected.length||m.volumes.some(v=>!expected.includes(v.id)||v.partId!=='TA98PART:ZA:Liver:'+v.id))return reject();
 const exclusionIds=['artery','portal','biliary','hepatic_vein','retained-native-parents'];
 if(new Set(m.exclusions.map(v=>v.id)).size!==5||m.exclusions.some(v=>!exclusionIds.includes(v.id)))return reject();
 const registry={byId:new Map(m.volumes.map(v=>[v.partId!,v])),exclusions:m.exclusions};cache.set(atlas,registry);return registry;
}

export function liverSectionGeometry(atlas:Atlas,part:Part,buffer:ArrayBuffer,offset:number,shared:(key:string,make:()=>T.BufferAttribute)=>T.BufferAttribute){
 const completion=liverCompletionRegistry(atlas),functional=(atlas as Atlas&{liverCompletionIntegration?:{functionalSections?:{schema:string;bankUrl:string;bankSha256:string;exclusions:Volume[]}}}).liverCompletionIntegration?.functionalSections;
 const validFunctional=!!completion?.functionalPartIds.includes(part.id)&&functional?.schema==='functional-liver-sections/v1'&&atlas.chunks[part.chunk]?.url===functional.bankUrl&&functional.exclusions.length===5;
 const registry=validFunctional?{exclusions:functional!.exclusions}:liverSectionRegistry(atlas);
 const volume=validFunctional?{id:part.id,positions:part.positions,indices:part.indices,vertexCount:part.vertexCount,indexCount:part.indexCount,origin:part.geometryOrigin??[0,0,0]} as Volume:liverSectionRegistry(atlas)?.byId.get(part.id);
 if(!registry||!volume)return undefined;
 const make=(v:Volume)=>{
  if(offset+v.positions+v.vertexCount*12>buffer.byteLength||offset+v.indices+v.indexCount*4>buffer.byteLength)throw new Error('Incomplete liver section volume.');
  const geometry=new T.BufferGeometry(),local=new T.Vector3(...v.origin).sub(new T.Vector3(...(part.geometryOrigin??[0,0,0])));
  // All classifiers in this bank share the organ's exact storage origin.
  if(local.lengthSq()!==0)throw new Error('Mismatched liver section coordinate frame.');
  const scale=v.geometryScale??part.geometryScale??1;
  if(!Number.isFinite(scale)||scale<=0)throw new Error('Invalid liver section scale.');
  geometry.setAttribute('position',shared(`section-position:${v.positions}:${v.vertexCount}:${scale}`,()=>{const points=new Float32Array(buffer,offset+v.positions,v.vertexCount*3);return new T.BufferAttribute(scale===1?points:points.map(x=>x*scale),3);}));
  geometry.setIndex(shared(`section-index:${v.indices}:${v.indexCount}`,()=>new T.BufferAttribute(new Uint32Array(buffer,offset+v.indices,v.indexCount),1)));
  geometry.computeBoundingBox();geometry.computeBoundingSphere();return geometry;
 };
 return {geometry:make(volume),exclusions:registry.exclusions.map(make)};
}
