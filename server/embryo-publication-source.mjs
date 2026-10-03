import {readFile} from 'node:fs/promises';
import {resolve} from 'node:path';
import {createHash} from 'node:crypto';
import {gzipSync} from 'node:zlib';

/** Preserve all CS23 source faces; HRA objects remain separately selectable. */
export async function originalEmbryoFiles(root,atlas,put,prefix){
 const placentalVessels=atlas.parts.find(p=>p.id==='HRA:VH_F_placenta_vessels');
 if(placentalVessels){placentalVessels.sourceGroups??=[...placentalVessels.groups];placentalVessels.system='cardiac';placentalVessels.groups=['Cardiovascular system','Placental vessels'];}
 const studyBytes=await readFile(resolve(root,'.local-models/embryo-cs23/source-study.json'));
 const plan=JSON.parse(await readFile(resolve(root,'reports/embryo-cs23/source-assembly-plan.json'),'utf8'));
 if(createHash('sha256').update(studyBytes).digest('hex')!==plan.diagnosticSha256)throw Error('Source study does not match recorded pre-cut source');
 const study=JSON.parse(studyBytes);
 const assembly=JSON.parse(await readFile(resolve(root,'.local-models/embryo-cs23/assembly.json'),'utf8'));
 const exemplar=study.objects.find(o=>o.name==='Notochord'),joined=assembly.objects.find(o=>o.node===exemplar.node);
 const signs=[-1,-1,1],translation=joined.v[0].map((v,a)=>v-3*signs[a]*exemplar.v[0][a]);
 const cloud=assembly.context.slice(0,2).flatMap(o=>o.v),min=[Infinity,Infinity,Infinity],max=[-Infinity,-Infinity,-Infinity];
 for(const v of cloud)for(let a=0;a<3;a++){min[a]=Math.min(min[a],v[a]);max[a]=Math.max(max[a],v[a]);}
 const anchor=min.map((v,a)=>(v+max[a])/2),shift=anchor.map((v,a)=>(2*v+translation[a])/3);
 // Model-local (0,0,0): spinal axis at the umbilical attachment level.
 // This is a coordinate reference, not a deformation of any source tissue.
 const attachment=JSON.parse(await readFile(resolve(root,'reports/embryo-cs23/attachment-assembly.json'),'utf8'));
 const attachmentPoint=attachment.sourceCommonRotationCenter.map((v,a)=>v+attachment.sourceCommonTranslation[a]);
 const spine=study.objects.find(o=>o.name==='Spinal cord');
 const extents=[0,1,2].map(a=>Math.max(...spine.v.map(v=>v[a]))-Math.min(...spine.v.map(v=>v[a]))),levelAxis=extents.indexOf(Math.max(...extents));
 const level=(attachmentPoint[levelAxis]+2*anchor[levelAxis])/3;
 const points=[];
 for(const face of spine.f)for(let edge=0;edge<3;edge++){
  const a=spine.v[face[edge]].map((v,k)=>v*signs[k]+shift[k]),b=spine.v[face[(edge+1)%3]].map((v,k)=>v*signs[k]+shift[k]);
  if((a[levelAxis]-level)*(b[levelAxis]-level)>0||a[levelAxis]===b[levelAxis])continue;
  const t=(level-a[levelAxis])/(b[levelAxis]-a[levelAxis]);if(t>=0&&t<=1)points.push(a.map((v,k)=>v+t*(b[k]-v)));
 }
 if(!points.length)throw Error('Umbilical-level spinal section does not intersect the source');
 const origin=[0,0,0];origin[levelAxis]=level;for(const axis of [0,1,2].filter(axis=>axis!==levelAxis))origin[axis]=(Math.min(...points.map(p=>p[axis]))+Math.max(...points.map(p=>p[axis])))/2;
 atlas.coordinateFrame={origin,landmark:'Spinal axis at the level of the umbilicus',method:'Center of the source spinal-cord section at the recorded cord-attachment level',status:'inferred source landmark; individual umbilical/body level correspondence remains provisional',sourceSpine:'CS23:9226:181',levelAxis,sectionIntersections:points.length};
 let preservedTriangles=0,maxFitError=0;const records=[];
 for(const original of study.objects){
  const part=atlas.parts.find(p=>p.id==='CS23:9226:'+original.node);if(!part)throw Error('Missing original source mesh');
  if(['Skin','Nipples'].includes(original.name)){part.sourceGroups=[...original.groups];part.groups=['Skin'];}
  const current=assembly.objects.find(o=>o.node===original.node);
  if(!['Skin','Umbilical arteries','Umbilical vein'].includes(original.name)){
   if(current.f.length!==original.f.length)throw Error('Unexpected altered source faces');
   for(let i=0;i<original.v.length;i++)for(let a=0;a<3;a++)maxFitError=Math.max(maxFitError,Math.abs(current.v[i][a]-(original.v[i][a]*3*signs[a]+translation[a])));
  }
  const count=original.v.length,indicesOffset=Math.ceil(count*18/4)*4,raw=Buffer.alloc(indicesOffset+original.f.length*12);
  const low=[Infinity,Infinity,Infinity],high=[-Infinity,-Infinity,-Infinity];
  for(let i=0;i<count;i++)for(let a=0;a<3;a++){
   const v=original.v[i][a]*signs[a]+shift[a];raw.writeFloatLE(v,(i*3+a)*4);const packed=raw.readFloatLE((i*3+a)*4);low[a]=Math.min(low[a],packed);high[a]=Math.max(high[a],packed);
   raw.writeInt16LE(Math.max(-32767,Math.min(32767,Math.round(original.n[i][a]*signs[a]*32767))),count*12+(i*3+a)*2);
  }
  for(let i=0;i<original.f.length;i++)for(let a=0;a<3;a++)raw.writeUInt32LE(original.f[i][a],indicesOffset+(i*3+a)*4);
  const hash=createHash('sha256').update(raw).digest('hex'),name=`source-${original.node}-${hash.slice(0,12)}.bin`,compressed=gzipSync(raw,{mtime:0});
  await put(prefix+name+'.gz',compressed);
  atlas.chunks[part.chunk]={url:'/'+prefix+name+'.gz',bytes:raw.length,gzip:'/'+prefix+name+'.gz',gzipBytes:compressed.length};
  Object.assign(part,{positions:0,normals:count*12,indices:indicesOffset,vertexCount:count,indexCount:original.f.length*3,bounds:[low,high]});
  const sourceColor=original.material?.pbrMetallicRoughness?.baseColorFactor??[.8,.7,.65,1];
  atlas.materials[part.material]={color:sourceColor.slice(0,3).map(v=>Math.round(v*255)),roughness:original.material?.pbrMetallicRoughness?.roughnessFactor??.8,metalness:original.material?.pbrMetallicRoughness?.metallicFactor??0,opacity:sourceColor[3]??1};
  part.provenance={...part.provenance,url:'https://www.3dembryoatlas.com/blank',license:'CC BY-NC-ND 4.0',detail:'Original source topology and material color retained. Browser packing and common rigid placement only. Separate HRA placental context is displayed alongside this source.'};
  if(['Skin','Umbilical arteries','Umbilical vein'].includes(original.name))part.coverageLimitation='Original source terminal surface retained; adjoining HRA context remains a separate object. Shared tissue walls and continuous vessel lumens are not established.';
  preservedTriangles+=original.f.length;records.push({node:original.node,vertices:count,triangles:original.f.length,geometrySha256:hash});
 }
 if(maxFitError>1e-12||preservedTriangles!==1072606)throw Error('Source pose/topology replay failed');
 atlas.triangles=atlas.parts.reduce((n,p)=>n+p.indexCount/3,0);
 return {sourceMeshes:records.length,sourceTriangles:preservedTriangles,allOriginalFacesRetained:true,maxUnchangedPoseFitError:maxFitError,sourceStudySha256:createHash('sha256').update(await readFile(resolve(root,'.local-models/embryo-cs23/source-study.json'))).digest('hex'),records};
}
