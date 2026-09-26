import {build} from 'esbuild';
import {readFileSync} from 'node:fs';
import assert from 'node:assert/strict';
import * as T from 'three';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';

const atlas=JSON.parse(readFileSync('.local-models/reference.json','utf8'));
const bundle=await build({entryPoints:['app/section-cap.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {sectionCapGeometry}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);
const y=Number(process.argv[2]??1.41827),plane=new T.Plane(new T.Vector3(0,-1,0),y);
const rows=[],buffers=new Map();
for(const part of atlas.parts){
 if(part.bounds[0][1]>=y||part.bounds[1][1]<=y)continue;
 const chunk=atlas.chunks[part.chunk];let buffer=buffers.get(part.chunk);
 if(!buffer){buffer=readFileSync('.'+chunk.url.replace('/local-models/','/.local-models/'));buffers.set(part.chunk,buffer);}
 const source=new T.BufferGeometry();
 source.setAttribute('position',new T.BufferAttribute(new Float32Array(buffer.buffer,buffer.byteOffset+part.positions,part.vertexCount*3),3));
 source.setIndex(new T.BufferAttribute(new Uint32Array(buffer.buffer,buffer.byteOffset+part.indices,part.indexCount),1));
 const thin=/fascia|pleura|peritone|omentum|mesenter|mesocolon|serosa|capsule|membrane/i.test(part.name)||part.system==='integumentary'||part.system==='fascia'||part.system==='attachments';
 const cut=sectionCapGeometry(source,plane,new T.Matrix4(),{thinShell:thin,width:part.system==='integumentary'?.0035:thin?.0012:part.system==='muscular'?.005:.002,outlineWidth:part.system==='skeletal'?.003:undefined,closureWidth:part.system==='skeletal'?.012:part.system==='muscular'?.016:undefined,closureFraction:part.system==='skeletal'||part.system==='muscular'?.3:undefined});
 let fill=0,total=0;
 if(cut){const p=cut.getAttribute('position'),ix=cut.getIndex(),range={start:0,count:cut.userData.fillIndexCount};
  for(let i=0;i<ix.count;i+=3){const a=new T.Vector3().fromBufferAttribute(p,ix.getX(i)),b=new T.Vector3().fromBufferAttribute(p,ix.getX(i+1)),c=new T.Vector3().fromBufferAttribute(p,ix.getX(i+2));const area=b.sub(a).cross(c.sub(a)).length()/2;total+=area;if(range&&i>=range.start&&i<range.start+range.count)fill+=area;}
 }
 rows.push({name:part.name,system:part.system,triangles:part.indexCount/3,fill:+fill.toFixed(6),total:+total.toFixed(6),thin,open:cut?.userData.openContours??0,closed:cut?.userData.closedContours??0});
 source.dispose();cut?.dispose();
}
const systems={};for(const row of rows){const s=systems[row.system]??={parts:0,filled:0,ribbonOnly:0,none:0,fillArea:0,totalArea:0};s.parts++;if(row.fill)s.filled++;else if(row.total)s.ribbonOnly++;else s.none++;s.fillArea+=row.fill;s.totalArea+=row.total;}
const groupedBones=[];
for(const concept of atlas.concepts){const pieces=atlas.parts.filter(p=>p.conceptId===concept.id&&p.system==='skeletal'&&p.bounds[0][1]<y&&p.bounds[1][1]>y);if(pieces.length<2)continue;
 const sources=pieces.map(part=>{const buffer=buffers.get(part.chunk),source=new T.BufferGeometry();source.setAttribute('position',new T.BufferAttribute(new Float32Array(buffer.buffer,buffer.byteOffset+part.positions,part.vertexCount*3),3));source.setIndex(new T.BufferAttribute(new Uint32Array(buffer.buffer,buffer.byteOffset+part.indices,part.indexCount),1));return source;});
 const combined=mergeGeometries(sources);const cap=sectionCapGeometry(combined,plane,new T.Matrix4(),{outlineWidth:.003,closureWidth:.012,closureFraction:.3});
 groupedBones.push({name:concept.name,parts:pieces.length,fillTriangles:(cap?.userData.fillIndexCount??0)/3,open:cap?.userData.openContours??0});
 for(const source of sources)source.dispose();combined?.dispose();cap?.dispose();
}
if(Math.abs(y-1.41827)<1e-5){
 const vertebra=groupedBones.find(bone=>bone.name==='Vertebra T3');
 assert.ok(vertebra&&vertebra.fillTriangles>0&&vertebra.open===0,'Adjacent bone and cartilage surfaces should form one closed vertebral cap');
 assert.ok(systems.muscular.filled>=60,'Most intersected thoracic muscles should have filled caps');
}
console.log(JSON.stringify({y,parts:rows.length,systems,groupedBones,largestRibbonOnly:rows.filter(r=>!r.fill&&r.total).sort((a,b)=>b.total-a.total).slice(0,30),largestFilled:rows.filter(r=>r.fill).sort((a,b)=>b.fill-a.fill).slice(0,20)},null,2));
