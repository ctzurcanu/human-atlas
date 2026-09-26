import {readFileSync,writeFileSync} from 'node:fs';

// Z-Anatomy supplies separate .o* (origin) and .e* (end/insertion) meshes.
// Resolve their measured contact with bone geometry once at import time. A
// bounding-box hit alone is insufficient near a joint or crossing bones.
const atlas=JSON.parse(readFileSync('public/models/atlas-male-complete.json','utf8'));
const bones=atlas.parts.filter(part=>part.system==='skeletal'&&/· Bone(?:-\d+)?$/i.test(part.name));
const markers=atlas.parts.filter(part=>part.system==='attachments'&&/\.([oe])\d*([lr])$/i.test(part.sourceId??part.name));
const chunks=new Map();
function positions(part){
 if(!chunks.has(part.chunk))chunks.set(part.chunk,readFileSync(`public${atlas.chunks[part.chunk].url}`));
 const buffer=chunks.get(part.chunk);
 return new Float32Array(buffer.buffer,buffer.byteOffset+part.positions,part.vertexCount*3);
}
function boxDistance(a,b){
 let squared=0;
 for(let axis=0;axis<3;axis++){
  const gap=Math.max(a.bounds[0][axis]-b.bounds[1][axis],b.bounds[0][axis]-a.bounds[1][axis],0);
  squared+=gap*gap;
 }
 return Math.sqrt(squared);
}
function vertexDistance(a,b){
 const left=positions(a),right=positions(b);
 let best=Infinity;
 for(let i=0;i<left.length;i+=3){
  for(let j=0;j<right.length;j+=3){
   const dx=left[i]-right[j],dy=left[i+1]-right[j+1],dz=left[i+2]-right[j+2];
   const squared=dx*dx+dy*dy+dz*dz;
   if(squared<best)best=squared;
  }
 }
 return Math.sqrt(best);
}
const entries=[];
for(const marker of markers){
 const source=marker.sourceId??marker.name;
 const match=/^(.+)\.([oe])\d*([lr])$/i.exec(source);
 if(!match)continue;
 const side=match[3].toLowerCase()==='l'?'left':'right';
 const candidates=bones.filter(bone=>{
  const boneSide=bone.name.match(/\((left|right)\)/i)?.[1]?.toLowerCase();
  return (!boneSide||boneSide===side)&&boxDistance(marker,bone)<.012;
 });
 const measured=new Map();
 for(const bone of candidates){
  const distance=vertexDistance(marker,bone);
  const prior=measured.get(bone.conceptId);
  if(!prior||distance<prior.distance)measured.set(bone.conceptId,{bone,distance});
 }
 const ranked=[...measured.values()].sort((a,b)=>a.distance-b.distance);
 if(!ranked.length||ranked[0].distance>.0015)continue;
 const matching=ranked.filter(entry=>entry.distance<=Math.min(.0015,ranked[0].distance+.0005));
 entries.push({marker:source,muscle:match[1],side,site:match[2].toLowerCase()==='o'?'origin':'insertion',bones:matching.map(({bone})=>bone.name.replace(/\s*·\s*Bone(?:-\d+)?$/i,''))});
}
entries.sort((a,b)=>a.marker.localeCompare(b.marker));
writeFileSync('app/generated-skeletal-attachments.json',JSON.stringify({source:'Z-Anatomy origin/end marker to bone contact',markerMaxDistanceMm:1.5,entries},null,2)+'\n');
console.log(`Mapped ${entries.length}/${markers.length} muscle attachment markers to bone surfaces.`);
