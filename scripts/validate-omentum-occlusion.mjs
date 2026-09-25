/** Check that the combined import does not overlay a second bowel model. */
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import path from 'node:path';
import * as THREE from 'three';

const filename=path.resolve(process.argv[2]??'public/models/atlas-male-complete.json');
const atlas=JSON.parse(readFileSync(filename,'utf8'));
const buffers=new Map();
function mesh(part){
 let buffer=buffers.get(part.chunk);
 if(!buffer){
  buffer=readFileSync(path.join(path.dirname(filename),path.basename(atlas.chunks[part.chunk].url)));
  buffers.set(part.chunk,buffer);
 }
 const geometry=new THREE.BufferGeometry();
 geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(buffer.buffer,buffer.byteOffset+part.positions,part.vertexCount*3),3));
 geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(buffer.buffer,buffer.byteOffset+part.indices,part.indexCount),1));
 if(part.sourceOffset)geometry.translate(...part.sourceOffset);
 geometry.computeBoundingSphere();
 return new THREE.Mesh(geometry,new THREE.MeshBasicMaterial({side:THREE.DoubleSide}));
}
const omentum=atlas.parts.find(part=>part.id==='ZA:Greater omentum');
const jejunum=atlas.parts.find(part=>part.id==='ZA:Jejunum');
assert.ok(omentum&&jejunum,'Missing Z-Anatomy omentum or jejunum');
assert.equal(atlas.parts.filter(part=>part.id.startsWith('FJ')&&/ileum|ileocecal junction|mesentery of small intestine/i.test(part.name)).length,0,
 'The separate BodyParts3D bowel representation must not overlay the Z-Anatomy jejunum');
assert.equal(atlas.parts.filter(part=>part.id.startsWith('FJ')&&/pancreatic duct tree|deferent duct|superficial dorsal vein of penis/i.test(part.name)).length,0,
 'BodyParts3D ducts and dorsal veins must not overlay the detailed Z-Anatomy equivalents');
assert.equal(atlas.parts.filter(part=>part.id.startsWith('O3M:')&&part.system==='skeletal'&&!/Sesamoid bones of hand/i.test(part.name)).length,0,
 'Open 3D Model hand bones must not overlay the detailed Z-Anatomy bones');
assert.equal(atlas.parts.filter(part=>part.id.startsWith('O3M:')&&/Brachiocephalic artery|Arm superficial vein-(Basilic|Cephalic|Median antebrachial|Median cubital) vein/i.test(part.name)).length,0,
 'Open 3D Model arm vessels must not overlay the detailed Z-Anatomy equivalents');
const cover=mesh(omentum),bowel=mesh(jejunum);
const ray=new THREE.Raycaster(),direction=new THREE.Vector3(0,0,-1);
let overlap=0,protrusions=0;
for(let y=.89;y<1.10;y+=.005)for(let x=-.095;x<.105;x+=.005){
 ray.set(new THREE.Vector3(x,y,.3),direction);
 const surface=ray.intersectObject(cover,false)[0];
 if(!surface)continue;
 const hit=ray.intersectObject(bowel,false)[0];
 if(!hit)continue;
 overlap++;
 if(hit.point.z>surface.point.z+.002)protrusions++;
}
assert.ok(overlap>500,`Too little omentum/jejunum overlap sampled: ${overlap}`);
assert.equal(protrusions,0,`${protrusions} bowel samples protrude through the greater omentum`);
console.log(`Validated ${overlap} omentum/jejunum overlap samples; no anterior protrusions over 2 mm or duplicate ileal meshes.`);
