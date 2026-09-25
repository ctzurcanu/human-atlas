import assert from 'node:assert/strict';
import {build} from 'esbuild';
import {readFileSync} from 'node:fs';
import * as THREE from 'three';

const bundle=await build({entryPoints:['app/section-plane.ts','app/viewer-state.ts','app/section-set.ts','app/section-stack.ts','app/section-opacity.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/section-validation'});
const load=async suffix=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(suffix)).text).toString('base64')}`);
const {sectionNormal,sectionPoint,sectionFraction,sectionPlaneEquation,sectionDot}=await load('section-plane.js');
const {readViewUrl,viewUrl,sectionPartVisible}=await load('viewer-state.js');
const {sectionSetBounds}=await load('section-set.js');
const {defaultSection,enabledSections}=await load('section-stack.js');
const {sectionCapOpacity}=await load('section-opacity.js');
assert.equal(defaultSection.axis,'axial');assert.equal(defaultSection.position,.38);
const initialCuts=enabledSections({section:defaultSection});assert.equal(initialCuts.length,0);
const simultaneous=[{...defaultSection,enabled:true},{enabled:true,axis:'sagittal',position:.5,flip:true}];
assert.deepEqual(enabledSections({sections:simultaneous,section:simultaneous[0]}).map(cut=>cut.index),[0,1]);
assert.deepEqual(enabledSections({sections:[simultaneous[0],{...simultaneous[1],enabled:false}],section:simultaneous[0]}).map(cut=>cut.index),[0]);
const bounds={min:[-1,0,-.5],max:[1,2,.5]};
const combined=enabledSections({sections:simultaneous,section:simultaneous[0]}).map(({section})=>sectionPlaneEquation(bounds,section));
const retained=point=>combined.every(plane=>sectionDot(plane.normal,point)+plane.constant>=0);
assert.equal(retained([.5,.5,0]),true,'a point kept by both sections should remain');
assert.equal(retained([-.5,.5,0]),false,'the sagittal cut should remove its opposite side');
assert.equal(retained([.5,1.5,0]),false,'the axial cut should remove its opposite side');
for(const [axis,expected] of [['axial',[0,-1,0]],['sagittal',[1,0,0]],['coronal',[0,0,-1]]])assert.deepEqual(sectionNormal({axis}),expected);
for(const section of [
 {axis:'axial',position:.2,flip:true},
 {axis:'sagittal',position:.7,flip:false},
 {axis:'coronal',position:.4,flip:true},
 {axis:'oblique',azimuth:37,elevation:-28,position:.63,flip:true},
]){
 const normal=sectionNormal(section);assert.ok(Math.abs(sectionDot(normal,normal)-1)<1e-12);
 const point=sectionPoint(bounds,section),fraction=sectionFraction(bounds,section,point);
 assert.ok(Math.abs(fraction-section.position)<1e-12,section.axis);
 const plane=sectionPlaneEquation(bounds,section);
 assert.ok(Math.abs(sectionDot(plane.normal,point)+plane.constant)<1e-12,section.axis);
 const kept=point.map((value,i)=>value+plane.normal[i]*.01),removed=point.map((value,i)=>value-plane.normal[i]*.01);
 assert.ok(sectionDot(plane.normal,kept)+plane.constant>0);
 assert.ok(sectionDot(plane.normal,removed)+plane.constant<0);
}
const atlas=JSON.parse(readFileSync('public/models/atlas-male-complete.json','utf8'));
assert.equal(readViewUrl('',atlas,{selected:[],visible:['skeletal'],view:'front'}).section.position,.38);
const humerus=atlas.parts.find(part=>part.id==='ZA:Humerus.l:Bone');
const femur=atlas.parts.find(part=>part.id==='ZA:Femur.l:Bone');
const sternum=atlas.parts.find(part=>part.id==='ZA:Body of sternum:Bone-7');
assert.ok(humerus&&femur&&sternum);
const skin=atlas.parts.find(part=>part.system==='integumentary');assert.ok(skin);
assert.equal(sectionCapOpacity(skin,{selected:[],contextOpacity:1,skinOpacity:.25}),.25,'skin cap follows skin opacity');
assert.equal(sectionCapOpacity(skin,{selected:[],contextOpacity:.4,skinOpacity:.25}),.1,'skin cap combines context and skin opacity');
assert.equal(sectionCapOpacity(humerus,{selected:[],contextOpacity:.2,skinOpacity:.25},.5),.1,'all caps follow context and source opacity');
assert.equal(sectionCapOpacity(humerus,{selected:[humerus.id],contextOpacity:.2,skinOpacity:.25},.5),1,'selected caps remain solid like selected surfaces');
const checked={selected:[sternum.id],visible:['skeletal'],hidden:[sternum.id],depthHidden:[],region:'all',skinOpacity:.1,isolate:false};
assert.equal(sectionPartVisible(sternum,checked),false,'selection must not add an unchecked mesh to a section');
assert.equal(sectionPartVisible(humerus,checked),true);
assert.equal(sectionPartVisible(femur,checked),true);
const checkedBounds=sectionSetBounds({parts:[humerus,femur,sternum]},checked);
for(let axis=0;axis<3;axis++){
 assert.equal(checkedBounds.min[axis],Math.min(humerus.bounds[0][axis],femur.bounds[0][axis]));
 assert.equal(checkedBounds.max[axis],Math.max(humerus.bounds[1][axis],femur.bounds[1][axis]));
}
const armOnly={...checked,region:'upper-left',selected:[]};
assert.equal(sectionPartVisible(humerus,armOnly),true);
assert.equal(sectionPartVisible(femur,armOnly),false);
assert.deepEqual(sectionSetBounds({parts:[humerus,femur,sternum]},armOnly),{min:humerus.bounds[0],max:humerus.bounds[1]});
const base={selected:[],visible:['skeletal','muscular'],explode:0,isolate:false,rotate:false,view:'front',reset:0};
const section={enabled:true,axis:'oblique',azimuth:-137,elevation:42,position:.37,flip:true};
const restored=readViewUrl(new URL(viewUrl('http://localhost:3016/','male-detail',{...base,section})).search,atlas,base);
assert.deepEqual(restored.section,section);
const sectionPair=[{enabled:true,axis:'axial',position:.2,flip:true},{enabled:true,axis:'oblique',azimuth:45,elevation:20,position:.7,flip:false}];
const restoredPair=readViewUrl(new URL(viewUrl('http://localhost:3016/','male-detail',{...base,section:sectionPair[1],sections:sectionPair,activeSection:1})).search,atlas,base);
assert.deepEqual(restoredPair.sections,sectionPair);assert.equal(restoredPair.activeSection,1);assert.deepEqual(restoredPair.section,sectionPair[1]);
const capBundle=await build({entryPoints:['app/section-cap.ts'],bundle:true,platform:'node',format:'esm',write:false,outfile:'/tmp/section-cap-validation.js'});
const capModule=await import(`data:text/javascript;base64,${Buffer.from(capBundle.outputFiles[0].text).toString('base64')}`);
const capArea=geometry=>{const positions=geometry.getAttribute('position'),indices=geometry.getIndex();let area=0;for(let i=0;i<indices.count;i+=3){const points=Array.from({length:3},(_,j)=>new THREE.Vector3().fromBufferAttribute(positions,indices.getX(i+j)));area+=new THREE.Vector3().crossVectors(points[1].clone().sub(points[0]),points[2].clone().sub(points[0])).length()/2;}return area;};
const cube=new THREE.BoxGeometry(1,1,1),plane=new THREE.Plane(new THREE.Vector3(0,1,0),0);
const cap=capModule.sectionCapGeometry(cube,plane,new THREE.Matrix4());
assert.ok(cap&&cap.index.count>=6,'cube cut should have a filled cap');
const positions=cap.getAttribute('position');for(let i=0;i<positions.count;i++)assert.ok(Math.abs(positions.getY(i)+.00015)<1e-5);
const boneCut=capModule.sectionCapGeometry(cube,plane,new THREE.Matrix4(),{outlineWidth:.002});
assert.deepEqual(boneCut.groups.map(group=>group.materialIndex),[0,1],'bone cuts have separate cancellous interior and cortical rim');
assert.equal(boneCut.getAttribute('uv').count,boneCut.getAttribute('position').count,'bone texture coordinates cover the whole cut');
boneCut.dispose();
const shifted=capModule.sectionCapGeometry(cube,new THREE.Plane(new THREE.Vector3(0,1,0),-2),new THREE.Matrix4().makeTranslation(0,2,0));
assert.ok(shifted&&shifted.index.count>=6,'translated meshes should cap in world space');
for(let i=0;i<shifted.getAttribute('position').count;i++)assert.ok(Math.abs(shifted.getAttribute('position').getY(i)-1.99985)<1e-5);
shifted.dispose();
cap.dispose();cube.dispose();
const sheet=new THREE.PlaneGeometry(1,1),sheetCut=capModule.sectionCapGeometry(sheet,plane,new THREE.Matrix4());
assert.ok(sheetCut,'an open surface should retain a visible narrow cut edge');
assert.ok(capArea(sheetCut)<.003,'an open surface must not become a broad filled cap');
sheetCut.dispose();sheet.dispose();
const actualCut=(id,y,options={})=>{
 const part=atlas.parts.find(item=>item.id===id);assert.ok(part,id);
 const chunk=readFileSync('public'+atlas.chunks[part.chunk].url),source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),1));
 if(part.sourceOffset)source.translate(...part.sourceOffset);
 const cut=capModule.sectionCapGeometry(source,new THREE.Plane(new THREE.Vector3(0,-1,0),y),new THREE.Matrix4(),options);
 source.dispose();assert.ok(cut,id+' should intersect the section');return cut;
};
for(const [id,y,minArea] of [
 ['ZA:Radius.l:Bone',1.058,.00015],
 ['ZA:Twelfth rib.l:Bone-4',1.058,.00001],
 ['ZA:Body of sternum:Bone-7',1.27,.0002],
 ['ZA:Vertebra T8:Bone-5',1.27,.0004],
 ['ZA:Vastus lateralis muscle.l:Extension',.603,.001],
 ['ZA:Vastus medialis muscle.l:Extension',.603,.0005],
 ['ZA:Adductor magnus.l:Adductor',.603,.001],
 ['ZA:Descending part of trapezius muscle.l:External rotation',1.27,.0001],
 ['ZA:External abdominal oblique muscle.l',1.058,.001],
]){
 const bone=id.includes(':Bone');
 const cut=actualCut(id,y,{outlineWidth:.003,closureWidth:bone?.012:.016,closureFraction:.3});
 assert.ok(cut.groups.some(group=>group.materialIndex===0),id+' needs a filled cut face');
 const fill=cut.groups.find(group=>group.materialIndex===0),indices=cut.getIndex(),position=cut.getAttribute('position');let filledArea=0;
 for(let i=fill.start;i<fill.start+fill.count;i+=3){const points=Array.from({length:3},(_,j)=>new THREE.Vector3().fromBufferAttribute(position,indices.getX(i+j)));filledArea+=new THREE.Vector3().crossVectors(points[1].sub(points[0]),points[2].sub(points[0])).length()/2;}
 assert.ok(filledArea>minArea,id+' cut face is too small');cut.dispose();
}
const skinCut=actualCut('ZA:Anterior region of thigh.l:Skin-1',.603,{thinShell:true,width:.0035});
assert.ok(capArea(skinCut)>.00015,'skin must retain a visible cut thickness');skinCut.dispose();
for(const [id,axis] of [['O3M:Pectoralis major.r',2],['ZA:Rectus abdominis muscle.l:Flexion',1]]){
 const part=atlas.parts.find(item=>item.id===id);
 assert.ok(part,`Missing source mesh ${id}`);
 const chunk=readFileSync('public'+atlas.chunks[part.chunk].url),source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),1));
 const normal=new THREE.Vector3().setComponent(axis,-1),middle=(part.bounds[0][axis]+part.bounds[1][axis])/2;
 const filled=capModule.sectionCapGeometry(source,new THREE.Plane(normal,middle),new THREE.Matrix4());
 assert.ok(filled?.index?.count>=6,`Open publisher mesh should have a filled ${id} cut`);
 filled.dispose();source.dispose();
}
const boundaryEdges=part=>{
 const chunk=readFileSync('public'+atlas.chunks[part.chunk].url),indices=new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),edges=new Map();
 for(let i=0;i<indices.length;i+=3)for(let side=0;side<3;side++){
  const a=indices[i+side],b=indices[i+(side+1)%3],key=a<b?`${a},${b}`:`${b},${a}`;
  edges.set(key,(edges.get(key)??0)+1);
 }
 return [...edges.values()].filter(count=>count===1).length;
};
const stomach=atlas.parts.find(part=>part.id==='ZA:Stomach');
assert.ok(stomach);
assert.ok(boundaryEdges(stomach)<=24,'the selected stomach must not contain the former anterior-wall opening');
for(const id of ['ZA:Pancreas','ZA:Spleen','ZA:Gallbladder','ZA:Kidney.l','ZA:Kidney.r','ZA:Urinary bladder']){
 const part=atlas.parts.find(item=>item.id===id);assert.ok(part,id);
 assert.equal(boundaryEdges(part),0,`${id} has an open surface`);
}
for(const [id,y,limit] of [['ZA:Pleura',1.27,.005],['ZA:Greater omentum',1.1,.002],['ZA:Hypochondriac region.l',1.1,.002],['ZA:Transverse colon',1.1,.003]]){
 const part=atlas.parts.find(item=>item.id===id);assert.ok(part,id);
 const chunk=readFileSync('public'+atlas.chunks[part.chunk].url),source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),1));
 const cut=capModule.sectionCapGeometry(source,new THREE.Plane(new THREE.Vector3(0,-1,0),y),new THREE.Matrix4(),{thinShell:id!=='ZA:Transverse colon',width:.0012});
 assert.ok(cut,id+' should retain a cut edge');assert.ok(capArea(cut)<limit,id+' generated a false broad cut face');
 cut.dispose();source.dispose();
}
for(const id of ['ZA:Kidney.l','ZA:Spleen']){
 const part=atlas.parts.find(item=>item.id===id);assert.ok(part,id);
 const chunk=readFileSync('public'+atlas.chunks[part.chunk].url),source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),1));
 const middle=(part.bounds[0][1]+part.bounds[1][1])/2;
 const cut=capModule.sectionCapGeometry(source,new THREE.Plane(new THREE.Vector3(0,-1,0),middle),new THREE.Matrix4());
 assert.ok(cut&&capArea(cut)>.0005,id+' should retain a solid cut face');
 cut.dispose();source.dispose();
}
console.log('Sections, shared URLs, closed organs, and open-surface cut geometry passed.');
