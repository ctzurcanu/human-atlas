import assert from 'node:assert/strict';
import {build} from 'esbuild';
import {existsSync,readFileSync} from 'node:fs';
import * as THREE from 'three';
import {mergeGeometries} from 'three/examples/jsm/utils/BufferGeometryUtils.js';

const bundle=await build({entryPoints:['app/section-plane.ts','app/viewer-state.ts','app/section-set.ts','app/section-stack.ts','app/section-opacity.ts','app/gastric-region-rim.ts','app/section-assembly.ts','app/stencil-section.ts','app/stencil-caps.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/section-validation'});
const load=async suffix=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(suffix)).text).toString('base64')}`);
const {sectionNormal,sectionPoint,sectionFraction,sectionPlaneEquation,sectionDot}=await load('section-plane.js');
const {readViewUrl,viewUrl,sectionPartVisible}=await load('viewer-state.js');
const {sectionSetBounds}=await load('section-set.js');
const {defaultSection,enabledSections,closeSections,sectionPlaneOnlyChange}=await load('section-stack.js');
const {sectionCapOpacity}=await load('section-opacity.js');
const planeState={selected:[],visible:['digestive'],view:'front',explode:0,reset:0,sections:[{enabled:true,axis:'axial',position:.4,flip:true}]};
assert.ok(sectionPlaneOnlyChange(planeState,{...planeState,sections:[{...planeState.sections[0],position:.5}]}),'dragging only a live plane avoids rebuilding anatomy');
assert.equal(sectionPlaneOnlyChange(planeState,{...planeState,selected:['different']}),false);
assert.equal(sectionPlaneOnlyChange(planeState,{...planeState,sections:[{...planeState.sections[0],enabled:false}]}),false);
assert.equal(defaultSection.axis,'axial');assert.equal(defaultSection.position,.38);
const initialCuts=enabledSections({section:defaultSection});assert.equal(initialCuts.length,0);
const simultaneous=[{...defaultSection,enabled:true},{enabled:true,axis:'sagittal',position:.5,flip:true}];
assert.deepEqual(enabledSections({sections:simultaneous,section:simultaneous[0]}).map(cut=>cut.index),[0,1]);
assert.deepEqual(enabledSections({sections:[simultaneous[0],{...simultaneous[1],enabled:false}],section:simultaneous[0]}).map(cut=>cut.index),[0]);
const camera=[.2,.3,.4,.5,.6,.7,.5,.5,.5,.5,1,.5];
for(const cuts of [[simultaneous[0]],simultaneous]){
 const original={sections:cuts,section:cuts[0],activeSection:cuts.length-1,camera,selected:['stomach'],explode:.3};
 const closed=closeSections(original);
 assert.equal(enabledSections(closed).length,0,'closing the panel disables every cut');
 assert.equal(closed.camera,camera,'closing sections preserves the camera');
 assert.deepEqual(closed.sections,cuts.map(cut=>({...cut,enabled:false})),'closing preserves plane settings');
 assert.equal(closed.section,closed.sections[closed.activeSection],'legacy section state matches the stack');
 assert.ok(cuts.every(cut=>cut.enabled),'closing does not mutate the original state');
}
assert.equal(enabledSections(closeSections({section:simultaneous[0],camera})).length,0,'legacy cuts close too');
assert.ok(!readFileSync('app/scene.tsx','utf8').includes('guideFill'),'the main scene has no section-plane guide');
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
assert.equal(sectionCapOpacity(skin,{selected:[],contextOpacity:.4,skinOpacity:.25}),.25,'unselected view uses skin opacity without ghosting');
assert.equal(sectionCapOpacity(humerus,{selected:[],contextOpacity:.2,skinOpacity:.25},.5),.5,'unselected view preserves source opacity without ghosting');
assert.equal(sectionCapOpacity(humerus,{selected:[humerus.id],contextOpacity:.2,skinOpacity:.25},.5),.5,'selection preserves source translucency');
assert.equal(sectionCapOpacity(skin,{selected:[humerus.id],contextOpacity:.2,skinOpacity:.25}),.05,'selection applies context to unselected skin');
assert.equal(sectionCapOpacity(humerus,{selected:[skin.id],contextOpacity:.2,skinOpacity:1,section:{enabled:true,axis:'axial',position:.38}}),1,'section context stays opaque');
assert.equal(sectionCapOpacity(skin,{selected:[skin.id],contextOpacity:.2,skinOpacity:.25}),.25,'selection preserves skin translucency');
const standard=JSON.parse(readFileSync('public/models/atlas.json','utf8'));
const standardSkin=standard.parts.find(part=>part.id==='FJ2810');assert.ok(standardSkin);
assert.equal(sectionPartVisible(standardSkin,{selected:[],visible:['integumentary'],hidden:[],depthHidden:[],region:'torso',skinOpacity:1,isolate:false}),true,'checked whole-body skin participates in a regional section');
const female=JSON.parse(readFileSync('public/models/atlas-hra-female.json','utf8'));
const femaleSkin=female.parts.find(part=>part.name==='Skin of body');assert.ok(femaleSkin);
assert.equal(sectionPartVisible(femaleSkin,{selected:[],visible:['integumentary'],hidden:[],depthHidden:[],region:'torso',skinOpacity:1,isolate:false}),true,'checked female whole-body skin participates in a regional section');
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
const displayedArmBounds=sectionSetBounds({parts:[humerus,femur,sternum]},armOnly,part=>({min:part.bounds[0].map(v=>v+2),max:part.bounds[1].map(v=>v+2)}));
assert.deepEqual(displayedArmBounds,{min:humerus.bounds[0].map(v=>v+2),max:humerus.bounds[1].map(v=>v+2)},'section range uses displayed anatomy coordinates when exploded');
const base={selected:[],visible:['skeletal','muscular'],explode:0,isolate:false,rotate:false,view:'front',reset:0};
const section={enabled:true,axis:'oblique',azimuth:-137,elevation:42,position:.37,flip:true};
const restored=readViewUrl(new URL(viewUrl('http://localhost:3016/','male-detail',{...base,section})).search,atlas,base);
assert.deepEqual(restored.section,section);
const sectionPair=[{enabled:true,axis:'axial',position:.2,flip:true},{enabled:true,axis:'oblique',azimuth:45,elevation:20,position:.7,flip:false}];
const restoredPair=readViewUrl(new URL(viewUrl('http://localhost:3016/','male-detail',{...base,section:sectionPair[1],sections:sectionPair,activeSection:1})).search,atlas,base);
assert.deepEqual(restoredPair.sections,sectionPair);assert.equal(restoredPair.activeSection,1);assert.deepEqual(restoredPair.section,sectionPair[1]);
const capBundle=await build({entryPoints:['app/section-cap.ts'],bundle:true,platform:'node',format:'esm',write:false,outfile:'/tmp/section-cap-validation.js'});
const capModule=await import(`data:text/javascript;base64,${Buffer.from(capBundle.outputFiles[0].text).toString('base64')}`);
const profileBundle=await build({entryPoints:['app/section-cap-profile.ts'],bundle:true,platform:'node',format:'esm',write:false,outfile:'/tmp/section-cap-profile-validation.js'});
const {sectionCapProfile}=await import(`data:text/javascript;base64,${Buffer.from(profileBundle.outputFiles[0].text).toString('base64')}`);
const tissueBundle=await build({entryPoints:['app/section-tissue.ts'],bundle:true,platform:'node',format:'esm',write:false,outfile:'/tmp/section-tissue-validation.js'});
const {sectionTissue}=await import(`data:text/javascript;base64,${Buffer.from(tissueBundle.outputFiles[0].text).toString('base64')}`);
const topologyBundle=await build({entryPoints:['app/mesh-topology.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {analyzeMeshTopology}=await import(`data:text/javascript;base64,${Buffer.from(topologyBundle.outputFiles[0].text).toString('base64')}`);
assert.equal(sectionCapProfile(standardSkin).outermostOnly,true,'paired standard skin has one exterior cut edge');
assert.equal(sectionCapProfile(femaleSkin).outermostOnly,false,'single female skin surface is not treated as a paired shell');
const tissueCase=(name,system)=>({name,system,bounds:[[0,0,0],[.1,.1,.1]],id:'test'});
for(const [name,system] of [['Jejunum','digestive'],['Stomach','digestive'],['Left atrium','cardiac'],['Nasal cavity','sensory'],['Trachea','respiratory'],['Abdominal aorta','arterial'],['Inferior vena cava','venous']])assert.ok(sectionCapProfile(tissueCase(name,system)).hollowWall,`${name} must have a wall cut and empty lumen`);
for(const [name,system,role] of [['Liver','digestive','liver'],['Superior lobe of right lung','respiratory','lung'],['White matter of spinal cord','nervous','cns'],['Spleen','lymphatic','spleen'],['Rectus abdominis muscle (left)','muscular','muscular']]){
 const part=tissueCase(name,system);assert.equal(sectionCapProfile(part).hollowWall,undefined,`${name} is solid tissue`);assert.equal(sectionTissue(part),role,`${name} cut texture`);
}
const capArea=geometry=>{const positions=geometry.getAttribute('position'),indices=geometry.getIndex();let area=0;for(let i=0;i<indices.count;i+=3){const points=Array.from({length:3},(_,j)=>new THREE.Vector3().fromBufferAttribute(positions,indices.getX(i+j)));area+=new THREE.Vector3().crossVectors(points[1].clone().sub(points[0]),points[2].clone().sub(points[0])).length()/2;}return area;};
const cube=new THREE.BoxGeometry(1,1,1),plane=new THREE.Plane(new THREE.Vector3(0,1,0),0);
const stencilBundle=await build({entryPoints:['app/stencil-caps.ts'],bundle:true,platform:'node',format:'esm',write:false,outfile:'/tmp/stencil-cap-validation.js'});
const {createStencilCaps}=await import(`data:text/javascript;base64,${Buffer.from(stencilBundle.outputFiles[0].text).toString('base64')}`);
const stencilScene=new THREE.Scene(),disabledPlane=new THREE.Plane(new THREE.Vector3(0,1,0),10000);
const stencil=createStencilCaps(stencilScene,[plane,disabledPlane],()=>new THREE.Texture());
const cubePart={id:'cube',name:'Cube',system:'skeletal',bounds:[[-.5,-.5,-.5],[.5,.5,.5]]};
stencil.update([{part:cubePart,index:0,cutIndex:0,geometry:cube,matrix:new THREE.Matrix4(),mode:'solid'}]);
assert.equal(stencil.pick(new THREE.Ray(new THREE.Vector3(0,2,0),new THREE.Vector3(0,-1,0)),()=>true)?.index,0,'GPU solid cap is pickable at its cut face');
assert.equal(stencil.pick(new THREE.Ray(new THREE.Vector3(.8,2,0),new THREE.Vector3(0,-1,0)),()=>true),undefined,'GPU cap does not pick empty space outside the shell');
stencil.update([{part:cubePart,index:1,cutIndex:0,geometry:cube,matrix:new THREE.Matrix4(),mode:'wall',wallWidth:.1}]);
assert.equal(stencil.pick(new THREE.Ray(new THREE.Vector3(0,2,0),new THREE.Vector3(0,-1,0)),()=>true),undefined,'GPU wall cap does not make its lumen selectable');
assert.equal(stencil.pick(new THREE.Ray(new THREE.Vector3(.45,2,0),new THREE.Vector3(0,-1,0)),()=>true)?.index,1,'GPU wall cap is pickable along the tissue wall');
stencil.dispose();
assert.equal(analyzeMeshTopology(cube).closed,true,'a UV-split cube still has a closed welded topology');
const openSheet=new THREE.PlaneGeometry(1,1),sheetTopology=analyzeMeshTopology(openSheet);
assert.ok(sheetTopology.boundaryRatio>.015,'an open sheet has source boundary edges');
assert.equal(sectionCapProfile(humerus,sheetTopology).thinShell,true,'open geometry cannot receive a broad solid cap');
openSheet.dispose();
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
const openTube=new THREE.CylinderGeometry(.2,.2,.6,64,1,true);
const tubeWall=capModule.sectionCapGeometry(openTube,plane,new THREE.Matrix4(),{hollowWall:true,width:.003});
const tubeSolid=capModule.sectionCapGeometry(openTube,plane,new THREE.Matrix4());
assert.ok(tubeWall&&tubeSolid&&tubeWall.userData.fillIndexCount===0,'an open canal needs a wall edge and an open lumen');
assert.ok(capArea(tubeWall)<capArea(tubeSolid)*.1,'an open canal must not fill across its lumen');
tubeWall.dispose();tubeSolid.dispose();openTube.dispose();
const outerWall=new THREE.CylinderGeometry(.24,.24,.6,64,1,true),innerWall=new THREE.CylinderGeometry(.16,.16,.6,64,1,true),walledTube=mergeGeometries([outerWall,innerWall]);
const annulus=capModule.sectionCapGeometry(walledTube,plane,new THREE.Matrix4(),{hollowWall:true,width:.003});
assert.ok(annulus&&annulus.userData.fillIndexCount>0,'paired inner and outer wall contours should produce an annular tissue cap');
assert.ok(capArea(annulus)<Math.PI*.24*.24*.75,'the annular tissue cap must preserve its central lumen');
annulus.dispose();walledTube.dispose();outerWall.dispose();innerWall.dispose();
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
for(const [id,y,limit] of [['ZA:Pleura',1.27,.005],['ZA:Greater omentum',1.1,.002],['ZA:Hypochondriac region.l',1.1,.002]]){
 const part=atlas.parts.find(item=>item.id===id);assert.ok(part,id);
 const chunk=readFileSync('public'+atlas.chunks[part.chunk].url),source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(chunk.buffer,chunk.byteOffset+part.positions,part.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(chunk.buffer,chunk.byteOffset+part.indices,part.indexCount),1));
 const cut=capModule.sectionCapGeometry(source,new THREE.Plane(new THREE.Vector3(0,-1,0),y),new THREE.Matrix4(),{thinShell:true,width:.0012});
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
if(existsSync('.local-models/reference.json')){
 const reference=JSON.parse(readFileSync('.local-models/reference.json'));
 const gluteus=reference.parts.find(part=>part.name==='Gluteus maximus muscle (left)');assert.ok(gluteus);
 assert.equal(sectionPartVisible(gluteus,{selected:[],visible:['muscular'],hidden:[],depthHidden:[],region:'torso',skinOpacity:1,isolate:false}),true,'torso sections include muscles crossing the hip boundary');
 const legSurface=reference.parts.find(part=>part.name==='Body surface (derived) · lower-left');assert.ok(legSurface);
 assert.equal(sectionPartVisible(legSurface,{selected:[],visible:['integumentary'],hidden:[],depthHidden:[],region:'torso',skinOpacity:1,isolate:false}),false,'torso sections exclude detached limb surface');
 const node=reference.parts.find(part=>part.name==='Juxta-intestinal mesenteric nodes');assert.ok(node);
 assert.equal(sectionTissue(node),'lymphatic','mesenteric nodes are lymph tissue, not serosa');
 assert.equal(sectionCapProfile(node).thinShell,false,'mesenteric nodes are solid, not covering sheets');
 const skin=reference.parts.find(part=>part.name==='Body surface (derived) · torso');assert.ok(skin);
 const chunk=reference.chunks[skin.chunk],buffer=readFileSync('.'+chunk.url.replace('/local-models/','/.local-models/'));
 const source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(buffer.buffer,buffer.byteOffset+skin.positions,skin.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(buffer.buffer,buffer.byteOffset+skin.indices,skin.indexCount),1));
 if(skin.sourceOffset)source.translate(...skin.sourceOffset);
 for(const y of [1.41827,1.05257]){
  const plane=new THREE.Plane(new THREE.Vector3(0,-1,0),y);
  const shell=capModule.sectionCapGeometry(source,plane,new THREE.Matrix4(),sectionCapProfile(skin));
  assert.ok(shell&&shell.userData.fillIndexCount===0,`Reference skin incorrectly fills the ${y} body interior`);
  assert.ok(capArea(shell)<.02,`Reference skin cut is thicker than its wall at ${y}`);
  if(y===1.05257){const allContours=capModule.sectionCapGeometry(source,plane,new THREE.Matrix4(),{thinShell:true,width:.0035});assert.ok(allContours&&capArea(shell)<capArea(allContours)*.8,'nested reference skin contours made a second cut ring');allContours.dispose();}
  shell.dispose();
 }
 assert.equal(sectionCapOpacity(skin,{selected:[],contextOpacity:1,skinOpacity:.25}),.25,'reference skin wall follows skin opacity');
 source.dispose();
 for(const name of ['Jejunum','Transverse colon','Stomach','Oesophagus']){
  const part=reference.parts.find(item=>item.name===name);assert.ok(part,name);
  const profile=sectionCapProfile(part);assert.ok(profile.hollowWall,`${name} lumen must remain open`);
  if(name==='Stomach')assert.ok(profile.width>=.007,'stomach cut needs a visibly thick wall');
  const data=readFileSync('.'+reference.chunks[part.chunk].url.replace('/local-models/','/.local-models/'));
  const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(data.buffer,data.byteOffset+part.positions,part.vertexCount*3),3));geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(data.buffer,data.byteOffset+part.indices,part.indexCount),1));
  const middle=(part.bounds[0][1]+part.bounds[1][1])/2,cut=capModule.sectionCapGeometry(geometry,new THREE.Plane(new THREE.Vector3(0,-1,0),middle),new THREE.Matrix4(),profile);
  assert.ok(cut,`${name} should have a visible wall at its middle section`);
  cut.dispose();geometry.dispose();
 }
 for(const name of ['Spleen','Suprarenal gland (left)','Suprarenal gland (right)']){
  const part=reference.parts.find(item=>item.name===name);assert.ok(part,name);
  const data=readFileSync('.'+reference.chunks[part.chunk].url.replace('/local-models/','/.local-models/'));
  const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(data.buffer,data.byteOffset+part.positions,part.vertexCount*3),3));geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(data.buffer,data.byteOffset+part.indices,part.indexCount),1));
  const middle=(part.bounds[0][1]+part.bounds[1][1])/2,cut=capModule.sectionCapGeometry(geometry,new THREE.Plane(new THREE.Vector3(0,-1,0),middle),new THREE.Matrix4(),sectionCapProfile(part,analyzeMeshTopology(geometry)));
  assert.ok(cut?.userData.fillIndexCount>0,`${name} needs a solid section cap`);
  cut.dispose();geometry.dispose();
 }
 assert.equal(sectionCapProfile(reference.parts.find(part=>part.name==='Pancreas')).hollowWall,undefined,'solid digestive glands retain solid cuts');
 assert.equal(sectionCapProfile(reference.parts.find(part=>part.name==='Anterior papillary muscle of right ventricle')).hollowWall,undefined,'solid cardiac muscles retain solid cuts');
}
console.log('Sections, shared URLs, closed organs, and open-surface cut geometry passed.');

const {stencilShell:portShell}=await load('stencil-section.js');
const openPortTube=new THREE.CylinderGeometry(.01,.01,.1,12,2,true),openPortCountingWall=portShell(openPortTube,.0015,false,true);
assert.equal(openPortCountingWall.userData.repairTriangles,0,'open conduits never receive counting end disks');
assert.equal(openPortCountingWall.userData.portRimTriangles,48,'both 12-edge mouths receive annular rims only');
assert.ok(analyzeMeshTopology(openPortCountingWall).closed,'paired conduit counting wall is closed tissue');
const lumenRay=new THREE.Raycaster(new THREE.Vector3(0,-.06,0),new THREE.Vector3(0,1,0));
assert.equal(lumenRay.intersectObject(new THREE.Mesh(openPortCountingWall,new THREE.MeshBasicMaterial({side:THREE.DoubleSide}))).length,0,'a ray traverses both mouths and the lumen without a barrier');
openPortTube.dispose();openPortCountingWall.dispose();

const {createStencilCaps:conduitCaps}=await load('stencil-caps.js');
for(const [system,name] of [['digestive','Jejunum'],['respiratory','Trachea'],['urinary','Ureter'],['reproductive','Vagina'],['arterial','Artery'],['venous','Vein']]){
 const source=new THREE.CylinderGeometry(.01,.01,.1,12,2,true),scene=new THREE.Scene(),caps=conduitCaps(scene,[new THREE.Plane(new THREE.Vector3(0,1,0),0),new THREE.Plane(new THREE.Vector3(1,0,0),10)],()=>new THREE.Texture());
 const part={id:`test:${name}`,name,system,bounds:[[-.01,-.05,-.01],[.01,.05,.01]]};
 caps.update([{part,index:0,cutIndex:0,geometry:source,matrix:new THREE.Matrix4(),mode:'wall',wallWidth:sectionCapProfile(part).width}]);
 assert.equal(scene.children.filter(mesh=>mesh.geometry.userData.openPortWall).length,4,`${name} receives annular counting across both cuts without a special source ID`);
 caps.dispose();source.dispose();
}
const {gastricRegionRim}=await load('gastric-region-rim.js');
const femalePath='.local-models/ta98-female-runtime.json';
if(existsSync(femalePath)){
 // Retained original regional wall remains a regression fixture after finite
 // tissues replace it in runtime. Validate it from the provenance catalogue.
 const runtimeFemale=JSON.parse(readFileSync(femalePath,'utf8')),female=runtimeFemale.parts.some(p=>p.pairedTissue)?JSON.parse(readFileSync('.local-models/ta98-female.json','utf8')):runtimeFemale,part=female.parts.find(part=>part.id==='TA98FEMALE:REGION:A05.5.01.014');
 const chunk=female.chunks[part.chunk],path=chunk.url.replace('/local-models/','.local-models/');
 const bytes=readFileSync(path),buffer=bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength),source=new THREE.BufferGeometry();
 source.setAttribute('position',new THREE.BufferAttribute(new Float32Array(buffer,part.positions,part.vertexCount*3),3));
 source.setIndex(new THREE.BufferAttribute(new Uint32Array(buffer,part.indices,part.indexCount),1));
 const duodenum=female.parts.find(p=>p.id==='LOCAL:female:intestine_duodenum:183'),raw=readFileSync(female.chunks[duodenum.chunk].url.replace('/local-models/','.local-models/')),data=raw.buffer.slice(raw.byteOffset,raw.byteOffset+raw.byteLength),duodenalSource=new THREE.BufferGeometry();
 duodenalSource.setAttribute('position',new THREE.BufferAttribute(new Float32Array(data,duodenum.positions,duodenum.vertexCount*3),3));duodenalSource.setIndex(new THREE.BufferAttribute(new Uint32Array(data,duodenum.indices,duodenum.indexCount),1));
 const originalPositions=Array.from(duodenalSource.getAttribute('position').array),originalIndices=Array.from(duodenalSource.index.array),width=sectionCapProfile(duodenum).width,wall=portShell(duodenalSource,width,false,true);
 assert.equal(width,.0015);assert.equal(wall.index.count/3,1200);assert.equal(wall.userData.portRimTriangles,48);assert.equal(wall.userData.repairTriangles,0);
 const top=analyzeMeshTopology(wall);assert.ok(top.closed);assert.equal(top.duplicateTriangles,0);assert.equal(top.degenerateTriangles,0);
 assert.deepEqual(Array.from(duodenalSource.getAttribute('position').array),originalPositions);assert.deepEqual(Array.from(duodenalSource.index.array),originalIndices);
 const {createStencilCaps}=await load('stencil-caps.js'),scene=new THREE.Scene(),gpu=createStencilCaps(scene,[new THREE.Plane(new THREE.Vector3(0,1,0),-1.186),new THREE.Plane(new THREE.Vector3(1,0,0),10)],()=>new THREE.Texture()),item={part:duodenum,index:183,cutIndex:0,geometry:duodenalSource,matrix:new THREE.Matrix4(),mode:'wall',wallWidth:width};
 gpu.update([item]);const counters=scene.children.filter(mesh=>mesh.geometry.userData.openPortWall);assert.equal(counters.length,4,'two counting passes per plane, without inset fan subtraction');
 assert.ok(gpu.planeError()<1e-12,'cap lies exactly on the cut plane without a geometric gap');
 const counted=counters[0].geometry;gpu.update([{...item,cutIndex:1}]);assert.ok(gpu.planeError()<1e-12,'second cap stays coincident while the plane changes');assert.equal(scene.children.find(mesh=>mesh.geometry.userData.openPortWall).geometry,counted,'cut changes reuse the paired counting wall');
 gpu.dispose();wall.dispose();duodenalSource.dispose();console.log('Female duodenum: paired 1200-triangle counting wall, annular mouths, immutable source and cached GPU passes.');
 const {createSectionAssemblyCache}=await load('section-assembly.js'),cache=createSectionAssemblyCache();
 const candidate={part,index:female.parts.indexOf(part),cutIndex:0,geometry:source,matrix:new THREE.Matrix4(),mode:'wall',wallWidth:part.sectionWall};
 const isolated=cache.collect([candidate]);assert.equal(isolated.caps[0].mode,'solid','paired pyloric wall counts actual tissue volume, not an inset sliver');
 const {stencilShell}=await load('stencil-section.js'),verified=stencilShell(isolated.caps[0].geometry,0,isolated.caps[0].closedWall);
 assert.equal(verified.getAttribute('position').array,isolated.caps[0].geometry.getAttribute('position').array,'verified wall counting shares immutable positions');
 assert.equal(verified.index.array,isolated.caps[0].geometry.index.array,'verified wall counting shares immutable indices');
 assert.equal(verified.userData.repairTriangles,0,'verified walls receive no invented fan caps');verified.dispose();
 const expectedRim=gastricRegionRim(source);assert.ok(expectedRim,'original pyloric regional fixture has an exposed tissue rim');assert.equal(isolated.caps[0].geometry.index.count,part.indexCount+expectedRim.getAttribute('position').count,'section counting includes only the actual current annular junction rim');expectedRim.dispose();
 const counting=isolated.caps[0].geometry;
 const moved=cache.collect([{...candidate,cutIndex:1}]);assert.equal(moved.caps[0].geometry,counting,'dragging or changing cut planes reuses cached wall geometry');
 const whole=female.parts.filter(p=>p.id.startsWith('TA98FEMALE:REGION:')&&p.sectionAssembly==='TA98:female-stomach-wall').map((part,index)=>{
  const raw=readFileSync(female.chunks[part.chunk].url.replace('/local-models/','.local-models/')),data=raw.buffer.slice(raw.byteOffset,raw.byteOffset+raw.byteLength),geometry=new THREE.BufferGeometry();
  geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(data,part.positions,part.vertexCount*3),3));geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(data,part.indices,part.indexCount),1));
  return {...candidate,part,index,geometry};
 });
 assert.equal(whole.length,4);assert.equal(cache.collect(whole).caps[0].mode,'solid','assembled female stomach counts its existing wall volume and lumen');
 for(let mask=1;mask<16;mask++){
  const subset=whole.filter((_,i)=>mask&(1<<i)),cap=cache.collect(subset).caps[0];
  assert.ok(cap.closedWall,`female gastric subset ${mask} must close only paired tissue boundaries`);
  assert.equal(cap.mode,'solid',`female gastric subset ${mask} counts its native paired wall`);
  const top=analyzeMeshTopology(cap.geometry);assert.ok(top.closed);assert.equal(top.duplicateTriangles,0);assert.equal(top.degenerateTriangles,0);
 }
 const separated=whole.map((item,i)=>({...item,matrix:new THREE.Matrix4().makeTranslation(i*.04,0,0)}));
 const separateCaps=cache.collect(separated);
 assert.equal(separateCaps.caps.length,4,'separated verified regional walls remain live GPU cap candidates');
 assert.equal(separateCaps.covered.size,4,'no separated gastric region falls back to delayed CPU work');
 for(const cap of separateCaps.caps){assert.ok(cap.closedWall);assert.equal(cap.mode,'solid');assert.ok(analyzeMeshTopology(cap.geometry).closed);}
 const separatedGeometry=separateCaps.caps.map(cap=>cap.geometry);
 const shifted=cache.collect(separated.map(item=>({...item,cutIndex:1,matrix:item.matrix.clone().premultiply(new THREE.Matrix4().makeTranslation(.01,0,0))})));
 assert.deepEqual(shifted.caps.map(cap=>cap.geometry),separatedGeometry,'movement and second cut reuse each separated wall geometry');
 console.log('All 15 nonempty female gastric regional subsets count paired tissue, without lumen fan disks.');
 for(const item of whole)item.geometry.dispose();cache.dispose();
 const before=Buffer.from(buffer).toString('base64'),rim=gastricRegionRim(source);assert.ok(rim,'native pyloric wall has paired boundaries');
 assert.equal(rim.getAttribute('position').count,isolated.caps[0].geometry.index.count-part.indexCount,'rim matches the currently approved regional seam, without assuming old edge counts');
 const edges=new Map(),pointKey=(p,i)=>[p.getX(i),p.getY(i),p.getZ(i)].map(v=>Math.round(v*1e7)).join(',');
 for(const geometry of [source,rim]){const p=geometry.getAttribute('position'),ix=geometry.index,count=ix?.count??p.count;
  for(let i=0;i<count;i+=3){const keys=[0,1,2].map(k=>pointKey(p,ix?ix.getX(i+k):i+k));
   for(let k=0;k<3;k++){const a=keys[k],b=keys[(k+1)%3],key=[a,b].sort().join('|'),edge=edges.get(key)??{count:0,winding:0};edge.count++;edge.winding+=a<b?1:-1;edges.set(key,edge);}}
 }
 assert.ok([...edges.values()].every(edge=>edge.count===2&&edge.winding===0),'rim closes wall boundaries with matching winding');
 assert.equal(Buffer.from(buffer).toString('base64'),before,'presentation rims do not modify source geometry');
 rim.dispose();source.dispose();console.log('Native pyloric annular rim: actual approved boundary, no open/nonmanifold edges or winding conflicts.');
}

if(existsSync('.local-models/ta98-male.json')){
 const male=JSON.parse(readFileSync('.local-models/ta98-male.json','utf8'));
 const {createSectionAssemblyCache}=await load('section-assembly.js'),{stencilShell}=await load('stencil-section.js'),cache=createSectionAssemblyCache();
 const chunks=new Map();
 const candidates=male.parts.flatMap((part,index)=>{
  if(part.suppressed||part.sectionAssembly!=='TA98:gastric-tract-wall')return [];
  if(!chunks.has(part.chunk)){const bytes=readFileSync(male.chunks[part.chunk].url.replace('/local-models/','.local-models/'));chunks.set(part.chunk,bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength));}
  const data=chunks.get(part.chunk),geometry=new THREE.BufferGeometry();
  geometry.setAttribute('position',new THREE.BufferAttribute(new Float32Array(data,part.positions,part.vertexCount*3),3));geometry.setIndex(new THREE.BufferAttribute(new Uint32Array(data,part.indices,part.indexCount),1));
  return [{part,index,cutIndex:0,geometry,matrix:new THREE.Matrix4(),mode:'wall',wallWidth:part.sectionWall}];
 });
 assert.equal(candidates.length,12,'male native gastric assembly includes both adjoining conduits and seam pieces');
 const first=cache.collect(candidates).caps[0];assert.ok(first.closedWall);assert.equal(first.mode,'solid');
 const count=stencilShell(first.geometry,0,first.closedWall);assert.equal(count.index.array,first.geometry.index.array);assert.equal(count.getAttribute('position').array,first.geometry.getAttribute('position').array);assert.equal(count.userData.repairTriangles,0);
 const after=cache.collect(candidates.map(item=>({...item,cutIndex:1}))).caps[0];assert.equal(after.geometry,first.geometry,'male wall is reused for the second/moving section');
 console.log(`Male gastric wall: ${count.index.count/3} existing triangles, zero added counting arrays or fan closures.`);
 count.dispose();cache.dispose();for(const candidate of candidates)candidate.geometry.dispose();
}
