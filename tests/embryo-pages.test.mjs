import assert from 'node:assert/strict';
import test from 'node:test';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {gunzipSync} from 'node:zlib';
import {build} from 'esbuild';
const visibility=await build({entryPoints:['app/viewer-state.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {partVisible}=await import('data:text/javascript;base64,'+Buffer.from(visibility.outputFiles[0].contents).toString('base64'));
const depthModule=await build({entryPoints:['app/depth-control.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {setDepthPosition}=await import('data:text/javascript;base64,'+Buffer.from(depthModule.outputFiles[0].contents).toString('base64'));
const directory='public/models/embryo-3month/';
test('Pages embryo retains the source collection, credits and self-contained geometry',async()=>{
 const atlas=JSON.parse(await readFile(directory+'atlas.json')),manifest=JSON.parse(await readFile(directory+'release-manifest.json'));
 assert.equal(atlas.coordinateFrame.landmark,'Spinal axis at the level of the umbilicus');
 assert.ok(atlas.coordinateFrame.sectionIntersections>0);
 const spine=atlas.parts.find(p=>p.id===atlas.coordinateFrame.sourceSpine);
 for(let axis=0;axis<3;axis++)assert.ok(atlas.coordinateFrame.origin[axis]>=spine.bounds[0][axis]&&atlas.coordinateFrame.origin[axis]<=spine.bounds[1][axis],'Origin lies within the source spinal envelope');
 assert.equal(atlas.displayName,'Embryo 3month');assert.equal(atlas.parts.length,174);
 assert.equal(manifest.sourcePreservation.sourceMeshes,166);assert.equal(manifest.sourcePreservation.sourceTriangles,1072606);assert.equal(manifest.sourcePreservation.allOriginalFacesRetained,true);
 for(const id of ['CS23:9226:48','CS23:9226:49','HRA:VH_F_amnion','HRA:VH_F_umbilical_cord'])assert.equal(atlas.parts.find(p=>p.id===id).groups[0],'Skin','All child coverings share the same Skin control');
 const internal=atlas.parts.find(p=>p.system==='skeletal'),skin=atlas.parts.find(p=>p.id==='CS23:9226:48');
 const state={selected:[],visible:[...new Set(atlas.parts.map(p=>p.system))],hidden:[],skinOpacity:1,explode:0,isolate:false,view:'front',rotate:false,reset:0};
 assert.equal(partVisible(internal,state,atlas.parts),false,'Intact Skin hides internal context');
 assert.equal(partVisible(internal,{...state,hidden:atlas.parts.filter(p=>p.groups?.[0]==='Skin').map(p=>p.id)},atlas.parts),true,'Hiding the unified Skin branch reveals anatomy');
 assert.equal(partVisible(internal,{...state,hidden:[skin.id]},atlas.parts),true,'Hiding body Skin alone reveals anatomy even with other coverings checked');
 assert.deepEqual(skin.sourceGroups,['Integumentary system'],'Retain source classification separately from display grouping');
 const cover=atlas.parts.find(p=>p.id==='HRA:VH_F_umbilical_cord');
 const vessels=atlas.parts.filter(p=>['HRA:VH_F_umbilical_artery_1','HRA:VH_F_umbilical_artery_2','HRA:VH_F_umbilical_vein','CS23:9226:157','CS23:9226:159'].includes(p.id));
 assert.equal(vessels.length,5);
 for(const vessel of vessels){assert.ok(['arterial','venous'].includes(vessel.system));assert.equal(vessel.groups[0],'Cardiovascular system');}
 assert.equal(atlas.parts.find(p=>p.id==='HRA:VH_F_placenta_vessels').groups[0],'Cardiovascular system');
 const uncovered=setDepthPosition(state,1/21,atlas.parts);
 assert.equal(partVisible(cover,uncovered,atlas.parts),false,'First depth stage removes cord covering');
 for(const vessel of vessels)assert.equal(partVisible(vessel,uncovered,atlas.parts),true,'Cord vessels remain visible after removing covering');
 const vesselsPeeled=setDepthPosition(uncovered,12/21,atlas.parts);
 for(const vessel of vessels)assert.equal(partVisible(vessel,vesselsPeeled,atlas.parts),false,'Later vessel stage removes internal cord vessels');
 assert.equal(atlas.parts.filter(p=>p.id.startsWith('CS23:')).reduce((n,p)=>n+p.indexCount/3,0),1072606);
 for(const record of manifest.files){const bytes=await readFile('public/'+record.path);assert.equal(createHash('sha256').update(bytes).digest('hex'),record.sha256);}
 for(const part of atlas.parts){const chunk=atlas.chunks[part.chunk];assert.ok(chunk.url.startsWith('/models/embryo-3month/'));
  const raw=gunzipSync(await readFile('public'+chunk.url));assert.ok(gunzipSync(await readFile('public'+chunk.gzip)).equals(raw));assert.equal(raw.length,chunk.bytes);
  for(let i=0;i<part.indexCount;i++)assert.ok(raw.readUInt32LE(part.indices+i*4)<part.vertexCount);
 }
 for(const file of ['ATTRIBUTION.md','public/ATTRIBUTION.md',directory+'ATTRIBUTION.md']){
  const text=await readFile(file,'utf8');assert.ok(text.includes('CC BY-NC-ND 4.0'));assert.ok(text.includes('Carnegie stage 23'));assert.ok(text.includes('University of Amsterdam'));assert.ok(text.includes('CC BY 4.0'));
 }
});
