import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const bundled=await build({stdin:{contents:`export {skinEyeHoleAt} from './app/skin-eye-appearance';export {partVisible} from './app/viewer-state';export {isSkinPart} from './app/depth-layers';`,resolveDir:process.cwd(),loader:'ts'},bundle:true,platform:'node',format:'esm',write:false});
const {skinEyeHoleAt,partVisible,isSkinPart}=await import(`data:text/javascript;base64,${Buffer.from(bundled.outputFiles[0].contents).toString('base64')}`);
const atlas=JSON.parse(readFileSync('.local-models/ta98-female-runtime.json'));
const state={selected:[],visible:['integumentary','sensory','nervous'],hidden:[],depthHidden:[],depth:0,skinOpacity:1,context:.18,isolate:false,explode:0,region:'all'};
test('zero-depth skin exposes ocular exterior, hides small eyelashes by default, and retains internal anatomy',()=>{
 const eyes=atlas.parts.filter(p=>/^LOCAL:female:[lr]_iris:/.test(p.id));assert.equal(atlas.parts.filter(p=>p.id.startsWith('LOCAL:female:eyeball:')).length,0);const lashes=atlas.parts.filter(p=>p.id.startsWith('LOCAL:female:eyelashes:'));assert.equal(eyes.length,2);assert.equal(lashes.length,2);
 const current={...state,hidden:atlas.defaultHidden??[]};
 for(const p of eyes)assert.ok(partVisible(p,current,atlas.parts),p.name);
 for(const p of lashes)assert.equal(partVisible(p,current,atlas.parts),false,p.name);
 for(const p of lashes)assert.ok(isSkinPart(p));
 const brain=atlas.parts.find(p=>p.system==='nervous'&&!p.suppressed&&!p.id.startsWith('HRA:'));assert.ok(brain);assert.equal(partVisible(brain,state,atlas.parts),true);
});
test('eye aperture permits passage only within the native eyelid footprint',()=>{
 const p=atlas.parts.find(p=>p.id==='LOCAL:female:skin:1005');const eyes=atlas.materials[p.material].skinEyes;
 for(const [x,y] of eyes.centers){assert.ok(skinEyeHoleAt(eyes,{x,y}));assert.equal(skinEyeHoleAt(eyes,{x:x+eyes.halfWidth*1.01,y}),false);assert.equal(skinEyeHoleAt(eyes,{x,y:y+eyes.halfHeight*1.01}),false);}
});
test('male and female checked internal anatomy remains present at opaque zero-depth skin',()=>{
 for(const path of ['.local-models/ta98-runtime.json','.local-models/ta98-female-runtime.json']){
  const model=JSON.parse(readFileSync(path));const current={...state,visible:[...new Set(model.parts.map(p=>p.system))],hidden:model.defaultHidden??[]};
  const internal=model.parts.filter(p=>!p.suppressed&&!current.hidden.includes(p.id)&&!isSkinPart(p));assert.ok(internal.length>100);
  for(const part of internal)assert.ok(partVisible(part,current,model.parts),`${path}: ${part.name}`);
 }
});
