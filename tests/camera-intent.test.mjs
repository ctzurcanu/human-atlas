import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/camera-intent.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {explicitCameraAction}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
test('Selection and display edits do not authorize camera movement; explicit camera actions do',()=>{
 const state={view:'three-quarter',reset:2,focus:0,selected:[],visible:[],hidden:[],explode:0,isolate:false,rotate:false};
 const previous={view:state.view,reset:state.reset,focus:state.focus};
 for(const change of [{selected:['organ']},{selected:[]},{isolate:true},{isolate:false},{hidden:['skin']},{visible:['arterial']},{region:'head-neck'},{depth:.5},{explode:.7},{section:{enabled:true,axis:'sagittal',position:.4}},{sections:[]},{labels:false}])assert.equal(explicitCameraAction({...state,...change},previous,false),false,JSON.stringify(change));
 for(const change of [{focus:1},{reset:3},{view:'back'}])assert.equal(explicitCameraAction({...state,...change},previous,false),true);
 assert.equal(explicitCameraAction(state,previous,true),true,'Explicit saved camera is applied');
});
