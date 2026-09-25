import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/viewer-state.ts','app/depth-layers.ts','app/depth-sort.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/depth-validation'});
const load=async file=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(output=>output.path.endsWith(file)).text).toString('base64')}`);
const {partVisible,readViewUrl,viewUrl}=await load('viewer-state.js');
const {DEPTH_LAYERS,depthLayerFor}=await load('depth-layers.js');
const {createDepthOrder}=await load('depth-sort.js');
assert.ok(DEPTH_LAYERS.length>=20);

for(const file of ['atlas-male-complete','atlas-hra-female','atlas']){
 const atlas=JSON.parse(readFileSync(`public/models/${file}.json`));
 const counts=new Map(DEPTH_LAYERS.map(layer=>[layer.id,0]));
 for(const part of atlas.parts){const id=depthLayerFor(part);assert.ok(counts.has(id),`${part.name} has no depth group`);counts.set(id,counts.get(id)+1);}
 assert.equal([...counts.values()].reduce((sum,count)=>sum+count,0),atlas.parts.length);
}

const atlas=JSON.parse(readFileSync('public/models/atlas-male-complete.json'));
const find=(system,term)=>{const part=atlas.parts.find(part=>part.system===system&&part.name.toLowerCase().includes(term));assert.ok(part,`Missing ${term}`);return part;};
for(const [system,term,expected] of [
 ['venous','basilic vein','superficial-veins'],['muscular','trapezius muscle','superficial-muscles'],
 ['muscular','pectoralis minor','second-muscles'],['muscular','flexor digitorum superficialis','intermediate-muscles'],
 ['muscular','flexor digitorum profundus','deep-muscles'],['muscular','pronator quadratus','deepest-muscles'],
 ['skeletal','sternum','thoracic-bones'],['skeletal','femur','limb-bones'],
 ['skeletal','atlas (c1)','spine'],['skeletal','ethmoid bone','skull'],
])assert.equal(depthLayerFor(find(system,term)),expected,term);

const base={selected:[],visible:['integumentary','muscular','skeletal'],skinOpacity:1,depthHidden:[],hidden:[],explode:0,isolate:false,view:'front',rotate:false,reset:0};
const muscle=find('muscular','trapezius muscle');
assert.equal(partVisible(muscle,base),false,'Opaque skin conceals internal anatomy');
const peeled={...base,depthHidden:['skin','superficial-muscles']};
assert.equal(partVisible(muscle,peeled),false,'Hidden muscle layer is removed');
assert.equal(partVisible(find('skeletal','femur'),peeled),true,'Removing skin reveals deeper layers');
assert.equal(partVisible(muscle,{...peeled,selected:[muscle.id]}),true,'Selection remains visible');
const url=viewUrl('http://localhost:3016/','male-full',peeled);
assert.deepEqual(readViewUrl(new URL(url).search,atlas,base).depthHidden,peeled.depthHidden);
assert.deepEqual(readViewUrl('?peel=3',atlas,base).depthHidden,['skin','investing-fascia','superficial-muscles','second-muscles','intermediate-muscles','deep-muscles','deepest-muscles']);

const syntheticPart=(id,x,y,z,system='muscular',groups=['upper limb'])=>({id,name:id,system,groups,bounds:[[x,y,z],[x,y,z]]});
const sortAtlas={parts:[
 syntheticPart('near wrist bone',.6,1.01,0,'skeletal'),
 syntheticPart('near arm bone',.3,1.01,0,'skeletal'),
 syntheticPart('upper core bone',0,1.2,0,'skeletal',['trunk']),
 syntheticPart('lower core bone',0,1,0,'skeletal',['trunk']),
]};
const entry=part=>({name:part.name,parts:[part]});
const compare=createDepthOrder(sortAtlas);
const ordered=[
 syntheticPart('deep arm',.3,1.01,.02),
 syntheticPart('superficial arm',.3,1.01,.1),
 syntheticPart('distal arm',.6,1.01,.02),
].map(entry).sort(compare).map(item=>item.name);
assert.deepEqual(ordered,['superficial arm','distal arm','deep arm'],'Depth precedes distal position');
const coreOrdered=[syntheticPart('upper core',0,1.2,.02,'muscular',['trunk']),syntheticPart('lower core',0,1,.02,'muscular',['trunk'])].map(entry).sort(compare).map(item=>item.name);
assert.deepEqual(coreOrdered,['lower core','upper core'],'Inferior position breaks otherwise equal ties');
console.log('Depth categories cover all anatomy, classify key structures, sort by depth and position, and preserve visibility in shared links.');
