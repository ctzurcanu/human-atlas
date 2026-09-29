import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';

const bundle=await build({entryPoints:['app/viewer-state.ts','app/depth-layers.ts','app/depth-sort.ts','app/depth-control.ts','app/section-opacity.ts','app/use-view-history.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/depth-validation'});
const load=async file=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(output=>output.path.endsWith(file)).text).toString('base64')}`);
const {partVisible,readViewUrl,viewUrl}=await load('viewer-state.js');
const {DEPTH_LAYERS,depthLayerFor}=await load('depth-layers.js');
const {depthPosition,partLayerOpacity,setDepthPosition}=await load('depth-control.js');
const {sectionCapOpacity}=await load('section-opacity.js');
const {historyReducer}=await load('use-view-history.js');
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
const skin=atlas.parts.find(part=>depthLayerFor(part)==='skin');
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

// Each full interval removes exactly one layer; its midpoint only fades that layer.
const allOn={...base,visible:[...new Set(atlas.parts.map(part=>part.system))]};
let current=allOn;
const n=DEPTH_LAYERS.length;
for(let index=0;index<n;index++){
 const layer=DEPTH_LAYERS[index],part=atlas.parts.find(part=>depthLayerFor(part)===layer.id);
 current=setDepthPosition(current,(index+.5)/n,atlas.parts);
 assert.equal(current.depthHidden.length,index);
 if(part){
  assert.ok(Math.abs(partLayerOpacity(part,current)-.5)<1e-10,layer.name+' fades halfway');
  assert.ok(Math.abs(sectionCapOpacity(part,{...current,sections:[{enabled:true}]} )-.5)<1e-10,'Section walls use the same fade');
  assert.equal(partVisible(part,{...current,selected:[part.id]}),true);
 }
 current=setDepthPosition(current,(index+1)/n,atlas.parts);
 assert.equal(current.depthHidden.length,index+1);
 if(part)assert.equal(partVisible(part,{...current,selected:[part.id]}),false,'A completed layer hides selected tissue too');
}
assert.equal(atlas.parts.filter(part=>partVisible(part,current)).length,0,'The last tick hides the entire body');
for(let index=n-1;index>=0;index--){
 const part=atlas.parts.find(part=>depthLayerFor(part)===DEPTH_LAYERS[index].id);
 // Save/reload at a tick must not turn the rounded prefix into independent hides.
 current=readViewUrl(viewUrl('http://localhost:3016/','male-full',current,undefined,atlas),atlas,base);
 current=setDepthPosition(current,(index+.5)/n,atlas.parts);
 assert.equal(current.depthHidden.length,index,'Reversing restores each completed layer');
 if(part)assert.ok(Math.abs(partLayerOpacity(part,current)-.5)<1e-10,'Reversed opacity is symmetric');
 current=setDepthPosition(current,index/n,atlas.parts);
}
assert.equal(current.depthHidden.length,0);
assert.ok(atlas.parts.every(part=>partLayerOpacity(part,current)===1),'Zero restores fully opaque layers');
const manuallyHidden=setDepthPosition({...allOn,depthHidden:['skull']},.5/n,atlas.parts);
assert.ok(manuallyHidden.depthHidden.includes('skull'),'Independent layer hides are retained');
assert.equal(partLayerOpacity(skin,{...base,skinOpacity:.23}),.23,'Old surface URLs retain opacity');
assert.ok(Math.abs(depthPosition({...base,skinOpacity:.23})-.77/n)<1e-10);
const halfway=setDepthPosition(allOn,1.5/n,atlas.parts),roundTrip=readViewUrl(viewUrl('http://localhost:3016/','male-full',halfway,undefined,atlas),atlas,base);
assert.ok(Math.abs(roundTrip.depth-halfway.depth)<.000001);
assert.deepEqual(roundTrip.depthHidden,halfway.depthHidden);
let h={present:allOn,past:[],future:[],time:0,key:''};
for(let i=1;i<=100;i++)h=historyReducer(h,{type:'set',time:i*10,update:setDepthPosition(h.present,i/100,atlas.parts)});
assert.equal(h.past.length,1,'A continuous depth drag makes one undo step');
console.log('Depth control fades all 21 layers in order and reverses symmetrically, including selected tissue, sections, saved views and undo.');

assert.ok(!new URL(viewUrl('http://localhost:3016/','male-full',halfway,undefined,atlas)).searchParams.has('depth'),'Sequential peeling uses dp alone, without long layer names');
const override={...halfway,depthHidden:[]},overrideReload=readViewUrl(viewUrl('http://localhost:3016/','male-full',override,undefined,atlas),atlas,base);
assert.deepEqual(overrideReload.depthHidden,[],'Manual layer overrides survive sharing');
assert.deepEqual(readViewUrl('?dp=1',atlas,base).depthHidden,DEPTH_LAYERS.map(layer=>layer.id));
