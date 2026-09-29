import assert from 'node:assert/strict';
import {labelAreaForPanels,layoutAnatomyLabels} from '../app/label-layout.ts';

const area={left:20,right:600,top:40,bottom:340};
const anchors=Array.from({length:18},(_,index)=>({id:String(index),x:index%2?390:220,y:170+index%3,text:`Selected anatomical structure ${index}`}));
const placed=layoutAnatomyLabels(anchors,area,text=>text.length*6);
assert.equal(placed.length,anchors.length);
assert.equal(placed.filter(label=>label.side==='left').length,9);
assert.equal(placed.filter(label=>label.side==='right').length,9);
for(const label of placed){
 assert.ok(label.left>=area.left&&label.right<=area.right);
 assert.ok(label.top>=area.top&&label.bottom<=area.bottom);
 assert.ok(label.label.length>0);
}
for(let a=0;a<placed.length;a++)for(let b=a+1;b<placed.length;b++){
 const first=placed[a],second=placed[b];
 assert.ok(first.right<=second.left||second.right<=first.left||first.bottom<=second.top||second.bottom<=first.top,`Labels ${a} and ${b} overlap`);
}
assert.deepEqual(layoutAnatomyLabels([],area,()=>0),[]);
for(const width of [640,1280,1920]){
 const panels=[{kind:'panel',left:5,right:205,top:54,bottom:715},{kind:'panel',left:width-205,right:width-5,top:54,bottom:715}];
 const bounds=labelAreaForPanels(width,720,panels),labels=layoutAnatomyLabels([{id:'left',x:width*.4,y:300,text:'Left anatomy'},{id:'right',x:width*.6,y:300,text:'Right anatomy'}],bounds,text=>text.length*6);
 assert.equal(labels.find(label=>label.side==='left').left-panels[0].right,5,'Left labels keep five pixel clearance at every viewport width.');
 assert.equal(panels[1].left-labels.find(label=>label.side==='right').right,5,'Right labels keep five pixel clearance at every viewport width.');
 const menu={kind:'panel',left:width-255,right:width-5,top:80,bottom:200};
 assert.equal(menu.left-labelAreaForPanels(width,720,[...panels,menu]).right,5,'Short popup panels also reserve five pixels.');
}
const mobileBounds=labelAreaForPanels(390,844,[{kind:'top',left:0,right:390,top:0,bottom:97},{kind:'panel',left:5,right:385,top:104,bottom:839}]);
assert.deepEqual(layoutAnatomyLabels(anchors,mobileBounds,text=>text.length*6),[],'Labels cannot fit behind a full width open panel.');
// Hip bone + its first relationship step: 50 pieces, 44 distinct names.
const hipAnchors=Array.from({length:44},(_,index)=>({id:`hip-${index}`,x:300,y:240,text:`Related hip structure ${index}`}));
const hipArea={left:280,right:930,top:82,bottom:606};
const hipLabels=layoutAnatomyLabels(hipAnchors,hipArea,text=>text.length*6);
assert.equal(hipLabels.length,44,'Every hip relationship label must fit, even when most anchors prefer one side.');
for(const label of hipLabels)assert.ok(label.top>=hipArea.top&&label.bottom<=hipArea.bottom);
for(let a=0;a<hipLabels.length;a++)for(let b=a+1;b<hipLabels.length;b++){
 const first=hipLabels[a],second=hipLabels[b];
 assert.ok(first.right<=second.left||second.right<=first.left||first.bottom<=second.top||second.bottom<=first.top,`Hip labels ${a} and ${b} overlap`);
}
console.log('Selected anatomy labels use both sides and remain inside nonoverlapping slots.');
