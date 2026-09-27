import assert from 'node:assert/strict';
import {layoutAnatomyLabels} from '../app/label-layout.ts';

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
console.log('Selected anatomy labels use both sides and remain inside nonoverlapping slots.');
