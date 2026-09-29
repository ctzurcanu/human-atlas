import assert from 'node:assert/strict';
import {dampMotion,explosionMarks,explosionToSlider,sliderToExplosion,smoothStep} from '../app/transition-motion.ts';
const near=(a,b)=>assert.ok(Math.abs(a-b)<1e-9,`${a} != ${b}`);
for(const steps of [1,2,3,4,7,17]){
 const marks=explosionMarks(steps);near(marks[0],0);near(marks.at(-1),1);
 assert.ok(marks.every((mark,i)=>i===0||mark>marks[i-1]),'Every hierarchy level retains usable slider space');
 if(steps>4){near(marks[1],2/steps);near(marks[2]-marks[1],2/steps);}
 for(let i=0;i<=1000;i++)near(sliderToExplosion(explosionToSlider(i/1000,steps),steps),i/1000);
 marks.forEach((mark,i)=>near(sliderToExplosion(mark,steps,true),i/steps));
 near(sliderToExplosion(marks[1]-.001,steps,true),1/steps);
 assert.notEqual(sliderToExplosion(marks[1]-.001,steps),1/steps,'Dragging remains precise instead of sticking to a tick');
}
const run=fps=>{
 let value=0,velocity=0;
 for(let frame=0;frame<fps;frame++){
  const next=dampMotion(value,velocity,1,.32,1/fps);
  assert.ok(next.value>=value&&next.value<=1,'Approach a mark without overshooting');
  value=next.value;velocity=next.velocity;
 }
 return {value,velocity};
};
const thirty=run(30),sixty=run(60),oneTwenty=run(120);near(thirty.value,sixty.value);near(sixty.value,oneTwenty.value);near(thirty.velocity,sixty.velocity);
assert.ok(dampMotion(0,0,1,.32,1/60).value<.01,'Important early levels start gently');
const before=dampMotion(0,0,1,.32,.2),reversed=dampMotion(before.value,before.velocity,0,.32,.001);
assert.ok(Math.abs(reversed.velocity-before.velocity)<.1,'Retargeting retains momentum instead of restarting');
near(smoothStep(0),0);near(smoothStep(1),1);near(smoothStep(.5),.5);
assert.ok(smoothStep(.001)<1e-7&&1-smoothStep(.999)<1e-7,'Geometry and camera curves have gentle endpoints');
console.log('Explode intervals double the first two widths, preserve URL values and precise dragging; motion is smooth and independent of frame rate.');
