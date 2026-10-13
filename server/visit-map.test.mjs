import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
const result=await build({entryPoints:['app/visit-map.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {visitMapSegments,visitRegions}=await import('data:text/javascript;base64,'+Buffer.from(result.outputFiles[0].contents).toString('base64'));
// Duck-typed bounds/curve keep this independent of bundled Three identities.
const bounds={getSize:v=>v.set(.6,1.8,.4)};
test('Map segments meet at the exact slider boundary without color interpolation across organs',()=>{
 const regions=visitRegions(bounds,'front');
 const curve={getPointAt:t=>({toArray:()=>[t,1-t,t*t],x:t,y:1-t,z:t*t})};
 const map=visitMapSegments(curve,regions);
 for(let i=0;i<map.positions.length;i+=6){
  const from=map.positions[i],to=map.positions[i+3];
  assert.ok(!regions.some(r=>r.start>from+1e-12&&r.start<to-1e-12));
  assert.deepEqual(map.colors.slice(i,i+3),map.colors.slice(i+3,i+6));
 }
 assert.deepEqual(map.positions.slice(0,3),[0,1,0]);
 assert.deepEqual(map.positions.slice(-3),[1,0,1]);
});
test('Orbit map uses the same four 90-degree colors and head-to-feet uses measured height',()=>{
 const orbit=visitRegions(bounds,'orbit');assert.deepEqual(orbit.map(r=>r.start),[0,.25,.5,.75]);
 assert.deepEqual(orbit.map(r=>r.name),['0–90°','90–180°','180–270°','270–360°']);
 const linear=visitRegions(bounds,'front');assert.equal(linear[0].name,'Head · 0–25 cm');assert.equal(linear.at(-1).name,'Feet · 169–180 cm');
 assert.deepEqual(orbit.map(r=>r.color),linear.slice(0,4).map(r=>r.color));
});
