import assert from 'node:assert/strict';
import {build} from 'esbuild';

const built=await build({entryPoints:['app/xr-orbit.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {advanceControllerOrbit,orbitAngleDelta}=await import(`data:text/javascript;base64,${Buffer.from(built.outputFiles[0].text).toString('base64')}`);
const near=(actual,expected)=>assert.ok(Math.abs(actual-expected)<1e-9,`${actual} ≠ ${expected}`);
near(orbitAngleDelta(-Math.PI+.1,Math.PI-.1),.2);
assert.deepEqual(advanceControllerOrbit(null,{yaw:1,pitch:1},{yaw:0,pitch:0}),{yaw:0,pitch:0},'pressing the trigger establishes a baseline without jumping');
const step=advanceControllerOrbit({yaw:.2,pitch:.1},{yaw:.3,pitch:.15},{yaw:0,pitch:0});
near(step.yaw,-.1);near(step.pitch,-.05);
assert.deepEqual(advanceControllerOrbit({yaw:0,pitch:0},{yaw:1,pitch:1},step),step,'tracking jumps do not move the body');
let orbit={yaw:0,pitch:0};for(let i=0;i<30;i++)orbit=advanceControllerOrbit({yaw:0,pitch:0},{yaw:0,pitch:.1},orbit);
near(orbit.pitch,-1.25);
console.log('Quest controller orbit stays continuous across angle wrap, ignores jumps, and limits pitch.');
