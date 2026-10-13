import assert from 'node:assert/strict';
import test from 'node:test';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export {createSharedTissueWall} from './app/shared-tissue-wall'; export {createSharedTissueCaps} from './app/shared-tissue-caps'; export {Scene,Matrix4,Plane,Vector3} from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {createSharedTissueWall,createSharedTissueCaps,Scene,Matrix4,Plane,Vector3}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
function fixture(){
 const faces=[],rims=[],n=12;
 for(let j=0;j<n;j++){const k=(j+1)%n;faces.push(j,k,n+j,k,n+k,n+j);rims.push(j,k,n+k,n+j);}
 return {interfaces:[1,.9,.7].map(r=>new Float32Array(Array.from({length:n*2},(_,i)=>[r*Math.cos(i%n*2*Math.PI/n),r*Math.sin(i%n*2*Math.PI/n),i<n?0:2]).flat())),faces:new Uint32Array(faces),rims:new Uint32Array(rims),layerIds:['outer','inner']};
}
function closed(g){
 const edges=new Map();for(let i=0;i<g.index.count;i+=3){const t=[0,1,2].map(j=>g.index.getX(i+j));for(let j=0;j<3;j++){const a=t[j],b=t[(j+1)%3],key=[a,b].sort((a,b)=>a-b).join(',');const row=edges.get(key)??[0,0];row[0]++;row[1]+=a<b?1:-1;edges.set(key,row);}}
 return [...edges.values()].every(([count,winding])=>count===2&&winding===0);
}
test('Peelable shells share immutable buffers, remain closed, and leave both tube mouths open',()=>{
 const input=fixture(),wall=createSharedTissueWall(input),bank=input.interfaces.length*24;
 for(const g of [...wall.full,...wall.displayed]){assert.equal(g.attributes.position,wall.position);assert.equal(g.attributes.normal,wall.normal);assert.ok(closed(g));}
 for(const g of wall.full){const start=input.faces.length*2;for(let i=start;i<g.index.count;i+=3){const zs=[0,1,2].map(k=>wall.position.getZ(g.index.getX(i+k)));assert.ok(zs.every(z=>z===zs[0]),'annular rim stays at an actual mouth');const p=[0,1,2].map(k=>g.index.getX(i+k)%24);assert.ok(new Set(p.map(x=>x%12)).size<=2,'no triangle spans a lumen disk');}}
 assert.equal(wall.position.count,bank);assert.equal(wall.normal.count,bank);
 for(let k=0;k<wall.full.length;k++){assert.equal(wall.full[k].index.getX(input.faces.length),(k+1)*24+input.faces[0]);assert.deepEqual(wall.full[k].groups.map(g=>g.materialIndex),[0,1,2]);}
 wall.dispose();
});
test('Visibility and peel update only draw groups, without allocating coordinates or indices',()=>{
 const wall=createSharedTissueWall(fixture()),indices=wall.displayed.map(g=>g.index),positions=wall.position.array;
 wall.update([true,true]);assert.equal(wall.displayed[0].groups.length,2);assert.equal(wall.displayed[1].groups.length,2);
 wall.update([false,true]);assert.equal(wall.displayed[0].drawRange.count,0);assert.equal(wall.displayed[1].groups.length,3);
 wall.update([true,true],[true]);assert.equal(wall.displayed[0].groups.length,3);assert.equal(wall.displayed[1].groups.length,3);
 assert.equal(wall.update([true,true],[true]),false);
 assert.equal(wall.position.array,positions);wall.displayed.forEach((g,k)=>assert.equal(g.index,indices[k]));wall.dispose();
});
test('Malformed fields fail before allocating tissue geometry',()=>{
 const f=fixture();assert.throws(()=>createSharedTissueWall({...f,interfaces:f.interfaces.slice(1)}));assert.throws(()=>createSharedTissueWall({...f,layerIds:['duplicate','duplicate']}));assert.throws(()=>createSharedTissueWall({...f,faces:new Uint32Array([0,1,1000])}));
});

test('Moving cuts follow the plane immediately and keep the closed counting buffers intact',()=>{
 const wall=createSharedTissueWall(fixture()),scene=new Scene(),plane=new Plane(new Vector3(0,0,-1),0);
 const caps=createSharedTissueCaps(scene,wall.full,[0,1].map(()=>({color:[.8,.4,.4],pattern:'circular',scale:.001,relief:.000001})),plane);
 const matrices=[new Matrix4(),new Matrix4().makeTranslation(.2,0,0)];
 const originals=caps.counters.map(m=>m.geometry.getAttribute('position'));
 for(let step=0;step<30;step++){plane.constant=step/20;caps.update([true,true],matrices,true);caps.caps.forEach(m=>{assert.ok(Math.abs(plane.distanceToPoint(m.position))<1e-12);assert.equal(m.visible,true);});caps.counters.forEach((m,i)=>assert.equal(m.geometry.getAttribute('position'),originals[i]));}
 caps.update([false,true],matrices,true);assert.equal(caps.caps[0].visible,false);assert.equal(caps.caps[1].visible,true);
 caps.update([true,true],matrices,false);assert.ok([...caps.caps,...caps.counters].every(m=>!m.visible));
 caps.dispose();assert.equal(scene.children.length,0);wall.dispose();
});
