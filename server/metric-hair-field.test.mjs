import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/metric-hair-field';",resolveDir:process.cwd(),loader:'ts'},bundle:true,platform:'node',format:'esm',write:false});
const {metricHairGeometry,MetricHairPicking}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].contents).toString('base64')}`);
const report=JSON.parse(readFileSync('reports/ta98-female/metric-hair-field.json'));const field=report.field;const bytes=readFileSync('.local-models/'+field.url.replace('/local-models/',''));const raw=bytes.buffer.slice(bytes.byteOffset,bytes.byteOffset+bytes.byteLength);
test('shared closed hair shaft has requested physical length and diameter with root-to-tip UVs',()=>{
 const g=metricHairGeometry(field,raw,2);assert.equal(g.instanceCount,field.count);assert.equal(field.length,.001);assert.equal(field.diameter,.0001);assert.ok(Math.abs(report.meanNearestNeighborMm-.6)<.01);
 const p=g.getAttribute('position'),uv=g.getAttribute('uv');let min=Infinity,max=-Infinity,radius=0;
 for(let i=0;i<p.count;i++){min=Math.min(min,p.getY(i));max=Math.max(max,p.getY(i));radius=Math.max(radius,Math.hypot(p.getX(i),p.getZ(i)));assert.ok(uv.getX(i)>=0&&uv.getX(i)<=1);assert.ok(uv.getY(i)>=0&&uv.getY(i)<=1);}
 assert.ok(Math.abs(min)<1e-10);assert.ok(Math.abs(max-.001)<1e-9);assert.ok(Math.abs(radius-.00005)<1e-10);assert.equal(g.index.count/3,24);
 const edges=new Map(),positions=new Map();let next=0;const ids=[];
 for(let i=0;i<p.count;i++){const key=[p.getX(i),p.getY(i),p.getZ(i)].map(v=>Math.round(v*1e10)).join(',');if(!positions.has(key))positions.set(key,next++);ids.push(positions.get(key));}
 for(let i=0;i<g.index.count;i+=3){const t=[0,1,2].map(k=>ids[g.index.getX(i+k)]);for(let k=0;k<3;k++){const key=[t[k],t[(k+1)%3]].sort((a,b)=>a-b).join(':');edges.set(key,(edges.get(key)||0)+1);}}
 assert.ok([...edges.values()].every(n=>n===2));g.dispose();
});

test('male and female share the requested hair dimensions and have continuous hairlines',()=>{
 for(const sex of ['male','female']){const r=JSON.parse(readFileSync(`reports/ta98-${sex}/metric-hair-field.json`));assert.equal(r.continuousHairline,true);assert.equal(r.field.length,.001);assert.equal(r.field.diameter,.0001);assert.ok(Math.abs(r.meanNearestNeighborMm-.6)<.025);}
 const male=JSON.parse(readFileSync('.local-models/ta98-runtime.json'));assert.ok(!male.parts.some(p=>p.id==='ZA:Hairs of head'));assert.equal(male.parts.filter(p=>p.hairField).length,1);assert.equal(male.parts.find(p=>p.hairField).id,'ZA:Scalp hair:metric');
 const iris=male.parts.find(p=>p.id==='ZA:Iris.l');assert.deepEqual(male.materials[iris.material].color,[54,32,23]);const cornea=male.parts.find(p=>p.name==='Cornea (left)');assert.ok(cornea);assert.ok(male.materials[cornea.material].opacity<.2);
 for(const side of ['l','r']){const sclera=male.parts.find(p=>p.id===`ZA:Sclera.${side}`);assert.ok(sclera,'actual native sclera must be present');const material=male.materials[sclera.material];assert.deepEqual(material.color,[255,255,255],'both actual sclera meshes must be white');assert.equal(material.map,undefined,'illustrative grain must not tint eye whites');}
});
