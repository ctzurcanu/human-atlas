import test from 'node:test';import assert from 'node:assert/strict';import {build} from 'esbuild';
const compiled=await build({entryPoints:['app/duodenal-gland-detail-geometry.ts'],bundle:true,write:false,platform:'node',format:'esm'});
const {duodenalGlandDetailGeometry,duodenalGlandField}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const d=duodenalGlandDetailGeometry();const sub=(a,b)=>a.map((v,i)=>v-b[i]),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0),p=i=>Array.from(d.positions.slice(i*3,i*3+3));
test('complete branched gland tissue is finite, closed, oriented and connected',()=>{
 const edges=new Map(),neighbors=new Map();let volume=0;
 for(let i=0;i<d.faces.length;i+=3){const f=Array.from(d.faces.slice(i,i+3)),[a,b,c]=f.map(p),n=cross(sub(b,a),sub(c,a));assert.ok(Math.hypot(...n)>1e-18,`degenerate face ${i/3}`);volume+=dot(a,n)/6;
  for(let j=0;j<3;j++){const a=f[j],b=f[(j+1)%3],key=[Math.min(a,b),Math.max(a,b)].join(':'),e=edges.get(key)??[0,0];e[0]++;e[1]+=a<b?1:-1;edges.set(key,e);if(!neighbors.has(a))neighbors.set(a,new Set());neighbors.get(a).add(b);neighbors.get(a).add(f[(j+2)%3]);}
 }
 for(const e of edges.values())assert.deepEqual(e,[2,0]);assert.ok(volume>0);const visited=new Set(),queue=[0];while(queue.length){const v=queue.pop();if(visited.has(v))continue;visited.add(v);queue.push(...neighbors.get(v));}assert.equal(visited.size,d.positions.length/3);assert.equal(d.positions.length/3-edges.size+d.faces.length/3,2);
 assert.ok(d.positions.byteLength+d.faces.byteLength<1_000_000);assert.equal(d.branches,7);
});
test('continuous duct and branch paths are empty, gland wall and specimen mouth preserve lumen',()=>{
 for(let i=0;i<80;i++){const y=d.mouthMeters-(d.mouthMeters+.00012)*i/79;assert.ok(duodenalGlandField([0,y,0])>0);}
 for(const side of [-1,1])for(let i=0;i<=40;i++){const t=i/40;assert.ok(duodenalGlandField([side*.00007*t,-.00012-.00008*t,side*.000016*t])>0);}
 assert.ok(duodenalGlandField([.00002,-.00005,0])<0);assert.ok(duodenalGlandField([0,.00002,0])>0);
 const mouth=[];for(let i=0;i<d.positions.length/3;i++)if(Math.abs(d.positions[i*3+1]-d.mouthMeters)<1e-11)mouth.push(p(i));assert.ok(mouth.length>10);for(const v of mouth){const r=Math.hypot(v[0],v[2]);assert.ok(r>.000015);assert.ok(r<.000035);}
});
test('both actual TA98 hierarchies retain missing duodenal glands and reject unsupported detail mapping',async()=>{
 const fs=await import('node:fs/promises'),r=await build({stdin:{contents:"export * from './app/intestine-detail-navigation';export * from './app/guest-hierarchy';export * from './app/guest-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false}),m=await import('data:text/javascript;base64,'+Buffer.from(r.outputFiles[0].contents).toString('base64'));
 const hierarchy=m.parseGuestHierarchy(JSON.parse(await fs.readFile('public/assets/ta98.json')));
 for(const sex of ['female','male']){const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`)),nodes=new Map(),walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};m.resolveGuestHierarchy(atlas,hierarchy).forEach(walk);const gland=nodes.get('SIW:duodenum:duodenal-glands'),parent=nodes.get('SIW:duodenum:submucosa');assert.ok(parent.children.includes(gland));const choice=m.guestNodeChoice('ta98',gland);assert.deepEqual(choice.elements,[]);assert.equal('duodenalGlandDetailUrl' in m,false);assert.equal(m.intestinalCryptDetailUrl(choice),null);}
 const inspector=await fs.readFile('app/selection-inspector.tsx','utf8'),viewer=await fs.readFile('web/small-intestine-layer-review.ts','utf8');assert.ok(!inspector.includes('duodenalGlandDetail'));assert.ok(!viewer.includes('/duodenal-gland-review.html'));
});
test('actual packed boundary keeps sampled duct paths empty and tissue finite',()=>{
 function winding(point){let total=0;for(let i=0;i<d.faces.length;i+=3){const [a,b,c]=Array.from(d.faces.slice(i,i+3)).map(id=>sub(p(id),point)),A=Math.hypot(...a),B=Math.hypot(...b),C=Math.hypot(...c);total+=2*Math.atan2(dot(a,cross(b,c)),A*B*C+dot(a,b)*C+dot(b,c)*A+dot(c,a)*B);}return total/(4*Math.PI);}
 for(const point of [[0,0,0],[0,-.00005,0],[0,-.00012,0],[-.00007,-.00020,-.000016],[.00007,-.00020,.000016],[-.00014,-.00028,-.000065],[.00014,-.00028,.000065]])assert.ok(Math.abs(winding(point))<1e-6,`packed boundary obstructs lumen at ${point}`);
 assert.ok(Math.abs(winding([.000025,-.00005,0])-1)<1e-6);assert.ok(Math.abs(winding([0,.00002,0]))<1e-6);
});
