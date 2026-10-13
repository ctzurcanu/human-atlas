import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/intestinal-lymphoid-detail-geometry';export * from './app/intestine-detail-navigation';export * from './app/shared-tissue-wall';export * from './app/guest-hierarchy';export * from './app/guest-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {intestinalLymphoidDetailGeometry,insideLymphoidEnvelope,intestinalLymphoidDetailUrl,createSharedTissueWall,parseGuestHierarchy,resolveGuestHierarchy,guestNodeChoice}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const sub=(a,b)=>a.map((v,i)=>v-b[i]),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
function boundary(p,indices){const edges=new Map(),adj=new Map();let volume=0;
 for(let i=0;i<indices.length;i+=3){const f=Array.from(indices.slice(i,i+3)),[a,b,c]=f.map(id=>Array.from(p.slice(id*3,id*3+3))),n=cross(sub(b,a),sub(c,a));assert.ok(Math.hypot(...n)>1e-15);volume+=dot(a,n)/6;for(let j=0;j<3;j++){const a=f[j],b=f[(j+1)%3],key=[Math.min(a,b),Math.max(a,b)].join(':'),e=edges.get(key)??[0,0];e[0]++;e[1]+=a<b?1:-1;edges.set(key,e);for(const [x,y] of [[a,b],[b,a]]){if(!adj.has(x))adj.set(x,new Set());adj.get(x).add(y);}}}
 for(const e of edges.values())assert.deepEqual(e,[2,0]);let components=0;const unseen=new Set(adj.keys());while(unseen.size){components++;const stack=[unseen.values().next().value];while(stack.length){const a=stack.pop();if(!unseen.delete(a))continue;for(const b of adj.get(a))if(unseen.has(b))stack.push(b);}}
 return {volume,components,euler:adj.size-edges.size+indices.length/3};
}
test('secondary follicle has convex connected solid core and disjoint finite mantle with valid continuous prisms',()=>{
 const d=intestinalLymphoidDetailGeometry(),wall=createSharedTissueWall({interfaces:[d.outer,d.core],faces:d.faces,rims:new Uint32Array(),layerIds:['mantle']});
 const outer=boundary(d.outer,d.faces),core=boundary(d.core,d.faces),mantle=boundary(wall.position.array,wall.full[0].index.array);
 assert.equal(outer.components,1);assert.equal(core.components,1);assert.equal(mantle.components,2);assert.equal(outer.euler,2);assert.equal(core.euler,2);assert.equal(mantle.euler,4);assert.ok(core.volume>0);assert.ok(mantle.volume>0);assert.ok(Math.abs(outer.volume-core.volume-mantle.volume)<1e-22);
 // Every vertex is within every oriented face half-space: the envelope is
 // convex, which rules out self-crossings and hidden tunnels independently.
 for(const bank of [d.outer,d.core])for(let i=0;i<bank.length;i+=3)assert.ok(insideLymphoidEnvelope(Array.from(bank.slice(i,i+3)),bank,d.faces));
 for(let i=0;i<d.faces.length;i+=3){const f=Array.from(d.faces.slice(i,i+3)),a=f.map(id=>Array.from(d.outer.slice(id*3,id*3+3))),b=f.map(id=>Array.from(d.core.slice(id*3,id*3+3))),gap=a.map((v,j)=>sub(b[j],v));for(const corner of gap){const det=t=>dot(cross(sub(a[1].map((v,j)=>v+t*gap[1][j]),a[0].map((v,j)=>v+t*gap[0][j])),sub(a[2].map((v,j)=>v+t*gap[2][j]),a[0].map((v,j)=>v+t*gap[0][j]))),corner),c0=det(0),end=det(1),half=det(.5),c2=2*(end+c0-2*half),c1=end-c0-c2,t=c2===0?0:-c1/(2*c2);assert.ok(Math.max(c0,end,t>0&&t<1?det(t):-Infinity)<0);}}
 assert.equal(wall.bytes.position,10032);assert.equal(wall.bytes.normal,5016);
 assert.ok(wall.bytes.position+wall.bytes.normal+wall.bytes.closedIndices+d.faces.length*2<31000);
});
test('cut picking distinguishes solid lymphoid center from its mantle and support-plane exterior',()=>{
 const d=intestinalLymphoidDetailGeometry(),contains=(p,bank)=>insideLymphoidEnvelope(p,bank,d.faces);
 assert.ok(contains([0,0,0],d.core));assert.ok(contains([.0004,0,0],d.outer));assert.equal(contains([.0004,0,0],d.core),false);assert.equal(contains([.0006,0,0],d.outer),false);
 // Every distinct occurrence has non-overlapping enclosing AABBs. There is
 // no fabricated inter-follicular bridge or shared lumen in this illustration.
 for(let a=0;a<d.centers.length;a++)for(let b=a+1;b<d.centers.length;b++)assert.ok(d.radii.some((r,j)=>Math.abs(d.centers[a][j]-d.centers[b][j])>2*r));
});
test('both TA98 trees link actual nodule occurrences without fabricating main geometry membership',async()=>{
 const hierarchy=parseGuestHierarchy(JSON.parse(await fs.readFile('public/assets/ta98.json')));
 for(const sex of ['female','male']){const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`)),nodes=new Map(),walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);for(const organ of ['duodenum','jejunum','ileum']){const choice=guestNodeChoice('ta98',nodes.get(`SIW:${organ}:solitary-nodules`));assert.deepEqual(choice.elements,[]);assert.equal(intestinalLymphoidDetailUrl(choice),`./intestinal-lymphoid-review.html?mode=solitary&organ=${organ}`);}const choice=guestNodeChoice('ta98',nodes.get('SIW:ileum:aggregated-nodules'));assert.deepEqual(choice.elements,[]);assert.equal(intestinalLymphoidDetailUrl(choice),'./intestinal-lymphoid-review.html?mode=aggregate&organ=ileum');}
 for(const [term,mode] of [['A05.6.01.013','solitary'],['A05.6.01.014','aggregate']])assert.equal(intestinalLymphoidDetailUrl({id:'node',elements:[],name:'nodule',ta98Term:term}),`./intestinal-lymphoid-review.html?mode=${mode}`);
 assert.equal(intestinalLymphoidDetailUrl({id:'node',elements:[],name:'lymph node',ta98Term:'A13.3.00.001'}),null);
});
