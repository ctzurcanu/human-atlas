import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/villus-detail-geometry';export * from './app/intestine-detail-navigation';export * from './app/guest-hierarchy';export * from './app/guest-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {villusDetailGeometry,insideVillusCarrier,intestinalVillusDetailUrl,parseGuestHierarchy,resolveGuestHierarchy,guestNodeChoice}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
function closed(vertices,faces){const edges=new Map();let volume=0;for(let i=0;i<faces.length;i+=3){const ids=faces.slice(i,i+3),[a,b,c]=ids.map(id=>vertices.slice(id*3,id*3+3));const u=b.map((v,j)=>v-a[j]),v=c.map((w,j)=>w-a[j]),n=[u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]];assert.ok(Math.hypot(...n)>1e-15);volume+=a.reduce((s,x,j)=>s+x*n[j],0)/6;for(let j=0;j<3;j++){const a=ids[j],b=ids[(j+1)%3],key=[Math.min(a,b),Math.max(a,b)].join(':'),prior=edges.get(key)??[0,0];prior[0]++;prior[1]+=a<b?1:-1;edges.set(key,prior);}}
 for(const [count,winding] of edges.values()){assert.equal(count,2);assert.equal(winding,0);}assert.ok(volume>0);return volume;
}
test('villus fragment has two closed, positively oriented finite cells and shared epithelium/core interface',()=>{
 const d=villusDetailGeometry(),n=d.count,bank=[...d.outer,...d.inner],shell=[...d.faces];for(let i=0;i<d.faces.length;i+=3)shell.push(n+d.faces[i],n+d.faces[i+2],n+d.faces[i+1]);for(let i=0;i<d.rims.length;i+=2){const a=d.rims[i],b=d.rims[i+1];shell.push(b,a,n+a,b,n+a,n+b);}
 const core=[...d.inner,0,0,0],v1=closed(bank,shell),v2=closed(core,Array.from(d.coreFaces));assert.ok(v1<v2);assert.equal(shell.length/3+d.coreFaces.length/3,864);assert.equal(d.heightMeters,.000750);
 for(let i=0;i<n;i++){const length=Math.hypot(...[0,1,2].map(j=>d.outer[i*3+j]-d.inner[i*3+j]));assert.ok(Math.abs(length-.000025)<1e-10);}assert.ok(d.outer.every(Number.isFinite));
});
test('cap picking distinguishes core, epithelial annulus and empty support plane',()=>{
 const d=villusDetailGeometry();for(const [point,outer,inner] of [[[0,.0004,0],true,true],[[.00008,.0004,0],true,false],[[.0002,.0004,0],false,false],[[0,.0008,0],false,false],[[0,-.0001,0],false,false]]){assert.equal(insideVillusCarrier(point,d.outer,d.sides),outer);assert.equal(insideVillusCarrier(point,d.inner,d.sides),inner);}
});
test('both actual TA98 trees link villi to an independent detail without mapping missing main geometry',async()=>{
 const hierarchy=parseGuestHierarchy(JSON.parse(await fs.readFile('public/assets/ta98.json')));
 for(const sex of ['female','male']){const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`));const nodes=new Map();const walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);for(const organ of ['duodenum','jejunum','ileum']){const choice=guestNodeChoice('ta98',nodes.get(`SIW:${organ}:villi`));assert.deepEqual(choice.elements,[]);assert.equal(intestinalVillusDetailUrl(choice),`./intestinal-villus-review.html?organ=${organ}`);assert.equal(intestinalVillusDetailUrl(guestNodeChoice('ta98',nodes.get(`SIW:${organ}:circular-folds`))),null);}}
 assert.equal(intestinalVillusDetailUrl({id:'selection-set',elements:[],name:'villi'}),null);
});
