import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';import {build} from 'esbuild';
const bundle=await build({stdin:{contents:"export * from './app/crypt-detail-geometry';export * from './app/villus-detail-geometry';export * from './app/intestine-detail-navigation';export * from './app/shared-tissue-wall';export * from './app/guest-hierarchy';export * from './app/guest-choice';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {cryptDetailGeometry,insideVillusCarrier,intestinalCryptDetailUrl,createSharedTissueWall,parseGuestHierarchy,resolveGuestHierarchy,guestNodeChoice}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const sub=(a,b)=>a.map((v,i)=>v-b[i]),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
test('crypt cells have closed oriented boundaries, finite positive volume and no inverted continuous prisms',()=>{
 const d=cryptDetailGeometry(),banks=[d.stroma,d.epithelium,d.lumen],wall=createSharedTissueWall({interfaces:banks,faces:d.faces,rims:d.rims,layerIds:['stroma','epithelium']});let triangles=0;
 for(let k=0;k<2;k++){
  const g=wall.full[k],p=g.attributes.position.array,ids=g.index.array,edges=new Map();let volume=0;
  for(let i=0;i<ids.length;i+=3){const f=Array.from(ids.slice(i,i+3)),[a,b,c]=f.map(id=>Array.from(p.slice(id*3,id*3+3))),n=cross(sub(b,a),sub(c,a));assert.ok(Math.hypot(...n)>1e-15);volume+=dot(a,n)/6;for(let j=0;j<3;j++){const a=f[j],b=f[(j+1)%3],key=[Math.min(a,b),Math.max(a,b)].join(':'),e=edges.get(key)??[0,0];e[0]++;e[1]+=a<b?1:-1;edges.set(key,e);}}
  for(const e of edges.values())assert.deepEqual(e,[2,0]);assert.ok(volume>0);triangles+=ids.length/3;
  for(let i=0;i<d.faces.length;i+=3){const f=Array.from(d.faces.slice(i,i+3)),a=f.map(id=>Array.from(banks[k].slice(id*3,id*3+3))),b=f.map(id=>Array.from(banks[k+1].slice(id*3,id*3+3))),gap=a.map((v,j)=>sub(b[j],v));
   for(const corner of gap){const determinant=t=>dot(cross(sub(a[1].map((v,j)=>v+t*gap[1][j]),a[0].map((v,j)=>v+t*gap[0][j])),sub(a[2].map((v,j)=>v+t*gap[2][j]),a[0].map((v,j)=>v+t*gap[0][j]))),corner),c0=determinant(0),end=determinant(1),half=determinant(.5),c2=2*(end+c0-2*half),c1=end-c0-c2,t=c2===0?0:-c1/(2*c2);let max=Math.max(c0,end);if(t>0&&t<1)max=Math.max(max,determinant(t));assert.ok(max<0,`inverted cell ${k}, triangle ${i/3}`);}
  }
 }
 assert.equal(triangles,896);
 for(let i=0;i<d.count;i++){const distance=(a,b)=>Math.hypot(...[0,1,2].map(j=>a[i*3+j]-b[i*3+j]));assert.ok(Math.abs(distance(d.epithelium,d.lumen)-.000015)<1e-10);assert.ok(Math.abs(distance(d.stroma,d.epithelium)-.00004)<1e-10);}
});
test('crypt mouth and section lumen stay empty, rounded blind base belongs to epithelial tissue',()=>{
 const d=cryptDetailGeometry(),contains=(p,b)=>insideVillusCarrier(p,b,d.sides);
 assert.equal(contains([0,0,0],d.lumen),true);assert.equal(contains([0,-.0001,0],d.lumen),true);
 assert.equal(contains([0,-.000247,0],d.lumen),false);assert.equal(contains([0,-.000247,0],d.epithelium),true);
 assert.equal(contains([.00007,-.0001,0],d.epithelium),false);assert.equal(contains([.00007,-.0001,0],d.stroma),true);
 assert.equal(contains([.0002,-.0001,0],d.stroma),false);assert.equal(contains([0,.00001,0],d.stroma),false);
 // Every mouth-rim bridge is annular, with no axial center vertex or disk.
 assert.equal(d.rims.length,32);for(const id of d.rims){assert.equal(Math.abs(d.lumen[id*3+1]),0);assert.ok(Math.hypot(d.lumen[id*3],d.lumen[id*3+2])>0);}
});
test('male/female gland hierarchy links preserve missing main geometry and exclude duodenal submucosal glands',async()=>{
 const hierarchy=parseGuestHierarchy(JSON.parse(await fs.readFile('public/assets/ta98.json')));
 for(const sex of ['female','male']){const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`)),nodes=new Map(),walk=n=>{nodes.set(n.id,n);n.children.forEach(walk);};resolveGuestHierarchy(atlas,hierarchy).forEach(walk);for(const organ of ['duodenum','jejunum','ileum']){const choice=guestNodeChoice('ta98',nodes.get(`SIW:${organ}:intestinal-glands`));assert.deepEqual(choice.elements,[]);assert.equal(intestinalCryptDetailUrl(choice),`./intestinal-crypt-review.html?organ=${organ}`);}}
 assert.equal(intestinalCryptDetailUrl({id:'anything',elements:[],name:'glands',ta98Term:'A05.6.02.017'}),null);
 assert.equal(intestinalCryptDetailUrl({id:'anything',elements:[],name:'glands',ta98Term:'A05.6.01.012'}),'./intestinal-crypt-review.html');
});
