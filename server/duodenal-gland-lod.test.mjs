import test from 'node:test';import assert from 'node:assert/strict';import fs from 'node:fs/promises';
const source=JSON.parse(await fs.readFile('.local-models/ta98-review/duodenal-gland-detail/source.json')),candidate=JSON.parse(await fs.readFile('.local-models/ta98-review/duodenal-gland-detail/candidate-0.75.json'));
const d={positions:new Float32Array(candidate.positions),faces:new Uint16Array(candidate.faces)},p=i=>Array.from(d.positions.slice(i*3,i*3+3)),sub=(a,b)=>a.map((v,i)=>v-b[i]),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a.reduce((s,v,i)=>s+v*b[i],0);
test('lighter gland retains original coordinates, exact outlet triangles and fewer resources',()=>{
 for(let i=0;i<candidate.sourceIds.length;i++)assert.deepEqual(p(i),source.positions.slice(candidate.sourceIds[i]*3,candidate.sourceIds[i]*3+3));
 const mouthFaces=(positions,faces,remap)=>{const out=[];for(let i=0;i<faces.length;i+=3){const f=faces.slice(i,i+3);if(f.every(id=>Math.abs(positions[id*3+1]-source.mouthMeters)<1e-11))out.push(f.map(id=>remap?remap[id]:id).join(':'));}return out.sort();};
 assert.deepEqual(mouthFaces(candidate.positions,candidate.faces,candidate.sourceIds),mouthFaces(source.positions,source.faces));assert.equal(d.faces.length/3,6088);assert.ok(d.positions.byteLength*2+d.faces.byteLength<120_000);
});
test('lighter actual boundary preserves sampled duct paths and finite tissue',()=>{
 function winding(point){let total=0;for(let i=0;i<d.faces.length;i+=3){const [a,b,c]=Array.from(d.faces.slice(i,i+3)).map(id=>sub(p(id),point)),A=Math.hypot(...a),B=Math.hypot(...b),C=Math.hypot(...c);total+=2*Math.atan2(dot(a,cross(b,c)),A*B*C+dot(a,b)*C+dot(b,c)*A+dot(c,a)*B);}return total/(4*Math.PI);}
 for(let i=0;i<=16;i++)assert.ok(Math.abs(winding([0,-.00012*i/16,0]))<1e-6);
 for(const side of [-1,1])for(const z of [-1,1])for(let i=0;i<=8;i++){const t=i/8;assert.ok(Math.abs(winding([side*(.00007+.00007*t),-.00020-.00008*t,side*(.000016+(.000065*z-.000016)*t)]))<1e-6);}
 assert.ok(Math.abs(winding([.000025,-.00005,0])-1)<1e-6);assert.ok(Math.abs(winding([0,.00002,0]))<1e-6);
});
test('section picking uses the retained boundary, excluding empty support quads and lumen',async()=>{
 const {build}=await import('esbuild'),r=await build({stdin:{contents:"export * from './app/closed-tissue-picking';export * from 'three';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false}),m=await import('data:text/javascript;base64,'+Buffer.from(r.outputFiles[0].contents).toString('base64'));
 const g=new m.BufferGeometry();g.setAttribute('position',new m.BufferAttribute(d.positions,3));g.setIndex(new m.BufferAttribute(d.faces,1));const mesh=new m.Mesh(g,new m.MeshBasicMaterial({side:m.DoubleSide}));mesh.updateMatrixWorld();
 for(const point of [[0,0,0],[0,-.00005,0],[0,-.00012,0],[0,.00002,0],[.0002,-.00005,0]])assert.equal(m.insideClosedTissue(new m.Vector3(...point),mesh),false);
 assert.equal(m.insideClosedTissue(new m.Vector3(.000025,-.00005,0),mesh),true);
 g.dispose();mesh.material.dispose();
});
