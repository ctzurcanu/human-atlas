import fs from 'node:fs/promises';import {build} from 'esbuild';import {MeshoptSimplifier} from 'meshoptimizer';
const r=await build({entryPoints:['app/duodenal-gland-detail-geometry.ts'],bundle:true,write:false,format:'esm',platform:'node'}),m=await import('data:text/javascript;base64,'+Buffer.from(r.outputFiles[0].contents).toString('base64')),source=m.duodenalGlandDetailGeometry();await MeshoptSimplifier.ready;
const folder='.local-models/ta98-review/duodenal-gland-detail';await fs.mkdir(folder,{recursive:true});
const locks=new Uint8Array(source.positions.length/3);for(let i=0;i<locks.length;i++)locks[i]=Number(source.positions[i*3+1]>=source.mouthMeters-.000012);
await fs.writeFile(folder+'/source.json',JSON.stringify({positions:Array.from(source.positions),faces:Array.from(source.faces),mouthMeters:source.mouthMeters}));
for(const microns of [.5,.75,1,2]){
 const [indices,error]=MeshoptSimplifier.simplifyWithAttributes(new Uint32Array(source.faces),source.positions,3,new Float32Array(),0,[],locks,3000*3,microns*.000001,['ErrorAbsolute']);
 const oldIds=Array.from(new Set(indices)).sort((a,b)=>a-b),remap=new Map(oldIds.map((id,i)=>[id,i])),positions=new Float32Array(oldIds.length*3);oldIds.forEach((id,i)=>positions.set(source.positions.subarray(id*3,id*3+3),i*3));const faces=new Uint16Array(Array.from(indices,id=>remap.get(id)));
 await fs.writeFile(folder+`/candidate-${microns}.json`,JSON.stringify({positions:Array.from(positions),faces:Array.from(faces),sourceIds:oldIds,estimatedErrorMeters:error,errorLimitMeters:microns*.000001,mouthMeters:source.mouthMeters}));console.log(microns+' µm:',oldIds.length,'controls,',indices.length/3,'triangles; estimated error',error);
}
