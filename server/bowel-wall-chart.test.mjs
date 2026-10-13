import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {gunzipSync} from 'node:zlib';

for(const sex of ['female','male'])test(`${sex}: shared compact bowel chart preserves course and has a continuous angular seam`,async()=>{
 const root=`.local-models/ta98-review/bowel-wall-main-prototype/${sex}`;
 const atlas=JSON.parse(await readFile(root+'/atlas.json','utf8')),bytes=await readFile(root+'/wall.bin');
 assert.deepEqual(gunzipSync(await readFile(root+'/wall.bin.gz')),bytes);
 const raw=Uint8Array.from(bytes).buffer;
 const sourceRoot='.local-models/ta98-review/small-intestine-layer-v1'+(sex==='male'?'-male':'');
 const meta=JSON.parse(await readFile(sourceRoot+'/manifest.json','utf8')),src=Uint8Array.from(await readFile(sourceRoot+'/wall.bin')).buffer;
 for(const organ of meta.organs){
  const parts=atlas.parts.filter(p=>p.bowelWallSource&&p.id.includes(`:BOWELWALL:${organ.id}:`));assert.equal(parts.length,12);
  for(const part of parts){
   const tissue=part.tissueAppearance;
   if(part.conceptId===`SIW:${organ.id}:epithelium`||part.conceptId===`BOWEL-FOLD:${organ.id}:epithelium`){assert.equal(tissue.pattern,'intestinal');assert.equal(tissue.scale,.00028);assert.equal(tissue.relief,.000012);assert.match(part.description,/lighting only/);}
   if(part.conceptId===`SIW:${organ.id}:lamina-propria`||part.conceptId===`BOWEL-FOLD:${organ.id}:lamina-propria`)assert.equal(tissue.pattern,'fibrous');
   assert.equal(atlas.materials[part.material].map,undefined);
  }
  assert.equal(new Set(parts.map(p=>p.tissueChart)).size,1);
  const n=organ.vertices,chart=new Int16Array(raw,parts[0].tissueChart,n*9*3),original=new Float32Array(src,organ.chart,n*2);
  const carrier=new Float32Array(src,organ.outer,n*3),bank=new Float32Array(raw,parts[0].positions,n*9*3);
  assert.deepEqual(bank.subarray(0,n*3),carrier);
  for(const part of parts.filter(p=>!/:epithelium$|:lamina-propria$|:muscularis-mucosae$|:submucosa$|:circular-folds(?:$|:)/.test(p.id)))assert.equal(part.indexCount,organ.triangles*6+organ.rims*6);
  const length=parts[0].tissueAppearance.courseLength;assert.ok(length>0);
  for(let k=0;k<9;k++)for(let i=0;i<n;i++){
   const at=(k*n+i)*3,s=chart[at]/32767*length,c=chart[at+1]/32767,t=chart[at+2]/32767,angle=original[i*2+1]*2*Math.PI;
   assert.ok(Math.abs(s-original[i*2])<=length/65534+1e-12);
   assert.ok(Math.abs(c-Math.cos(angle))<=1/65534+1e-12);
   assert.ok(Math.abs(t-Math.sin(angle))<=1/65534+1e-12);
   assert.ok(Math.abs(Math.hypot(c,t)-1)<.00003);
   if(k)assert.deepEqual(chart.subarray(at,at+3),chart.subarray(i*3,i*3+3));
  }
  // Exported controls need not retain ring order. Inspect actual triangle edges
  // crossing the scalar-angle seam, rather than assuming consecutive rings.
  const faces=new Uint32Array(src,organ.faces,organ.triangles*3);let seams=0;
  for(let j=0;j<faces.length;j+=3)for(let d=0;d<3;d++){
   const a=faces[j+d],b=faces[j+(d+1)%3],delta=Math.abs(original[a*2+1]-original[b*2+1]);
   if(delta<=.5)continue;
   const chord=Math.hypot((chart[a*3+1]-chart[b*3+1])/32767,(chart[a*3+2]-chart[b*3+2])/32767);
   assert.ok(Math.abs(chord-2*Math.abs(Math.sin(delta*Math.PI)))<.00005);seams++;
  }
  assert.ok(seams>0);
 }
});
