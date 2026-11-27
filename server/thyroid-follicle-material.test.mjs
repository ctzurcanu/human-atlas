import test from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
import {Matrix4,ShaderLib} from 'three';

const compiled=await build({stdin:{contents:"export * from './app/tissue-wall-material';",resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false});
const {tissueWallMaterial,applyPairedNormalSigns}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
const appearance=(extra={})=>({color:[.64,.28,.21],pattern:'thyroid-follicles',scale:.0005,relief:.000003,...extra});
function shaderFor(a=appearance(),frame,inner=false){
 const result=tissueWallMaterial(a,frame,inner),shader={uniforms:{},vertexShader:ShaderLib.standard.vertexShader,fragmentShader:ShaderLib.standard.fragmentShader};
 result.material.onBeforeCompile(shader,{});return {...result,shader};
}
const fract=x=>x-Math.floor(x),mix=(a,b,t)=>a+(b-a)*t;
function smoothstep(a,b,x){const t=Math.max(0,Math.min(1,(x-a)/(b-a)));return t*t*(3-2*t);}
const hash=id=>{
 const h=id.map((v,k)=>fract(v*[.1031,.1030,.0973][k]));
 const dot=h[0]*(h[1]+33.33)+h[1]*(h[0]+33.33)+h[2]*(h[2]+33.33),p=h.map(v=>v+dot);
 return [fract((p[0]+p[1])*p[2]),fract((p[0]+p[2])*p[1]),fract((p[1]+p[2])*p[0])];
};
function packing(shader){
 const source=shader.fragmentShader,constant=name=>Number(source.match(new RegExp(`const float ${name}=([^;]+);`))[1]);
 const minimum=constant('minimum'),maximum=constant('maximum'),shell=constant('shell'),padding=constant('padding');
 const frame=source.match(/mat3 packingFrame=mat3\(([\s\S]*?)\);/)[1].match(/vec3\(([^)]+)\)/g).map(s=>s.slice(5,-1).split(',').map(Number));
 const forward=p=>[0,1,2].map(r=>frame.reduce((sum,c,i)=>sum+c[r]*p[i],0)),reverse=p=>frame.map(c=>c.reduce((sum,v,i)=>sum+v*p[i],0));
 const stepLength=maximum/Math.SQRT2;
 const site=id=>{
  const large=((id[0]+id[1]+id[2])%2+2)%2===0,container=large?maximum*.5:stepLength-maximum*.5-padding;
  const maxDiameter=2*(container-padding);if(maxDiameter<minimum)return null;
  const random=hash(id),sizeRandom=hash(id.map((v,k)=>v+[19.13,3.73,11.17][k]));
  const radius=mix(minimum,maxDiameter,sizeRandom[0]**.3)*.5,direction=random.map(v=>v*2-1),norm=Math.max(Math.hypot(...direction),.0001);
  const travel=Math.max(0,container-radius-padding)*sizeRandom[1]**(1/3),centre=id.map((v,k)=>v*stepLength+direction[k]/norm*travel);
  return {id,large,container,radius,travel,centre};
 };
 return {minimum,maximum,shell,padding,stepLength,forward,reverse,site};
}

// Numerical oracle uses the physical spherical packing contract. Read emitted
// shader constants so a metre/scale mismatch cannot be concealed by test inputs.
function fieldOracle(shader,point,footprint=0){
 const source=shader.fragmentShader,{maximum,shell,stepLength,forward,site}=packing(shader);
 const average=source.match(/const vec3 average=vec3\(([^)]+)\)/)[1].split(',').map(Number);
 const detail=1-smoothstep(maximum*.25,maximum*.9,footprint);
 if(detail<=.001)return {field:average,queries:0};
 const q=forward(point),base=q.map(v=>Math.floor(v/stepLength));
 let nearest=1e6,radius=maximum*.5,queries=0;
 for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++){
  const id=[base[0]+x,base[1]+y,base[2]+z];
  const follicle=site(id);if(!follicle)continue;
  const {radius:r,centre}=follicle,d=Math.hypot(...q.map((v,k)=>v-centre[k]))/r;queries++;
  if(d<nearest){nearest=d;radius=r;}
 }
 const localDetail=detail*(1-smoothstep(radius*.5,radius*1.8,footprint));
 const edge=Math.min(.25,Math.max(.002,footprint/radius*.65)),tissue=1-smoothstep(1-edge,1+edge,nearest);
 const core=Math.min(tissue,1-smoothstep(1-shell-edge,1-shell+edge,nearest)),ring=Math.max(0,tissue-core);
 return {field:[mix(average[0],core,localDetail),mix(average[1],ring,localDetail),(.22*core+.5*ring)*localDetail],queries};
}

test('Only sourced follicle diameters and finite visual shell fractions are accepted',()=>{
 for(const thyroidFollicles of [{diameterMin:.000049},{diameterMax:.000501},{diameterMin:.0004,diameterMax:.0003},
  {diameterMin:NaN},{diameterMax:Infinity},{epithelialFraction:0},{epithelialFraction:.201}])
  assert.throws(()=>tissueWallMaterial(appearance({thyroidFollicles})),/Invalid thyroid follicle appearance/);
 assert.throws(()=>tissueWallMaterial(appearance({courseLength:.1})),/Invalid thyroid follicle appearance/,'A chart cannot replace the shared object field');
 assert.throws(()=>tissueWallMaterial(appearance({pattern:'solid',thyroidFollicles:{}})),/require thyroid-follicles/);
 const {material}=shaderFor();assert.equal(material.map,null);assert.equal(material.normalMap,null);material.dispose();
});

test('Physical diameters remain 50–500 µm independently of the appearance scale',()=>{
 const variants=[.0001,.0005,.01].map(scale=>shaderFor(appearance({scale}))),physicalPoints=[[.0001,.0002,.0003],[.002,.005,.003],[0,0,0]];
 for(const point of physicalPoints){
  const values=variants.map(({shader},i)=>fieldOracle(shader,point.map(v=>v/[.0001,.0005,.01][i]),.000001/[.0001,.0005,.01][i]).field);
  for(const other of values.slice(1))for(let k=0;k<3;k++)assert.ok(Math.abs(values[0][k]-other[k])<1e-7);
 }
 for(let i=0;i<variants.length;i++){
  const {shader,material}=variants[i],scale=[.0001,.0005,.01][i];
  assert.ok(Math.abs(Number(shader.fragmentShader.match(/const float minimum=([^;]+);/)[1])*scale-.00005)<1e-12);
  assert.ok(Math.abs(Number(shader.fragmentShader.match(/const float maximum=([^;]+);/)[1])*scale-.0005)<1e-12);material.dispose();
 }
});

test('A spherical cut has separate colloid core, epithelial ring and interfollicular tissue',()=>{
 const {shader,material}=shaderFor(),{site,reverse}=packing(shader),{centre,radius}=site([0,0,0]);
 const at=fraction=>fieldOracle(shader,reverse([centre[0]+radius*fraction,centre[1],centre[2]])).field;
 assert.deepEqual(at(.3),[1,0,.22]);assert.deepEqual(at(.97),[0,1,.5]);assert.deepEqual(at(1.15),[0,0,0]);
 assert.ok(shader.fragmentShader.includes('vec3 epithelium=color*vec3(.72,.60,.52)'));
 assert.ok(shader.fragmentShader.includes('vec3 colloid=mix(color*.95,vec3(.62,.54,.42),.3)'));
 assert.ok(shader.fragmentShader.includes('stroma*(1.0-field.x-field.y)'));
 material.dispose();
});

test('Subpixel follicles return stable mean color and zero relief before hash work',()=>{
 const {shader,material}=shaderFor();
 assert.ok(shader.fragmentShader.indexOf('if(detail<=.001)return average;')<shader.fragmentShader.indexOf('for(int x=-1'));
 const fields=[[0,0,0],[123.8,995,-335],[.4,.2,1.3]].map(p=>fieldOracle(shader,p,2));
 for(const value of fields){assert.equal(value.queries,0);assert.deepEqual(value.field,fields[0].field);assert.equal(value.field[2],0);}
 const [core,ring]=fields[0].field;assert.ok(core>0&&ring>0&&core+ring<1);
 // FCC sphere volume fixes the average; shell fraction splits core and lining.
 assert.ok(Math.abs(core/(core+ring)-.94**3)<1e-9);material.dispose();
});

test('Derivative filtering suppresses relief smoothly and preserves bounded tissue weights',()=>{
 const {shader,material}=shaderFor(),{site,reverse}=packing(shader),centre=reverse(site([0,0,0]).centre);
 let previous=Infinity;
 for(let width=0;width<=1;width+=.001){
  const {field}=fieldOracle(shader,centre,width);
  assert.ok(field.every(Number.isFinite));assert.ok(field[0]>=0&&field[1]>=0&&field[0]+field[1]<=1+1e-12);
  assert.ok(field[2]<=previous+1e-12);assert.ok(Math.abs(previous-field[2])<.012||previous===Infinity);previous=field[2];
 }
 assert.equal(previous,0);assert.ok(shader.fragmentShader.includes('max(length(dFdx(point)),length(dFdy(point)))'));material.dispose();
});

test('Variable-size follicles wander inside disjoint containers, including smaller void fillers',()=>{
 const {shader,material}=shaderFor(),p=packing(shader),sites=[];
 for(let x=-3;x<=3;x++)for(let y=-3;y<=3;y++)for(let z=-3;z<=3;z++)sites.push(p.site([x,y,z]));
 for(const s of sites){
  assert.ok(s.radius*2>=p.minimum-1e-12&&s.radius*2<=p.maximum+1e-12);
  assert.ok(s.travel+s.radius+p.padding<=s.container+1e-12);
 }
 for(let i=0;i<sites.length;i++)for(let j=i+1;j<sites.length;j++){
  const a=sites[i],b=sites[j],distance=Math.hypot(...a.centre.map((v,k)=>v-b.centre[k]));
  assert.ok(distance>=a.radius+b.radius-1e-12,'Actual follicle cores cannot cross');
 }
 const small=sites.filter(s=>!s.large),large=sites.filter(s=>s.large);
 assert.ok(small.length>150&&large.length>150);
 assert.ok(Math.max(...small.map(s=>s.radius))<Math.max(...large.map(s=>s.radius))*.5);
 assert.ok(Math.max(...large.map(s=>s.travel))>p.maximum*.25,'Small-radius follicles can leave the grid sites substantially');
 const axes=[[1,0,0],[0,1,0],[0,0,1]].map(v=>p.forward(v));
 for(const a of axes)assert.ok(Math.abs(Math.hypot(...a)-1)<1e-12);
 for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)assert.ok(Math.abs(axes[i].reduce((sum,v,k)=>sum+v*axes[j][k],0))<1e-12);
 material.dispose();
});

test('Surface and live section caps share the same object field and updatable coordinate frame',()=>{
 const frame=new Matrix4().makeTranslation(-.02,-1.55,.03),surface=shaderFor(),cap=shaderFor(appearance(),frame),inner=shaderFor(appearance(),undefined,true);
 const body=s=>s.fragmentShader.slice(s.fragmentShader.indexOf('vec3 thyroidFollicleHash'),s.fragmentShader.indexOf('varying vec3 tissuePoint;'));
 assert.equal(body(surface.shader),body(cap.shader));assert.equal(body(surface.shader),body(inner.shader));
 for(const value of [surface,cap,inner]){
  assert.ok(value.shader.vertexShader.includes('(tissueFrame*modelMatrix*vec4(position,1.0)).xyz:position)+tissueLocalOrigin'));
  assert.ok(value.shader.fragmentShader.includes('thyroidFollicleField(tissuePoint/0.00050000000)'));
  assert.ok(!value.shader.vertexShader.includes('attribute vec3 tissueChart'));assert.ok(!body(value.shader).includes('texture('));
 }
 assert.equal(surface.shader.uniforms.tissueWorld.value,0);assert.equal(cap.shader.uniforms.tissueWorld.value,1);
 assert.equal(cap.shader.uniforms.tissueFrame.value,frame);
 const key=cap.material.customProgramCacheKey();frame.makeTranslation(-.03,-1.57,.04);assert.equal(cap.material.customProgramCacheKey(),key,'Moving a cap updates its uniform without shader rebuilding');
 assert.ok(inner.shader.vertexShader.includes('objectNormal=-objectNormal;'));
 [surface,cap,inner].forEach(v=>v.material.dispose());
});

test('Changing follicle diameter or visual lining recompiles, while the other organ patterns keep their field',()=>{
 const values=[shaderFor(),shaderFor(appearance({thyroidFollicles:{diameterMin:.0001,diameterMax:.0004}})),shaderFor(appearance({thyroidFollicles:{epithelialFraction:.1}}))];
 assert.equal(new Set(values.map(v=>v.material.customProgramCacheKey())).size,3);
 const plain=shaderFor(appearance({pattern:'fibrous'}));assert.ok(!plain.shader.fragmentShader.includes('thyroidFollicleField'));
 assert.ok(plain.shader.fragmentShader.includes('.5+.30*sin(p.y*4.0+p.x*1.3+.2*sin(p.z*3.0))'));
 for(const {material} of [...values,plain])material.dispose();
});

test('Shared-bank rebasing retains the original tissue field and plain vascular owner normals',()=>{
 const origin=[0,1.5625,0],result=tissueWallMaterial(appearance(),undefined,false,origin);
 const shader={uniforms:{},vertexShader:ShaderLib.standard.vertexShader,fragmentShader:ShaderLib.standard.fragmentShader};
 result.material.onBeforeCompile(shader,{});
 assert.deepEqual(shader.uniforms.tissueLocalOrigin.value.toArray(),origin);
 const cap=tissueWallMaterial(appearance(),new Matrix4().makeTranslation(0,-1.5625,0),false,origin);
 const cut={uniforms:{},vertexShader:ShaderLib.standard.vertexShader,fragmentShader:ShaderLib.standard.fragmentShader};cap.material.onBeforeCompile(cut,{});
 assert.deepEqual(cut.uniforms.tissueLocalOrigin.value.toArray(),origin);
 assert.ok(cut.vertexShader.includes(':position)+tissueLocalOrigin;'));cap.material.dispose();
 const world=[.013,1.571,.009],local=world.map((v,i)=>v-origin[i]);
 assert.deepEqual(local.map((v,i)=>v+shader.uniforms.tissueLocalOrigin.value.toArray()[i]),world);
 assert.ok(shader.vertexShader.includes(':position)+tissueLocalOrigin;'));
 applyPairedNormalSigns(result.material);
 const signed={uniforms:{},vertexShader:ShaderLib.standard.vertexShader,fragmentShader:ShaderLib.standard.fragmentShader};
 result.material.onBeforeCompile(signed,{});
 assert.ok(signed.vertexShader.includes('objectNormal*=normalSign;'));
 assert.ok(result.material.customProgramCacheKey().endsWith(':paired-normal-signs'));
 result.material.dispose();
 const plain=tissueWallMaterial({color:[.7,.2,.15],pattern:'solid',scale:.001,relief:0}).material;
 plain.onBeforeCompile=()=>{};applyPairedNormalSigns(plain);
 const vessel={uniforms:{},vertexShader:ShaderLib.standard.vertexShader,fragmentShader:ShaderLib.standard.fragmentShader};
 plain.onBeforeCompile(vessel,{});
 assert.ok(vessel.vertexShader.includes('attribute float normalSign;'));
 assert.ok(vessel.vertexShader.includes('objectNormal*=normalSign;'));
 assert.equal(vessel.fragmentShader,ShaderLib.standard.fragmentShader);
 plain.dispose();
});
