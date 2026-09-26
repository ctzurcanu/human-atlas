/** Deterministic, exactly tileable cut-face illustrations. Run: npm run build:tissue-textures */
import {mkdir} from 'node:fs/promises';
import path from 'node:path';
import sharp from 'sharp';

const size=256,root=path.resolve('public/models');
await mkdir(root,{recursive:true});
const fract=x=>x-Math.floor(x),clamp=(v,a=0,b=1)=>Math.max(a,Math.min(b,v));
const hash=(x,y,seed)=>fract(Math.sin(x*127.1+y*311.7+seed*74.7)*43758.5453);
function noise(u,v,cells,seed){
 const x=u*cells,y=v*cells,ix=Math.floor(x),iy=Math.floor(y),fx=fract(x),fy=fract(y),sx=fx*fx*(3-2*fx),sy=fy*fy*(3-2*fy);
 const h=(a,b)=>hash((a+cells)%cells,(b+cells)%cells,seed);
 return (h(ix,iy)*(1-sx)+h(ix+1,iy)*sx)*(1-sy)+(h(ix,iy+1)*(1-sx)+h(ix+1,iy+1)*sx)*sy;
}
function cellDistance(u,v,cells,seed){
 let nearest=2;
 const x=u*cells,y=v*cells,ix=Math.floor(x),iy=Math.floor(y);
 for(let dy=-1;dy<=1;dy++)for(let dx=-1;dx<=1;dx++){
  const gx=ix+dx,gy=iy+dy,cx=gx+.25+.5*hash((gx+cells)%cells,(gy+cells)%cells,seed),cy=gy+.25+.5*hash((gx+cells)%cells,(gy+cells)%cells,seed+1);
  nearest=Math.min(nearest,Math.hypot(x-cx,y-cy));
 }
 return nearest;
}
const specs={
 'muscle-surface-seamless.png':(u,v)=>{
  const n=noise(u,v,8,2),f=Math.sin(2*Math.PI*(u*38+.16*Math.sin(v*2*Math.PI*5))),band=.5+.5*f;
  return [154+29*band+13*n,72+19*band+11*n,70+18*band+9*n];
 },
 'muscle-axial-seamless.png':(u,v)=>{
  const d=cellDistance(u,v,18,3),ring=Math.exp(-(((d-.37)/.085)**2)),core=Math.exp(-((d/.24)**2)),n=noise(u,v,15,4);
  return [190-39*ring+13*core+13*n,104-36*ring+12*core+9*n,99-31*ring+10*core+8*n];
 },
 'muscle-coronal-seamless.png':(u,v)=>{
  const n=noise(u,v,9,5),stripe=.5+.5*Math.sin(2*Math.PI*(u*34+.11*Math.sin(v*2*Math.PI*4))),streak=.5+.5*Math.sin(2*Math.PI*(u*11+v*3));
  return [164+37*stripe+11*streak+9*n,78+34*stripe+8*streak+8*n,77+31*stripe+7*streak+6*n];
 },
 'muscle-sagittal-seamless.png':(u,v)=>{
  const n=noise(u,v,8,6),stripe=.5+.5*Math.sin(2*Math.PI*(u*25+v*5+.12*Math.sin(v*2*Math.PI*3))),streak=.5+.5*Math.sin(2*Math.PI*(u*9-v*2));
  return [163+34*stripe+10*streak+10*n,80+28*stripe+9*streak+8*n,80+25*stripe+8*streak+7*n];
 },
 'lung-cut-seamless.png':(u,v)=>{
  const d=cellDistance(u,v,23,7),edge=Math.exp(-(((d-.39)/.08)**2)),air=Math.exp(-((d/.22)**2)),n=noise(u,v,11,8);
  return [214-35*edge+12*air+11*n,146-43*edge+13*air+10*n,151-39*edge+13*air+11*n];
 },
 'liver-cut-seamless.png':(u,v)=>{
  const a=noise(u,v,7,9),b=noise(u,v,19,10),sinus=.5+.5*Math.sin(2*Math.PI*(u*8+v*5+.13*a));
  return [125+27*a+16*sinus+13*b,57+17*a+7*sinus+8*b,63+15*a+7*sinus+7*b];
 },
 'spleen-cut-seamless.png':(u,v)=>{
  const a=noise(u,v,8,17),b=noise(u,v,21,18),trabecula=.5+.5*Math.sin(2*Math.PI*(u*9+v*4+.15*a));
  return [204+24*a+9*trabecula+7*b,165+24*a+10*trabecula+7*b,92+17*a+7*trabecula+6*b];
 },
 'cns-cut-seamless.png':(u,v)=>{
  const a=noise(u,v,8,11),b=noise(u,v,22,12),tract=.5+.5*Math.sin(2*Math.PI*(u*7+v*3+.16*a));
  return [176+33*a+16*tract+6*b,165+32*a+17*tract+6*b,153+32*a+17*tract+6*b];
 },
 'organ-cut-seamless.png':(u,v)=>{
  const a=noise(u,v,9,13),b=noise(u,v,26,14),grain=.5+.5*Math.sin(2*Math.PI*(u*13+v*8+.13*a));
  return [177+24*a+9*grain+7*b,113+23*a+8*grain+6*b,113+22*a+8*grain+6*b];
 },
};
for(const [filename,sample] of Object.entries(specs)){
 const raw=Buffer.alloc(size*size*3);
 for(let y=0;y<size;y++)for(let x=0;x<size;x++){
  const color=sample(x/(size-1),y/(size-1));
  for(let channel=0;channel<3;channel++)raw[(y*size+x)*3+channel]=Math.round(clamp(color[channel],0,255));
 }
 // Fix floating-point noise at the seams as well as the periodic sampling.
 for(let y=0;y<size;y++)for(let channel=0;channel<3;channel++)raw[(y*size+size-1)*3+channel]=raw[y*size*3+channel];
 for(let x=0;x<size;x++)for(let channel=0;channel<3;channel++)raw[((size-1)*size+x)*3+channel]=raw[x*3+channel];
 await sharp(raw,{raw:{width:size,height:size,channels:3}}).png().toFile(path.join(root,filename));
 console.log(filename);
}
