import type {Texture} from 'three';
const masks=new WeakMap<Texture,{image:object;width:number;height:number;alpha:Uint8Array}>();
/** Cache alpha only; transparent texture gaps let picking reach actual skin. */
export function textureHitAlpha(texture:Texture,uv:{x:number;y:number}):number|undefined{
 const image=texture.image as HTMLImageElement|undefined;if(!image?.width||!image.height)return;
 let mask=masks.get(texture);
 if(!mask||mask.image!==image){
  const canvas=document.createElement('canvas');canvas.width=image.width;canvas.height=image.height;
  const context=canvas.getContext('2d',{willReadFrequently:true});if(!context)return;
  context.drawImage(image,0,0);const rgba=context.getImageData(0,0,canvas.width,canvas.height).data;
  const alpha=new Uint8Array(canvas.width*canvas.height);for(let i=0;i<alpha.length;i++)alpha[i]=rgba[i*4+3];
  mask={image,width:canvas.width,height:canvas.height,alpha};masks.set(texture,mask);
 }
 const u=uv.x-Math.floor(uv.x),v=uv.y-Math.floor(uv.y);const x=Math.min(mask.width-1,Math.floor(u*mask.width)),y=Math.min(mask.height-1,Math.floor((texture.flipY?1-v:v)*mask.height));
 return mask.alpha[y*mask.width+x]/255;
}
