export interface TextureImage {
 image:ImageBitmap|HTMLImageElement;
 flipY:boolean;
 dispose:()=>void;
}
/** Decode fetched images without blob: URLs, which host sandboxes can block. */
export async function decodeTextureImage(blob:Blob,flipY:boolean):Promise<TextureImage>{
 if(typeof createImageBitmap==='function'){
  const image=await createImageBitmap(blob,{imageOrientation:flipY?'flipY':'none',premultiplyAlpha:'none',colorSpaceConversion:'none'});
  // Bitmap orientation is already baked in. Picking must use that same image.
  return {image,flipY:false,dispose:()=>image.close()};
 }
 const dataUrl=await new Promise<string>((resolve,reject)=>{
  const reader=new FileReader();reader.onload=()=>resolve(reader.result as string);reader.onerror=()=>reject(reader.error);reader.readAsDataURL(blob);
 });
 const image=await new Promise<HTMLImageElement>((resolve,reject)=>{
  const image=new Image();image.onload=()=>resolve(image);image.onerror=()=>reject(new Error('Texture image could not be decoded.'));image.src=dataUrl;
 });
 return {image,flipY,dispose:()=>{}};
}
