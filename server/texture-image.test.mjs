import {test} from 'node:test';
import assert from 'node:assert/strict';
import {build} from 'esbuild';
const bundle=await build({entryPoints:['app/texture-image.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {decodeTextureImage}=await import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles[0].text).toString('base64')}`);

test('bitmap decoding preserves both UV orientations and releases decoded image memory',async()=>{
 const original=globalThis.createImageBitmap,calls=[];
 globalThis.createImageBitmap=async(blob,options)=>{
  const image={width:2,height:2,closed:false,close(){this.closed=true;}};
  calls.push({blob,options,image});return image;
 };
 try{
  const blob=new Blob(['image']);
  for(const flipY of [true,false]){
   const decoded=await decodeTextureImage(blob,flipY),call=calls.at(-1);
   assert.equal(call.blob,blob);
   assert.equal(call.options.imageOrientation,flipY?'flipY':'none');
   assert.equal(call.options.premultiplyAlpha,'none');
   assert.equal(call.options.colorSpaceConversion,'none');
   assert.equal(decoded.flipY,false,'Picking uses the bitmap orientation already applied during decode');
   assert.equal(decoded.image,call.image);decoded.dispose();assert.equal(call.image.closed,true);
  }
 }finally{globalThis.createImageBitmap=original;}
});

test('older browsers decode a data image without creating a blob URL',async()=>{
 const saved={createImageBitmap:globalThis.createImageBitmap,FileReader:globalThis.FileReader,Image:globalThis.Image};
 let imageUrl;
 globalThis.createImageBitmap=undefined;
 globalThis.FileReader=class {readAsDataURL(){this.result='data:image/png;base64,iVBORw0KGgo=';this.onload();}};
 globalThis.Image=class {set src(value){imageUrl=value;this.onload();}};
 try{
  const decoded=await decodeTextureImage(new Blob(['image']),true);
  assert.ok(imageUrl.startsWith('data:image/png;'));assert.equal(decoded.flipY,true);decoded.dispose();
 }finally{Object.assign(globalThis,saved);}
});
