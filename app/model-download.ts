import {gunzipSync} from 'fflate';
import {fetchModelAsset,forgetAsset,isCachedAsset,rememberDecodedModel} from './asset-cache';

/** Static hosts may serve .gz as a compressed response or as a gzip file.
 * Fetch already decodes Content-Encoding; inspect the payload to avoid decoding twice.
 */
export async function decodeModelResponse(response:Response,expectedBytes:number,compressed:boolean):Promise<ArrayBuffer>{
 if(!response.ok)throw new Error('An anatomy file could not be loaded.');
 const payload=await response.arrayBuffer(),signature=new Uint8Array(payload,0,Math.min(2,payload.byteLength));
 const gzip=compressed&&signature[0]===0x1f&&signature[1]===0x8b;
 const buffer=gzip?(typeof DecompressionStream!=='undefined'?await new Response(new Blob([payload]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer():gunzipSync(new Uint8Array(payload)).slice().buffer as ArrayBuffer):payload;
 if(buffer.byteLength!==expectedBytes)throw new Error('An anatomy file was incomplete. Please reload the viewer.');
 return buffer;
}

/** Discard a damaged cached chunk and retry it once from the source. */
export async function loadModelBuffer(url:string,expectedBytes:number,compressed:boolean,signal:AbortSignal):Promise<ArrayBuffer>{
 const response=await fetchModelAsset(url,{signal});
 if(response instanceof ArrayBuffer){if(response.byteLength===expectedBytes)return response;await forgetAsset(url);const fresh=await fetchModelAsset(url,{signal,cache:'reload'});if(fresh instanceof ArrayBuffer)throw new Error('An anatomy file was incomplete. Please reload the viewer.');const buffer=await decodeModelResponse(fresh,expectedBytes,compressed);if(signal.aborted)throw signal.reason??new DOMException('Aborted','AbortError');await rememberDecodedModel(url,buffer);return buffer;}
 try{const buffer=await decodeModelResponse(response,expectedBytes,compressed);if(signal.aborted)throw signal.reason??new DOMException('Aborted','AbortError');if(response.headers.get('x-atlas-decoded')!=='1')await rememberDecodedModel(url,buffer);return buffer;}
 catch(error){
  if(signal.aborted)throw error;
  await forgetAsset(url);
  if(!isCachedAsset(response))throw error;
  const fresh=await fetchModelAsset(url,{signal,cache:'reload'});
  if(fresh instanceof ArrayBuffer)throw new Error('An anatomy file could not be reloaded.');
  try{const buffer=await decodeModelResponse(fresh,expectedBytes,compressed);if(signal.aborted)throw signal.reason??new DOMException('Aborted','AbortError');await rememberDecodedModel(url,buffer);return buffer;}catch(error){await forgetAsset(url);throw error;}
 }
}
