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

/** Retry a damaged cache entry or a transient fetch failure once, without resetting the view. */
export async function loadModelBuffer(url:string,expectedBytes:number,compressed:boolean,signal:AbortSignal):Promise<ArrayBuffer>{
 for(let attempt=0;attempt<2;attempt++){
  let cached=false;
  try{
   if(signal.aborted)throw signal.reason??new DOMException('Aborted','AbortError');
   const response=await fetchModelAsset(url,{signal,...(attempt?{cache:'reload' as const}:{})});
   cached=response instanceof ArrayBuffer||isCachedAsset(response);
   const buffer=response instanceof ArrayBuffer?response:await decodeModelResponse(response,expectedBytes,compressed);
   if(buffer.byteLength!==expectedBytes)throw new Error('An anatomy file was incomplete. Please reload the viewer.');
   if(signal.aborted)throw signal.reason??new DOMException('Aborted','AbortError');
   if(!(response instanceof ArrayBuffer)&&response.headers.get('x-atlas-decoded')!=='1')await rememberDecodedModel(url,buffer);
   return buffer;
  }catch(error){
   if(signal.aborted)throw error;
   await forgetAsset(url);
   // Fetch rejects with TypeError in Chromium/WebKit and NetworkError in Firefox.
   // A body stream can fail after the headers, so decoding belongs in this retry too.
   const network=error instanceof TypeError||error instanceof Error&&error.name==='NetworkError';
   if(!attempt&&(cached||network))continue;
   throw new Error(`Could not load anatomy geometry: ${url} (${error instanceof Error?error.message:String(error)})`,{cause:error});
  }
 }
 throw new Error(`Could not load anatomy geometry: ${url}`);
}
