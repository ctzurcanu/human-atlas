/** Keep saved-view indices while delivering the exact lossless runtime meshes.
 * Suppressed historical records retain their original buffer references but
 * never cause source chunks to be requested by the on-demand viewer. */
export function orderedRuntimeCatalogue(source,runtime,sourceSha256) {
 if(runtime.delivery?.mode!=='on-demand'||runtime.delivery.sourceAtlasSha256!==sourceSha256)throw Error('Runtime catalogue is not bound to the current source');
 const packed=new Map(runtime.parts.map(p=>[p.id,p]));
 if(packed.size!==runtime.parts.length)throw Error('Duplicate runtime part');
 const originalIds=new Set(source.parts.map(p=>p.id));
 if([...packed.keys()].some(id=>!originalIds.has(id)))throw Error('Runtime contains an unrelated mesh');
 if(JSON.stringify(source.materials)!==JSON.stringify(runtime.materials))throw Error('Runtime materials differ');
 const result=structuredClone(source),chunks=structuredClone(runtime.chunks),historicalChunks=new Map();
 result.parts=source.parts.map(p=>{
  const delivered=packed.get(p.id);
  if(delivered){if(p.suppressed||delivered.suppressed)throw Error('Runtime contains suppressed anatomy');return structuredClone(delivered);}
  if(!p.suppressed)throw Error('Runtime omitted visible anatomy '+p.id);
  let chunk=historicalChunks.get(p.chunk);
  if(chunk===undefined){chunk=chunks.length;chunks.push(structuredClone(source.chunks[p.chunk]));historicalChunks.set(p.chunk,chunk);}
  return {...structuredClone(p),chunk};
 });
 result.chunks=chunks;result.delivery=structuredClone(runtime.delivery);
 return result;
}
