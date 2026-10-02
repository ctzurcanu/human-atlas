interface ChunkRetention {maxIdleBytes:number;bytes:(id:number)=>number;release:(id:number)=>void}
/** Bounded geometry queue with optional byte-limited inactive retention. */
export class OnDemandChunks {
 readonly loaded=new Set<number>();
 private required:number[]=[];
 private active=new Map<number,AbortController>();
 private failed=new Set<number>();
 private stopped=false;
 private limit:number;
 private load:(id:number,signal:AbortSignal)=>Promise<void>;
 private changed:()=>void;
 private error:(error:unknown)=>void;
 private retention?:ChunkRetention;
 private used=new Map<number,number>();
 private clock=0;
 constructor(limit:number,load:(id:number,signal:AbortSignal)=>Promise<void>,changed:()=>void,error:(error:unknown)=>void,retention?:ChunkRetention){
  this.limit=limit;this.load=load;this.changed=changed;this.error=error;
  this.retention=retention;
 }
 update(ids:number[]){
  if(this.stopped)return;
  this.required=[...new Set(ids)];
  for(const id of this.required)this.used.set(id,++this.clock);
  for(const [id,controller] of this.active)if(!this.required.includes(id))controller.abort();
  this.prune();
  this.pump();this.changed();
 }
 get pending(){return this.required.filter(id=>!this.loaded.has(id)).length;}
 get total(){return this.required.length;}
 get idleBytes(){return this.retention?[...this.loaded].filter(id=>!this.required.includes(id)).reduce((sum,id)=>sum+this.retention!.bytes(id),0):0;}
 private prune(){
  if(!this.retention)return;
  let bytes=this.idleBytes;
  const idle=[...this.loaded].filter(id=>!this.required.includes(id)).sort((a,b)=>(this.used.get(a)??0)-(this.used.get(b)??0));
  for(const id of idle){if(bytes<=this.retention.maxIdleBytes)break;bytes-=this.retention.bytes(id);this.retention.release(id);this.loaded.delete(id);this.used.delete(id);}
 }
 stop(){this.stopped=true;for(const controller of this.active.values())controller.abort();}
 private pump(){
  if(this.stopped)return;
  for(const id of this.required){
   if(this.active.size>=this.limit)break;
   if(this.loaded.has(id)||this.active.has(id)||this.failed.has(id))continue;
   const controller=new AbortController();this.active.set(id,controller);
   void this.load(id,controller.signal).then(()=>{
    // A cancelled fetch may still finish decoding. Never install it as ready.
    if(!controller.signal.aborted&&!this.stopped){this.loaded.add(id);this.used.set(id,++this.clock);this.prune();}
    else if(!this.stopped)this.retention?.release(id);
   }).catch(error=>{if(!controller.signal.aborted&&!this.stopped){this.failed.add(id);this.error(error);}}).finally(()=>{
    this.active.delete(id);if(!this.stopped){this.changed();this.pump();}
   });
  }
 }
}
