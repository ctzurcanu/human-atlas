import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {gallbladderRelations} from './gallbladder-relations';
import {gallbladderWholeSource,gallbladderWallInterface} from './gallbladder-wall';
import {nativePartFrame} from './ta98-modeled-parent';

const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const independent=(p:Part)=>!!p.sourceOffset||!!(p.provenance&&'mirrored' in p.provenance&&p.provenance.mirrored===true);
/** PMID26240761: cystic arterial branches supply the duct.
 * PMC4319126: hepatic plexus routes accompany it.
 * NBK459288: the cystic node is a hepatocystic-triangle landmark.
 * These references do not create intramural geometry or verified attachments. */
export function cysticFunctionalRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(duct:Part,target:Part,kind:ResolvedRelation['kind'],reverse:ResolvedRelation['kind'],note:string,via?:string[])=>{
  for(const [a,b,k] of [[duct,target,kind],[target,duct,reverse]] as const){
   const rows=graph.get(a.id)??[];rows.push({kind:k,target:b,note,...(via?.length?{via,viaModeled:false}: {})});graph.set(a.id,rows);
  }
 };
 for(const sex of ['male','female'] as const){
  const duct=parts.find(p=>p.id===(sex==='male'?'FJ3080':'LOCAL:female:cystic_duct:authored'));
  if(!duct)continue;
  const coats=parts.filter(p=>!p.suppressed&&p.tissueWallAssembly===`gallbladder:${sex}`);
  const source=gallbladderWholeSource(coats,parts);
  if(!source||coats.filter(p=>gallbladderWallInterface(p,duct)).length!==4)continue;
  const reference=parts.map(p=>p===source?{...p,suppressed:false}:p);
  for(const row of gallbladderRelations(reference).get(source.id)??[]){
   if(independent(row.target))continue;
   if(row.kind==='arterial')add(duct,row.target,'arterial','supplies',
    `${row.target.name} is modeled. Its small arterial branch into the cystic-duct wall is not modeled.`,
    [...(row.via??[]),'Arterial branch to cystic duct']);
   if(row.kind==='innervation')add(duct,row.target,'innervation','innervates',
    `${row.target.name} is modeled. The hepatic plexus and its cystic-duct fibers are not modeled.`, ['Hepatic plexus']);
   if(row.kind==='lymphaticDrainage')add(duct,row.target,'adjacent','adjacent',
    'The cystic lymph-node mesh is present. Its position relative to this duct has not been verified; the connecting lymphatic vessels are not modeled.');
  }
  const frame=nativePartFrame(source);
  for(const nerve of parts){
   if(nerve.suppressed||nerve.system!=='nervous'||independent(nerve)||!nerve.sourceId||nativePartFrame(nerve)!==frame)continue;
   const native=nerve.id===`ZA:${nerve.sourceId}`||frame?.startsWith('LOCAL:')&&nerve.id.startsWith(`${frame}${nerve.sourceId}:`);
   if(!native||!/^Sympathetic nerves \((left|right)\)$/.test(nerve.name))continue;
   add(duct,nerve,'innervation','innervates',`${nerve.name} is modeled. The coeliac plexus, hepatic plexus and duct-wall fibers are not modeled.`, ['Coeliac plexus','Hepatic plexus']);
  }
 }
 cache.set(parts,graph);return graph;
}
