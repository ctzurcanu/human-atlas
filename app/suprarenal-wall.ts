import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';

type Gland={side:'left'|'right';sourcePartId:string;wholeConceptId:string;partIds:string[]};
type Marker={schema:string;sex:string;partIds:string[];glands:Gland[]};
type Registry={marker:Marker;glands:{definition:Gland;parts:Part[]}[]};
const registries=new WeakMap<Atlas,Registry|null>();
const same=(a:readonly string[],b:readonly string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));
const roles=['capsule','cortex','medulla'];

/** Recover the native whole's semantic identity only for its complete, explicit
 * replacement. A cortex or medulla leaf does not inherit whole-organ routes. */
export function suprarenalWholeSource(selected:Part[],parts:Part[]):Part|undefined{
 if(selected.length!==3||new Set(selected.map(p=>p.id)).size!==3)return;
 const first=selected[0],sourceId=(first as Part&{adrenalSourcePartId?:string}).adrenalSourcePartId;
 const match=/^TA98PART:(female|male):Suprarenal:(left|right):(?:capsule|cortex|medulla)$/.exec(first.id);
 if(!match||!sourceId)return;
 const [,sex,side]=match,expected=roles.map(role=>`TA98PART:${sex}:Suprarenal:${side}:${role}`),source=parts.find(p=>p.id===sourceId);
 const native=sex==='female'?`LOCAL:female:${side[0]}_adrenal_gland:${side==='left'?'190':'191'}`:`ZA:Suprarenal gland.${side[0]}`;
 if(source?.id!==native||!source.suppressed||source.system!=='endocrine'||source.sourceOffset||!same(source.replacementParts??[],expected))return;
 if(!same(selected.map(p=>p.id),expected)||selected.some(p=>p.suppressed||p.sourceOffset||p.system!=='endocrine'
  ||(p as Part&{adrenalSourcePartId?:string}).adrenalSourcePartId!==sourceId||!p.pairedTissue||!p.sectionSolid
  ||p.tissueWallAssembly!==`suprarenal:${sex}:${side}`||p.normalSigns===undefined
  ||p.id!==`TA98PART:${sex}:Suprarenal:${side}:${p.tissueWallLayer}`||p.chunk!==first.chunk||p.positions!==first.positions
  ||p.normals!==first.normals||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return;
 return source;
}

export function suprarenalWallRegistry(atlas:Atlas):Registry|null{
 if(registries.has(atlas))return registries.get(atlas)!;
 const reject=()=>{registries.set(atlas,null);return null;};
 const marker=(atlas as Atlas&{finiteSuprarenalIntegration?:Marker}).finiteSuprarenalIntegration;
 if(marker?.schema!=='suprarenal-local-tissues-v223'||marker.sex!==atlas.sex||!['male','female'].includes(marker.sex)
  ||!Array.isArray(marker.glands)||marker.glands.length!==2||!same(marker.glands.map(g=>g.side),['left','right'])
  ||!Array.isArray(marker.partIds)||!same(marker.partIds,marker.glands.flatMap(g=>g.partIds??[])))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c])),glands=[];
 for(const g of marker.glands){
  const expected=roles.map(role=>`TA98PART:${marker.sex}:Suprarenal:${g.side}:${role}`);
  if(!Array.isArray(g.partIds)||!same(g.partIds,expected))return reject();
  const members=g.partIds.map(id=>byId.get(id)),first=members[0],source=byId.get(g.sourcePartId),whole=concepts.get(g.wholeConceptId);
  if(!first||!source?.suppressed||!same(source.replacementParts??[],expected)||!whole||!same(whole.elements,expected))return reject();
  if(members.some(p=>!p||p.suppressed||p.sourceOffset||p.system!=='endocrine'||!p.pairedTissue||!p.sectionSolid||!p.tissueAppearance
   ||p.normalSigns===undefined||p.tissueWallAssembly!==`suprarenal:${marker.sex}:${g.side}`||!roles.includes(p.tissueWallLayer??'')
   ||p.id!==`TA98PART:${marker.sex}:Suprarenal:${g.side}:${p.tissueWallLayer}`||p.chunk!==first.chunk||p.positions!==first.positions
   ||p.normals!==first.normals||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return reject();
  const parts=members as Part[];
  for(const p of parts){
   let end=0;
   for(const group of p.tissueDrawGroups??[]){
    if(group.start!==end||!Number.isInteger(group.count)||group.count<=0||group.count%3)return reject();
    end+=group.count;
    if(group.neighbor){const peer=parts.find(q=>q.id===group.neighbor);if(!peer||peer.id===p.id
     ||peer.tissueDrawGroups?.filter(q=>q.neighbor===p.id).reduce((n,q)=>n+q.count,0)!==group.count)return reject();}
   }
   if(end!==p.indexCount)return reject();
  }
  glands.push({definition:g,parts});
 }
 const value={marker,glands};registries.set(atlas,value);return value;
}

/** Only complete registered gland tissues peel; vessels are separate anatomy. */
export function suprarenalTissuePeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const r=suprarenalWallRegistry(atlas);if(!r||choice.id==='selection-set')return null;
 const selected=new Set(choice.elements),glands=r.glands.filter(g=>g.parts.every(p=>selected.has(p.id)));
 if(!glands.length||selected.size!==glands.length*3)return null;
 return roles.map(role=>{const parts=glands.flatMap(g=>g.parts.filter(p=>p.tissueWallLayer===role));return {name:parts.length===1?parts[0].name:role==='capsule'?'Fibrous capsules of suprarenal glands':role==='cortex'?'Suprarenal cortex':'Suprarenal medulla',ids:parts.map(p=>p.id)};});
}
