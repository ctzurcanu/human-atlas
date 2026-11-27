import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';

type Point=[number,number,number];
type Boundary={name:string;closed:boolean;points:Point[];samplePartIds:string[];faceIds:number[];barycentric:Point[];sourceBankSha256:string};
type TermGroup={term:string;conceptId:string;partIds:string[];representation:'whole-organ'|'tissue-assembly'|'finite-tissue'|'boundary-reference';boundaries?:Boundary[]};
type Marker={schema:string;partIds:string[];bankSha256:string;sourcePlanSha256:string;termGroups:TermGroup[]};
type Registry={marker:Marker;byId:Map<string,Part>};
const cache=new WeakMap<Atlas,Registry|null>();
const sha=(s:unknown)=>typeof s==='string'&&/^[a-f0-9]{64}$/.test(s);
const same=(a:readonly string[],b:readonly string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));
const referenceTerms=new Set(['A04.4.02.010','A04.4.02.011','A04.4.02.014','A04.4.02.015','A04.4.02.016']);
const terms=Array.from({length:17},(_,i)=>`A04.4.02.${String(i+1).padStart(3,'0')}`);
const lumbarTerms=new Set(terms.slice(1,7));

/** Openings and muscular triangles are source-bound references. Their tissue
 * hosts provide selectable context; no solid is inserted into a passage. */
export function diaphragmControlRegistry(atlas:Atlas):Registry|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const m=(atlas as Atlas&{finiteDiaphragmIntegration?:Marker}).finiteDiaphragmIntegration;
 if(!m||m.schema!=='diaphragm-complete-local-v7'||!sha(m.bankSha256)||!sha(m.sourcePlanSha256)
  ||!Array.isArray(m.partIds)||!m.partIds.length||!Array.isArray(m.termGroups)||!same(m.termGroups.map(g=>g?.term),terms))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 if(new Set(m.partIds).size!==m.partIds.length)return reject();
 for(const id of m.partIds){const p=byId.get(id);
  if(!p||p.suppressed||!p.pairedTissue||!p.sectionSolid||p.sourceOffset||!p.normalSigns
   ||(atlas.chunks[p.chunk]as{sha256?:string}|undefined)?.sha256!==m.bankSha256)return reject();
 }
 for(const g of m.termGroups){const c=concepts.get(g.conceptId);
  if(!c||c.ta98Kind!=='exact'||c.ta98Term!==g.term||!Array.isArray(g.partIds)||!g.partIds.length
   ||!same(c.elements,g.partIds)||g.partIds.some(id=>!m.partIds.includes(id)))return reject();
  if(g.term==='A04.4.02.001'){
   if(g.representation!=='whole-organ'||!same(g.partIds,m.partIds))return reject();
  }else if(g.term==='A04.4.02.002'){
   const lumbar=m.partIds.filter(id=>lumbarTerms.has(byId.get(id)?.ta98Term??'')||id==='TA98PART:ZA:Diaphragm:Tendon:crural-attachment');
   if(g.representation!=='tissue-assembly'||g.boundaries||!same(g.partIds,lumbar))return reject();
  }else if(!referenceTerms.has(g.term)){
   if(g.representation!=='finite-tissue'||g.boundaries||g.partIds.some(id=>byId.get(id)?.ta98Term!==g.term))return reject();
  }else{
   if(g.representation!=='boundary-reference'||!Array.isArray(g.boundaries)||!g.boundaries.length
    ||(['A04.4.02.015','A04.4.02.016'].includes(g.term)&&!same(g.boundaries.map(b=>b.name),['left','right'])))return reject();
   for(const b of g.boundaries){const n=b.points?.length;
    if(!Number.isInteger(n)||n<(b.closed?3:2)||n>20000||b.sourceBankSha256!==m.bankSha256
     ||[b.samplePartIds,b.faceIds,b.barycentric].some(a=>!Array.isArray(a)||a.length!==n))return reject();
    for(let i=0;i<n;i++){const p=byId.get(b.samplePartIds[i]),q=b.points[i],w=b.barycentric[i],face=b.faceIds[i];
     if(!p||!g.partIds.includes(p.id)||!Number.isInteger(face)||face<0||face>=p.indexCount/3
      ||!Array.isArray(q)||q.length!==3||q.some((x,k)=>!Number.isFinite(x)||x<p.bounds[0][k]-1e-12||x>p.bounds[1][k]+1e-12)
      ||!Array.isArray(w)||w.length!==3||w.some(x=>!Number.isFinite(x)||x<0||x>1)||Math.abs(w.reduce((a,b)=>a+b,0)-1)>1e-12)return reject();
    }
   }
  }
 }
 const r={marker:m,byId};cache.set(atlas,r);return r;
}

export function diaphragmTermParts(atlas:Atlas):Map<string,Part[]>{
 const r=diaphragmControlRegistry(atlas);
 return new Map(r?.marker.termGroups.map(g=>[g.term,g.partIds.map(id=>r.byId.get(id)!)])??[]);
}

/** Whole-organ controls select the stored assembly in every browser. A named
 * tissue, including the crural attachment with the whole-organ TA98 code,
 * remains separately selectable. */
export function diaphragmAssemblyChoice(atlas:Atlas,choice:HierarchyChoice):HierarchyChoice{
 const r=diaphragmControlRegistry(atlas),whole=r?.marker.termGroups.find(g=>g.term==='A04.4.02.001');
 if(!whole||choice.id!==whole.conceptId&&(!choice.children||(choice.ta98Term??choice.terminology?.ta98)!==whole.term))return choice;
 return {...choice,elements:[...whole.partIds]};
}

export function diaphragmLabelReference(atlas:Atlas,choice:HierarchyChoice):{partId:string;point:Point}|undefined{
 const r=diaphragmControlRegistry(atlas);if(!r)return;
 const term=choice.ta98Term??choice.terminology?.ta98;
 const g=r.marker.termGroups.find(g=>g.representation==='boundary-reference'&&(g.conceptId===choice.id||g.term===term));
 if(!g||!same(choice.elements,g.partIds))return;
 const b=g.boundaries![0];
 // Anchor to a real rim point, never to a centroid floating in the opening.
 const i=b.points.reduce((best,p,k)=>p[2]>b.points[best][2]?k:best,0);
 return {partId:b.samplePartIds[i],point:[...b.points[i]]};
}
