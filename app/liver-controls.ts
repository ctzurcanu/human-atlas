import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';

type LiverGroup={term:string;conceptId:string;parts:string[]};
type LiverMarker={schema:string;partIds:string[];wholeConceptId:string;bankSha256:string;bankUrl:string;termGroups:LiverGroup[]};
type Registry={marker:LiverMarker;byId:Map<string,Part>};
const cache=new WeakMap<Atlas,Registry|null>();
const same=(a:string[],b:string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));

/** Explicit local assembly. Open hilar references are not certified solid
 * tissue, and named functional segments are not inferred from the gross lobes. */
export function liverControlRegistry(atlas:Atlas):Registry|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const m=(atlas as Atlas&{liverIntegration?:LiverMarker}).liverIntegration;
 const sex=atlas.sex,wholeId=sex==='female'?'LOCAL:female:liver':'ZA:Liver';
 const expanded=m?.schema===`${sex}-liver-lobes-capsule/v2`;
 const names=['right','left','quadrate','caudate','capsule',...(expanded?['serosa','subserosa','fibrous-appendix','papillary-process','caudate-process']:[])];
 if(!['male','female'].includes(sex??'')||!m||![`${sex}-liver-lobes-capsule/v1`,`${sex}-liver-lobes-capsule/v2`].includes(m.schema)||m.wholeConceptId!==wholeId
  ||!Array.isArray(m.partIds)||!same(m.partIds,names.map(n=>'TA98PART:ZA:Liver:'+n))
  ||!/^[a-f0-9]{64}$/.test(m.bankSha256)||!Array.isArray(m.termGroups)||m.termGroups.length!==(expanded?11:6))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 for(const id of m.partIds){const p=byId.get(id);
  if(!p||p.suppressed||p.defaultHidden||p.system!=='digestive'||p.conceptId!==id||p.sourceOffset
   ||p.sourceId!=='Liver:'+id.split(':').at(-1)||p.tissueWallAssembly!==`liver:${sex}:lobes-v1`
   ||!p.pairedTissue||p.normalSigns===undefined||atlas.chunks[p.chunk]?.url!==m.bankUrl)return reject();
 }
 const expected=['A05.8.01.001','A05.8.01.026','A05.8.01.027','A05.8.01.029','A05.8.01.030','A05.8.01.054',...(expanded?[28,31,32,52,53].map(n=>`A05.8.01.${String(n).padStart(3,'0')}`):[])];
 if(!same(m.termGroups.map(g=>g.term),expected))return reject();
 for(const g of m.termGroups){const c=concepts.get(g.conceptId);
  if(!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||!same(c.elements,g.parts)||g.parts.some(id=>!m.partIds.includes(id)))return reject();
  const children=expanded&&g.term==='A05.8.01.027'?['left','fibrous-appendix']:expanded&&g.term==='A05.8.01.030'?['caudate','papillary-process','caudate-process']:undefined;
  if(g.term==='A05.8.01.001'?!same(g.parts,m.partIds):children?!same(g.parts,children.map(n=>'TA98PART:ZA:Liver:'+n)):g.parts.length!==1||byId.get(g.parts[0])?.ta98Term!==g.term)return reject();
 }
 const registry={marker:m,byId};cache.set(atlas,registry);return registry;
}

export function liverTermParts(atlas:Atlas):Map<string,Part[]>{
 const r=liverControlRegistry(atlas);
 return new Map(r?.marker.termGroups.map(g=>[g.term,g.parts.map(id=>r.byId.get(id)!)])??[]);
}

export function isInstalledLiverPart(p:Part):boolean{
 const key=p.id.replace('TA98PART:ZA:Liver:','');
 return ['right','left','quadrate','caudate','capsule','serosa','subserosa','fibrous-appendix','papillary-process','caudate-process'].includes(key)&&!p.suppressed&&p.system==='digestive'
  &&p.conceptId===p.id&&p.sourceId==='Liver:'+key&&/^liver:(?:male|female):lobes-v1$/.test(p.tissueWallAssembly??'');
}
export const isInstalledLiverLobe=(p:Part)=>isInstalledLiverPart(p)&&!['capsule','serosa','subserosa'].includes(p.tissueWallLayer??'');

/** The capsule and the parenchyma are the two stored tissue stages. Gross
 * lobes share the same depth; peeling must never remove them one by one. */
export function liverTissuePeelSteps(atlas:Atlas,choice:HierarchyChoice){
 const r=liverControlRegistry(atlas);
 if(!r||!r.marker.partIds.every(id=>choice.elements.includes(id)))return null;
 return [
  ...['serosa','subserosa'].flatMap(layer=>{const ids=r.marker.partIds.filter(id=>r.byId.get(id)!.tissueWallLayer===layer);return ids.length?[{name:layer==='serosa'?'Tunica serosa':'Tela subserosa',ids}]:[];}),
  {name:'Fibrous capsule of liver',ids:r.marker.partIds.filter(id=>r.byId.get(id)!.tissueWallLayer==='capsule')},
  {name:'Liver parenchyma',ids:r.marker.partIds.filter(id=>isInstalledLiverLobe(r.byId.get(id)!))},
 ];
}
