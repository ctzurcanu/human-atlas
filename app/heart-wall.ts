import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';

export type CardiacTermGroup={key:string;name:string;partIds:string[];term:string;conceptId:string;parent:string|null;kind:'exact'|'parent'};
export type CardiacIntegration={schema:'cardiac-local-finite-v1';sex:'male'|'female';wallPartIds:string[];valvePartIds:string[];papillaryPartIds:string[];chordalPartIds?:string[];rootPartIds?:string[];wholeConceptId:string;termGroups?:CardiacTermGroup[];coupledBankSha256?:string};
export function cardiacIntegration(atlas:Atlas):CardiacIntegration|undefined{
 const marker=(atlas as Atlas&{finiteCardiacIntegration?:CardiacIntegration}).finiteCardiacIntegration;
 return marker?.schema==='cardiac-local-finite-v1'&&['male','female'].includes(marker.sex)&&Array.isArray(marker.wallPartIds)?marker:undefined;
}
/** Peel only the stored, coordinated finite wall; no camera or selection edits. */
export function cardiacWallPeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const marker=cardiacIntegration(atlas);if(!marker||choice.id==='selection-set'||!choice.elements.length)return null;
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),walls=marker.wallPartIds.map(id=>byId.get(id));
 if(walls.some(p=>!p||p.suppressed||!p.pairedTissue||!p.sectionSolid||p.normalSigns===undefined||p.tissueWallAssembly!==`heart:${marker.sex}`))return null;
 const wall=walls as Part[],first=wall[0];if(!first||wall.some(p=>p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)||JSON.stringify(p.geometryInstance)!==JSON.stringify(first.geometryInstance)))return null;
 const allowed=new Set([...marker.wallPartIds,...marker.valvePartIds,...marker.papillaryPartIds,...marker.chordalPartIds??[],...marker.rootPartIds??[]]);
 if(choice.elements.some(id=>!allowed.has(id)||!byId.get(id)||byId.get(id)!.suppressed))return null;
 const selected=new Set(choice.elements),parts=wall.filter(p=>selected.has(p.id));
 const steps:WallPeelStep[]=[
  {name:'Epicardium',ids:parts.filter(p=>p.tissueWallLayer==='epicardium').map(p=>p.id)},
  {name:'Myocardium and septal tissue',ids:parts.filter(p=>['myocardium','fibrous-septum'].includes(p.tissueWallLayer??'')).map(p=>p.id)},
  {name:'Endocardium',ids:parts.filter(p=>p.tissueWallLayer==='endocardium').map(p=>p.id)},
 ];
 if(steps.some(step=>!step.ids.length))return null;
 const apparatus=[...marker.valvePartIds,...marker.papillaryPartIds,...marker.chordalPartIds??[]].filter(id=>selected.has(id));
 if(apparatus.length)steps.push({name:marker.chordalPartIds?.length?'Valves, chordae and papillary muscles':'Valves and papillary muscles',ids:apparatus});
 const roots=(marker.rootPartIds??[]).filter(id=>selected.has(id));
 if(roots.length)steps.push({name:'Vessel roots',ids:roots});
 return steps.flatMap(step=>step.ids).length===selected.size?steps:null;
}

/** Catalogue-declared controls reference stored tissues, never extra meshes. */
export function cardiacTermParts(atlas:Atlas):Map<string,Part[]>{
 const marker=cardiacIntegration(atlas),result=new Map<string,Part[]>();if(!marker?.termGroups)return result;
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 const allowed=new Set([...marker.wallPartIds,...marker.valvePartIds,...marker.papillaryPartIds,...marker.chordalPartIds??[],...marker.rootPartIds??[]]);
 for(const g of marker.termGroups){
  const c=concepts.get(g.conceptId);if(g.kind!=='exact'||!c||c.ta98Term!==g.term||c.ta98Kind!=='exact'||!g.partIds.length||g.partIds.some(id=>!allowed.has(id)||!byId.has(id)||byId.get(id)!.suppressed)||c.elements.length!==g.partIds.length||c.elements.some(id=>!g.partIds.includes(id)))continue;
  result.set(g.term,[...new Map([...(result.get(g.term)??[]),...g.partIds.map(id=>byId.get(id)!)].map(p=>[p.id,p])).values()]);
 }
 return result;
}
