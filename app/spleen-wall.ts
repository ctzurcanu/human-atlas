import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';

/** Only the explicitly installed finite tissues provide this peel sequence.
 * Pulp is a parenchymal cell, not an invented hollow-organ wall layer. */
export function spleenTissuePeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const marker=(atlas as Atlas&{finiteSpleenIntegration?:{schema:string;sex:string;partIds:string[];sourcePartId:string}}).finiteSpleenIntegration;
 if(marker?.schema!=='spleen-local-tissues-v203'||!['male','female'].includes(marker.sex)||choice.id==='selection-set'||!Array.isArray(marker.partIds)||marker.partIds.length!==4||new Set(marker.partIds).size!==4)return null;
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),parts=marker.partIds.map(id=>byId.get(id));
 if(parts.some(p=>!p||p.suppressed||!p.pairedTissue||!p.sectionSolid||p.tissueWallAssembly!=='splenic:'+marker.sex))return null;
 const cells=parts as Part[],first=cells[0],selected=new Set(choice.elements);
 if(selected.size!==4||!cells.every(p=>selected.has(p.id))||cells.some(p=>p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return null;
 const order=['serosa','capsule','pulp','trabeculae'];
 const steps=order.map(layer=>{const tissue=cells.filter(p=>p.tissueWallLayer===layer);return {name:tissue[0]?.name??layer,ids:tissue.map(p=>p.id)};});
 return steps.every(step=>step.ids.length===1)?steps:null;
}
