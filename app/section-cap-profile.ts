import type {Part} from './anatomy';
import type {MeshTopology} from './mesh-topology';
import {isSkinPart} from './depth-layers';

// A single surface of a hollow organ does not describe solid tissue across its
// lumen. Draw a narrow cut edge for its wall instead of filling its outline.
export function hollowSectionWall(part:Part):number|undefined{
 const name=part.name.toLowerCase();
 const span=Math.min(...[0,1,2].map(axis=>part.bounds[1][axis]-part.bounds[0][axis]).filter(value=>value>0));
 const vesselWall=Math.max(.0007,Math.min(.004,span*.09));
 if(part.system==='arterial'||part.system==='venous')return vesselWall;
 if(/\b(?:oral|bucal|buccal|nasal|pharyngeal|laryngeal|paranasal) cavit(?:y|ies)\b|\b(?:maxillary|frontal|sphenoidal|ethmoidal) sinus\b/.test(name))return .002;
 if(part.system==='digestive'){
  // Digestive glands are parenchyma; coverings are thin sheets. All remaining
  // canal and cavity pieces expose a cut wall, never a solid lumen.
  if(/\b(liver|pancreas|salivary|parotid|submandibular|sublingual|tongue|tooth|teeth|spleen)\b/.test(name))return undefined;
  if(/\b(omentum|mesenter\w*|mesocolon|peritone\w*|serosa|fascia|ligament|mucosa|gingiva)\b/.test(name))return undefined;
  return /\bstomach\b/.test(name)?.007:/colon|rectum|gallbladder/.test(name)?.004:.003;
 }
 if(part.system==='respiratory'&&/\b(bronch\w*|trachea|larynx|nasal passage|nasopharynx)\b/.test(name))return .0025;
 if(part.system==='urinary'&&/\b(renal pelvis|ureter|urethra|urinary bladder)\b/.test(name))return .003;
 if(part.system==='reproductive'&&/\b(ductus deferens|ejaculatory duct|uterine cavity|vagina)\b/.test(name))return .002;
 if(part.system==='cardiac'&&/\b(heart|atri\w*|ventric\w*|auricle|cardiac chamber)\b/.test(name)&&!/\b(papillary muscle|myocardium|septum|valve|chorda)\b/.test(name))return .004;
 return undefined;
}

export function sectionCapProfile(part:Part,topology?:MeshTopology){
 if(isSkinPart(part))return {thinShell:true,width:.0035,outermostOnly:part.name.startsWith('Body surface (derived)')||part.id==='FJ2810'};
 const wall=hollowSectionWall(part);
 if(wall)return {hollowWall:true,width:wall};
 if(part.system==='fascia'||part.system==='attachments'||part.system!=='lymphatic'&&/fascia|pleura|peritone|omentum|mesentery|mesocolon|serosa|capsule|membrane/i.test(part.name))return {thinShell:true,width:.0012};
 // A mesh with many source boundary edges is a sheet, even if its anatomical
 // name sounds solid. Do not bridge its cut contour into invented volume.
 if(topology&&topology.boundaryRatio>.015)return {thinShell:true,width:part.system==='muscular'?.002:.0015};
 return {thinShell:false,width:part.system==='muscular'?.005:.002,closureWidth:topology?(topology.closed ? .002 : topology.boundaryRatio<.003 ? .012 : .004):undefined,closureFraction:topology?(topology.closed ? .04 : .08):undefined};
}
