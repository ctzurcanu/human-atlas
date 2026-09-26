import type {Part} from './anatomy';
import {isSkinPart} from './depth-layers';

// A single surface of a hollow organ does not describe solid tissue across its
// lumen. Draw a narrow cut edge for its wall instead of filling its outline.
export function hollowSectionWall(part:Part):number|undefined{
 const name=part.name.toLowerCase();
 if(part.system==='digestive'){
  if(/\b(oesophagus|esophagus|\w*pharynx|stomach|duodenum|jejunum|ileum|intestin\w*|colon|rectum|anal canal|appendix|gallbladder)\b/.test(name))return /stomach|colon|rectum|gallbladder/.test(name)?.004:.003;
  if(/\b(duct|gingiva|mucosa)\b/.test(name))return .0015;
 }
 if(part.system==='respiratory'&&/\b(bronch\w*|trachea|larynx)\b/.test(name))return .0025;
 if(part.system==='urinary'&&/\b(renal pelvis|ureter|urethra|urinary bladder)\b/.test(name))return .003;
 if(part.system==='reproductive'&&/\b(ductus deferens|ejaculatory duct)\b/.test(name))return .002;
 if(part.system==='cardiac'&&/^(left|right) (atrium|ventricle)$/.test(name))return .004;
 return undefined;
}

export function sectionCapProfile(part:Part){
 if(isSkinPart(part))return {thinShell:true,width:.0035,outermostOnly:part.name.startsWith('Body surface (derived)')};
 const wall=hollowSectionWall(part);
 if(wall)return {hollowWall:true,width:wall};
 if(part.system==='fascia'||part.system==='attachments'||/fascia|pleura|peritone|omentum|mesenter|mesocolon|serosa|capsule|membrane/i.test(part.name))return {thinShell:true,width:.0012};
 return {thinShell:false,width:part.system==='muscular'?.005:.002};
}
