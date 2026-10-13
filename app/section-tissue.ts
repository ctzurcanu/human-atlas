import type {Part} from './anatomy';
import {isSkinPart} from './depth-layers';

export type SectionTissue='integumentary'|'serosa'|'fascia'|'adipose'|'cartilage'|'tendon'|'lung'|'liver'|'spleen'|'cns'|'organ'|'gastric-wall'|'mucosa'|'submucosa'|'muscular'|'cardiac'|'skeletal'|'arterial'|'venous'|'nervous'|'lymphatic'|'connective'|'attachments';

export function sectionTissue(part:Part):SectionTissue{
 const name=part.name.toLowerCase(),material=(part.material??'').toLowerCase();
 if(isSkinPart(part))return 'integumentary';
 if(part.tissueWallLayer==='mucosal-folds')return 'mucosa';
 if(/\b(submucosa|submucosal|tela submucosa)\b/.test(name))return 'submucosa';
 if(/\b(muscularis mucosae|tunica muscularis|muscular coat|longitudinal layer|circular layer|stratum longitudinale|stratum circulare|fibrae obliquae|oblique fibres)\b/.test(name))return 'muscular';
 if(/\b(mucosa|mucosal)\b/.test(name))return 'mucosa';
 if(part.system==='lymphatic'&&!/\bspleen\b/.test(name))return 'lymphatic';
 if(/\b(pleura|peritone\w*|pericardi\w*|omentum|mesenter\w*|mesocolon|serosa)\b/.test(name))return 'serosa';
 if(part.system==='fascia'||/\b(fascia|aponeurosis|membrane)\b/.test(name))return 'fascia';
 if(/\b(liver|hepatic parenchyma)\b/.test(name))return 'liver';
 if(/\bspleen\b/.test(name))return 'spleen';
 if(part.system==='respiratory'&&/\b(lung|pulmonary lobe|alveol\w*)\b/.test(name))return 'lung';
 if(part.system==='nervous'&&/\b(brain|cerebr\w*|cerebell\w*|spinal cord|pons|medulla|thalam\w*|hypothalam\w*|cortex|white matter|gray matter|grey matter|hippocamp\w*)\b/.test(name))return 'cns';
 if(part.system==='integumentary')return 'adipose';
 if(part.system==='connective'&&/cartilage/.test(material+' '+name))return 'cartilage';
 if(part.system==='connective'&&/tendon/.test(material+' '+name))return 'tendon';
 if(part.system==='skeletal')return 'skeletal';
 if(part.system==='muscular')return 'muscular';
 if(part.system==='cardiac')return 'cardiac';
 if(part.system==='arterial'||part.system==='venous'||part.system==='nervous'||part.system==='connective'||part.system==='attachments')return part.system;
 if(part.system==='digestive'&&(/stomach|gastric/.test(part.sectionAssembly??'')||/\b(stomach|cardia|pyloric part)\b/.test(name)))return 'gastric-wall';
 return 'organ';
}
