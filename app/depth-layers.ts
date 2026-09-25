import {structureName,type Part} from './anatomy';

// Keep this order from the outside of the body inward. Every source part has
// exactly one home, including annotations and structures the source cannot rank.
export const DEPTH_LAYERS=[
 {id:'skin',name:'Skin'},
 {id:'superficial-veins',name:'Superficial veins'},
 {id:'investing-fascia',name:'Investing fascia'},
 {id:'superficial-muscles',name:'Superficial muscles'},
 {id:'second-muscles',name:'Second muscle group'},
 {id:'intermediate-muscles',name:'Intermediate muscles'},
 {id:'deep-muscles',name:'Deep muscles'},
 {id:'deepest-muscles',name:'Deepest muscles'},
 {id:'visceral-coverings',name:'Visceral coverings'},
 {id:'anterior-organs',name:'Anterior organs'},
 {id:'deep-organs',name:'Deep organs'},
 {id:'deep-vessels',name:'Deep vessels'},
 {id:'lymphatic',name:'Lymph nodules & vessels'},
 {id:'deep-nerves',name:'Deep nerves'},
 {id:'ligaments',name:'Ligaments'},
 {id:'other',name:'Muscle attachments & other annotations'},
 {id:'thoracic-bones',name:'Thoracic cage bones, capsules & cartilages'},
 {id:'limb-bones',name:'Bones of the limbs'},
 {id:'other-bones',name:'Other bones'},
 {id:'spine',name:'Spine'},
 {id:'skull',name:'Skull'},
] as const;
export type DepthLayerId=typeof DEPTH_LAYERS[number]['id'];

const matches=(value:string,pattern:RegExp)=>pattern.test(value);
const depthCache=new WeakMap<Part,DepthLayerId>();
export function depthLayerFor(part:Part):DepthLayerId{
 const cached=depthCache.get(part);if(cached)return cached;
 const result=classifyDepthLayer(part);depthCache.set(part,result);return result;
}
function classifyDepthLayer(part:Part):DepthLayerId{
 const name=structureName(part.name).toLowerCase();
 const groups=(part.groups??[]).join(' ').toLowerCase();
 switch(part.system){
  case 'integumentary':case 'regions':case 'cell-boundary':return 'skin';
  case 'fascia':return matches(name,/pleura|peritone|pericardi|mesenter|omentum/)?'visceral-coverings':'investing-fascia';
  case 'venous':return matches(name,/\b(superficial|cephalic|basilic|saphenous|dorsal venous|cutaneous vein|median cubital)\b/)?'superficial-veins':'deep-vessels';
  case 'arterial':return 'deep-vessels';
  case 'lymphatic':return 'lymphatic';
  case 'nervous':return 'deep-nerves';
  case 'attachments':case 'schematic':return 'other';
  case 'muscular':{
   if(matches(name,/\b(bursa|ligament|tendon|tendinous|sheath|retinaculum|aponeurosis|iliotibial tract)\b/))return 'ligaments';
   if(matches(name,/\b(gluteus medius)\b/))return 'second-muscles';
   if(matches(name,/\b(gluteus minimus)\b/))return 'deep-muscles';
   if(matches(name,/\b(multifidus|rotatores|pronator quadratus|popliteus|obturator internus|obturator externus|gemellus|quadratus femoris|transversus abdominis|interspinal|intertransvers)\b/))return 'deepest-muscles';
   if(matches(groups,/superficial (?:gluteal )?muscles/)||matches(name,/\b(trapezius|latissimus dorsi|deltoid|pectoralis major|gluteus maximus|rectus abdominis|external oblique|sternocleidomastoid|biceps brachii|triceps brachii|gastrocnemius|sartorius|gracilis|platysma)\b/))return 'superficial-muscles';
   if(matches(name,/\b(pectoralis minor|gluteus medius|brachialis|soleus|semimembranosus|semitendinosus|biceps femoris|rectus femoris|vastus|serratus anterior|internal oblique|temporalis|masseter|brachioradialis)\b/))return 'second-muscles';
   if(matches(name,/\b(flexor digitorum superficialis|pronator teres|supinator|adductor longus|adductor brevis|fibularis|peroneus|levator scapulae|rhomboid|splenius|erector spinae|longissimus|iliocostalis)\b/))return 'intermediate-muscles';
   if(matches(groups,/deep gluteal muscles/)||matches(name,/\b(flexor digitorum profundus|subscapularis|supraspinatus|infraspinatus|teres minor|gluteus minimus|iliopsoas|psoas|quadratus lumborum|interossei|interosseous|adductor magnus)\b/))return 'deep-muscles';
   // The source lacks a depth tag for many small muscles. Keep these in the
   // middle rather than claiming a precise superficial or deepest plane.
   return 'intermediate-muscles';
  }
  case 'connective':{
   if(matches(name,/\b(capsul\w*|cartilag\w*|synovi\w*|menisc\w*|articular disc)\b/))return 'thoracic-bones';
   if(matches(name,/\b(pleura|peritone|pericardi|mesenter|omentum|tunica|serosa)\b/))return 'visceral-coverings';
   return 'ligaments';
  }
  case 'skeletal':{
   if(matches(name,/\b(skull|crani\w*|mandib\w*|maxill\w*|zygomat\w*|occipital|parietal|frontal bone|temporal bone|sphenoid|ethmoid|nasal concha bone|nasal bone|lacrimal bone|palatine bone|vomer|facial bone|orbital bone|canine|incisor|molar|premolar|tooth|teeth)\b/))return 'skull';
   if(matches(name,/\b(vertebr\w*|spine|spinal|sacrum|sacral|coccyx|lumbar|cervical vertebra|thoracic vertebra)\b/)||/^(atlas \(c1\)|axis \(c2\))/.test(name))return 'spine';
   if(matches(name,/\b(rib|sternum|costal|manubrium|xiphoid|thoracic cage|cartilage)\b/))return 'thoracic-bones';
   if(matches(name,/\b(femur|fibula|tibia|humerus|radius|ulna|clavicle|scapula|patella|calcaneus|talus|tarsal|metatars\w*|metacarp\w*|carpal|phalan\w*|navicular|cuboid bone|cuneiform|pisiform|trapezium|trapezoid|hamate|capitate|lunate|scaphoid|triquetrum|sesamoid bones of (?:hand|foot))\b/))return 'limb-bones';
   return 'other-bones';
  }
  case 'respiratory':return matches(name,/\b(pleura|membrane)\b/)?'visceral-coverings':matches(name,/\b(lung|trachea|larynx|bronch|nasal cavity)\b/)?'anterior-organs':'deep-organs';
  case 'digestive':return matches(name,/\b(omentum|mesenter|mesocolon|peritone|serosa)\b/)?'visceral-coverings':matches(name,/\b(stomach|liver|colon|intestin|jejunum|ileum|tongue|oesophagus|esophagus|pharynx|gallbladder)\b/)?'anterior-organs':'deep-organs';
  case 'cardiac':return 'anterior-organs';
  case 'urinary':case 'endocrine':return 'deep-organs';
  case 'reproductive':case 'pregnancy':return 'anterior-organs';
  case 'sensory':return matches(name,/\b(nerve|nuclei)\b/)?'deep-nerves':'deep-organs';
  default:return 'other';
 }
}
