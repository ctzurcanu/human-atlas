import {structureName,type Part,type Concept} from './anatomy';
import {taExactEntityForConcept} from './anatomical-terminology';
import skeletalAttachments from './generated-skeletal-attachments.json';

export type RelationKind='before'|'after'|'innervation'|'arterial'|'venous'|'innervates'|'supplies'|'drains'|'articulates'|'connects'|'connectedBy'|'joint'|'continuous'|'covers'|'coveredBy'|'cartilages'|'bones'|'tendons'|'fascia'|'muscles'|'origin'|'insertion'|'partOf'|'contains'|'originFor'|'insertionFor'|'originSites'|'insertionSites'|'counterpart';
export interface ResolvedRelation {kind:RelationKind;target:Part;via?:string[];viaModeled?:boolean;note?:string}
type Relations=Partial<Record<RelationKind,string[]>>;

// These are anatomical paths, not spatial guesses. Names are resolved against
// leaves in the active atlas; omitted intermediate structures are shown as text.
// Sources: OpenStax A&P 2e, ch. 20/23; NCBI Bookshelf NBK482334,
// NBK448070, NBK556044, NBK470256, NBK470197, NBK482390.
const digestiveRoute=[
 ['Oropharynx'],['Laryngopharynx'],['Oesophagus','Esophagus'],['Stomach'],
 ['Duodenum','Superior part of duodenum'],['Jejunum','Proximal part of jejunum'],['Middle part of jejunum'],
 ['Distal part of jejunum'],['Ileum','Proximal part of ileum'],['Middle part of ileum'],['Distal part of ileum'],
 ['Caecum','Cecum'],['Ascending colon'],['Hepatic flexure of colon','Right colic flexure'],
 ['Transverse colon'],['Splenic flexure of colon','Left colic flexure'],
 ['Descending colon'],['Sigmoid colon'],['Rectum'],['Anal canal'],
];
const airwayRoute=[
 ['Nasopharynx'],['Oropharynx'],['Laryngopharynx'],['Larynx'],['Trachea'],
 ['Right main bronchus','Right main bronchus proper'],
];
const leftAirwayRoute=[['Trachea'],['Left main bronchus']];
const arterialRoute=[
 ['Ascending aorta'],['Aortic arch','Arch of aorta'],['Thoracic aorta'],
 ['Abdominal aorta'],['Common iliac artery (left)','Common iliac artery (right)'],
];
const legArterialRoute=[
 ['Common iliac artery (left)','Common iliac artery (right)'],
 ['External iliac artery (left)','External iliac artery (right)'],
 ['Femoral artery (left)','Femoral artery (right)'],
 ['Popliteal artery (left)','Popliteal artery (right)'],
];
const legVenousRoute=[
 ['Popliteal vein (left)','Popliteal vein (right)'],
 ['Femoral vein (left)','Femoral vein (right)'],
 ['External iliac vein (left)','External iliac vein (right)'],
 ['Common iliac vein (left)','Common iliac vein (right)'],
 ['Inferior vena cava (abdominal part)','Inferior vena cava'],
];
const armArterialRoute=[
 ['Subclavian artery (left)','Subclavian artery (right)'],
 ['Axillary artery (left)','Axillary artery (right)'],
 ['Brachial artery (left)','Brachial artery (right)'],
];
const armVenousRoute=[
 ['Radial veins (left)','Radial veins (right)','Ulnar veins (left)','Ulnar veins (right)'],
 ['Brachial veins (left)','Brachial veins (right)'],
 ['Axillary vein (left)','Axillary vein (right)'],
 ['Subclavian vein (left)','Subclavian vein (right)'],
];
const portalVenousRoute=[
 ['Inferior mesenteric vein'],['Splenic vein'],['Hepatic portal vein'],
];
const portalMesentericRoute=[['Superior mesenteric vein'],['Hepatic portal vein']];
const greatSaphenousRoute=[['Great saphenous vein (left)','Great saphenous vein (right)'],['Femoral vein (left)','Femoral vein (right)']];
const smallSaphenousRoute=[['Small saphenous vein (left)','Small saphenous vein (right)'],['Popliteal vein (left)','Popliteal vein (right)']];
const rightHeartRoute=[['Superior vena cava'],['Right atrium'],['Right ventricle'],['Pulmonary trunk']];
const inferiorCavaRoute=[['Inferior vena cava (abdominal part)'],['Inferior vena cava (thoracic part)'],['Right atrium']];
const leftHeartRoute=[
 ['Left superior pulmonary vein','Left inferior pulmonary vein','Right superior pulmonary vein','Right inferior pulmonary vein'],
 ['Left atrium'],['Left ventricle'],['Ascending aorta'],
];
// Urine flows through the renal pelvis and ureter to the bladder; sperm passes
// through the epididymis and ductus deferens before entering the urethra.
// Sources: OpenStax A&P 2e 25.2 and 27.1; NCBI Bookshelf NBK562291.
const urinaryRoute=[
 ['Kidney'],['Renal pelvis'],['Ureter'],['Urinary bladder'],['Urethra'],
];
const spermRoute=[
 ['Testis'],['Epididymis'],['Ductus deferens','Vas deferens'],
 ['Ejaculatory duct'],['Prostatic urethra','Urethra'],
];
const airwayBranches:[string[],string[][]][]=[
 [['Trachea'],[['Right main bronchus','Right main bronchus proper'],['Left main bronchus']]],
 [['Right main bronchus','Right main bronchus proper'],[['Right superior lobar bronchus'],['Intermediate bronchus (right)']]],
 [['Intermediate bronchus (right)'],[['Middle lobar bronchus (right)'],['Right inferior lobar bronchus']]],
 [['Left main bronchus'],[['Left superior lobar bronchus'],['Left inferior lobar bronchus']]],
 [['Right superior lobar bronchus'],[['Apical segmental bronchus of right lung (BI)'],['Posterior segmental bronchus of right lung (BII)'],['Anterior segmental bronchus of right lung (BIII)']]],
 [['Middle lobar bronchus (right)'],[['Lateral segmental bronchus of right lung (BIV)'],['Medial segmental bronchus of right lung (BV)']]],
 [['Right inferior lobar bronchus'],[['Superior segmental bronchus of right lung (BVI)'],['Medial basal segmental bronchus of right lung (BVII)'],['Anterior basal segmental bronchus of right lung (BVIII)'],['Lateral basal segmental bronchus of right lung (BIX)'],['Posterior basal segmental bronchus of right lung (BX)']]],
 [['Left superior lobar bronchus'],[['Apicoposterior segmental bronchus of left lung (BI+BII)'],['Anterior segmental bronchus of left lung (BIII)'],['Superior lingular segmental bronchus of left lung (BIV)'],['Inferior lingular segmental bronchus of left lung (BV)']]],
 [['Left inferior lobar bronchus'],[['Superior segmental bronchus of left lung (BVI)'],['Medial basal segmental bronchus of left lung (BVII)'],['Anterior basal segmental bronchus of left lung (BVIII)','(Anteromedial basal segmental bronchus of left lung)'],['Lateral basal segmental bronchus of left lung (BIX)'],['Posterior basal segmental bronchus of left lung (BX)']]],
];
const arterialBranches:[string[],string[][]][]=[
 [['Pulmonary trunk'],[['Right pulmonary artery'],['Left pulmonary artery']]],
 [['Coeliac trunk','Celiac trunk'],[['Left gastric artery'],['Common hepatic artery'],['Splenic artery']]],
 [['Common hepatic artery'],[['Proper hepatic artery'],['Gastroduodenal artery']]],
 [['Superior mesenteric artery'],[['Ileocolic artery'],['Right colic artery'],['Middle colic artery']]],
 [['Inferior mesenteric artery'],[['Left colic artery'],['Sigmoid arteries'],['Superior rectal artery','Superior anorectal artery']]],
 [['Common iliac artery (left)','Common iliac artery (right)'],[['External iliac artery (left)','External iliac artery (right)'],['Internal iliac artery (left)','Internal iliac artery (right)']]],
 [['Femoral artery (left)','Femoral artery (right)'],[['Deep femoral artery (left)','Deep femoral artery (right)'],['Popliteal artery (left)','Popliteal artery (right)']]],
 [['Popliteal artery (left)','Popliteal artery (right)'],[['Anterior tibial artery (left)','Anterior tibial artery (right)'],['Posterior tibial artery (left)','Posterior tibial artery (right)']]],
 [['Brachial artery (left)','Brachial artery (right)'],[['Radial artery (left)','Radial artery (right)'],['Ulnar artery (left)','Ulnar artery (right)']]],
];
const nerveBranches:[string[],string[][]][]=[
 [['Lateral cord of brachial plexus'],[['Lateral root of median nerve'],['Musculocutaneous nerve'],['Lateral pectoral nerve']]],
 [['Medial cord of brachial plexus'],[['Medial root of median nerve'],['Ulnar nerve'],['Medial pectoral nerve'],['Medial antebrachial cutaneous nerve']]],
 [['Posterior cord of brachial plexus'],[['Radial nerve'],['Axillary nerve'],['Thoracodorsal nerve']]],
 [['Lateral root of median nerve'],[['Median nerve']]],
 [['Medial root of median nerve'],[['Median nerve']]],
 [['Sciatic nerve'],[['Tibial nerve'],['Common fibular nerve']]],
 [['Common fibular nerve'],[['Deep fibular nerve'],['Superficial fibular nerve'],['Sural communicating branch of common fibular nerve']]],
 [['Tibial nerve'],[['Medial plantar nerve'],['Lateral plantar nerve']]],
 [['Radial nerve'],[['Deep branch of radial nerve','Radial nerve (deep branch)'],['Superficial branch of radial nerve','Radial nerve (superficial br)']]],
 [['Ulnar nerve'],[['Deep branch of ulnar nerve'],['Superficial branch of ulnar nerve'],['Dorsal branch of ulnar nerve'],['Palmar branch of ulnar nerve']]],
 [['Median nerve'],[['Anterior interosseous nerve of forearm'],['Palmar branch of median nerve'],['Common palmar digital branches of median nerve']]],
 [['Femoral nerve'],[['Saphenous nerve'],['Anterior cutaneous branches of femoral nerve']]],
 [['Obturator nerve'],[['Anterior branch of obturator nerve'],['Posterior branch of obturator nerve']]],
 [['Trigeminal nerve (V)'],[['Ophthalmic nerve'],['Maxillary nerve'],['Anterior division of mandibular nerve'],['Posterior division of mandibular nerve']]],
];

// Joint labels come from the atlas's ligament/capsule groups. Each mapping is
// restricted to the bones forming that named joint; digit/rib groups without a
// specific bone number are intentionally excluded.
const jointBones:Record<string,string[]>={
 'acromioclavicular joint':['Clavicle','Scapula'],
 'sternoclavicular joint':['Clavicle','Manubrium of sternum'],
 'glenohumeral joint':['Scapula','Humerus'],
 'elbow joint':['Humerus','Ulna','Radius'],
 'distal radio-ulnar joint':['Radius','Ulna'],
 'radiocarpal joint':['Radius','Scaphoid bone','Lunate bone'],
 'hip joint':['Hip bone','Femur'],
 'knee joint':['Femur','Tibia','Patella'],
 'superior tibiofibular joint':['Tibia','Fibula'],
 'sacro-iliac joint':['Sacrum','Hip bone'],
 'temporomandibular joint':['Temporal bone','Mandible'],
 'subtalar joint':['Talus','Calcaneus'],
 'calcaneocuboid joint':['Calcaneus','Cuboid bone'],
 'cuboidonavicular joint':['Cuboid bone','Navicular bone'],
 'cuneocuboid joint':['Lateral cuneiform bone','Cuboid bone'],
 'cuneonavicular joint':['Navicular bone','Medial cuneiform bone','Intermediate cuneiform bone','Lateral cuneiform bone'],
 'talocalcaneonavicular joint':['Talus','Calcaneus','Navicular bone'],
 'pisiform joint':['Pisiform bone','Triquetrum bone'],
};
const fascialContinuity=[
 ['Deltoid fascia','Brachial fascia'],
 ['Brachial fascia','Antebrachial fascia'],
 ['Brachial fascia','Lateral intermuscular septum of arm'],
 ['Brachial fascia','Medial intermuscular septum of arm'],
 ['Antebrachial fascia','Palmar aponeurosis'],
 ['Antebrachial fascia','Dorsal fascia of hand'],
 ['Fascia lata','Crural fascia'],
 ['Fascia lata','Medial femoral intermuscular septum'],
 ['Fascia lata','Lateral femoral intermuscular septum'],
 ['Crural fascia','Anterior intermuscular septum of leg'],
 ['Crural fascia','Posterior intermuscular septum of leg'],
 ['Crural fascia','Transverse intermuscular septum of leg'],
 ['Crural fascia','Popliteal fascia'],
 ['Anterior layer of thoracolumbar fascia','Middle layer of thoracolumbar fascia'],
 ['Middle layer of thoracolumbar fascia','Posterior layer of thoracolumbar fascia'],
];
const fascialCoverings:Record<string,string[]>={
 'deltoid fascia':['Deltoid muscle'],
 'piriformis fascia':['Piriformis muscle'],
 'masseteric fascia':['Deep part of masseter','Superficial part of masseter'],
 'iliopsoas fascia':['Iliopsoas muscle'],
 'pectoral fascia':['Pectoralis major'],
 'clavipectoral fascia':['Pectoralis minor muscle','Subclavius muscle'],
 'diaphragmatic fascia':['Diaphragm'],
 'transversalis fascia':['Transversus abdominis muscle'],
 'superficial layer of temporal fascia':['Temporalis muscle'],
};
const namedTendons:Record<string,string[]>={
 'calcaneal tendon':['Calcaneus','Lateral head of gastrocnemius','Medial head of gastrocnemius','Soleus muscle'],
 'common tendon of biceps brachii':['Radius','Biceps brachii muscle'],
 'common tendon of triceps brachii':['Ulna','Triceps brachii muscle'],
 // The quadriceps tendon inserts on the superior patella; the source model
 // has no distal .e marker for rectus femoris or the vasti.
 'rectus femoris muscle':['Patella'],
 'vastus lateralis muscle':['Patella'],
 'vastus medialis muscle':['Patella'],
 'vastus intermedius muscle':['Patella'],
 // Distal slips vary, so these are the documented attachment sites rather
 // than a claim that every individual has every slip.
 'tibialis posterior muscle':['Navicular bone','Medial cuneiform bone','Intermediate cuneiform bone','Lateral cuneiform bone','Cuboid bone','Second metatarsal bone','Third metatarsal bone','Fourth metatarsal bone'],
 'ascending part of trapezius muscle':['Scapula'],
 'transverse part of trapezius muscle':['Scapula'],
};
// NCBI Bookshelf NBK519516: palmaris longus ends in aponeurosis/retinaculum,
// not a distal bone. Show that interface in both directions.
const tendonFascialAttachments:Record<string,string[]>={
 'palmaris longus muscle':['Palmar aponeurosis'],
};
// Named fascial attachments are added only where the attachment is described
// anatomically. Spatial overlap alone does not establish fascial continuity.
// Sources: NCBI Bookshelf NBK557497 (fascia lata), NBK518984 (thoracolumbar);
// PMC10625514 (crural), PMC10426691 (antebrachial).
const fascialBoneAttachments:Record<string,string[]>={
 'fascia lata':['Hip bone','Sacrum','Femur','Tibia','Fibula'],
 'iliotibial tract':['Hip bone','Tibia'],
 'crural fascia':['Tibia','Fibula'],
 'antebrachial fascia':['Humerus','Ulna'],
 'anterior layer of thoracolumbar fascia':['Hip bone','Sacrum','Vertebra L5'],
 'middle layer of thoracolumbar fascia':['Hip bone','Sacrum','Vertebra L5'],
 'posterior layer of thoracolumbar fascia':['Hip bone','Sacrum','Vertebra L5'],
};
// Bone blood and periosteal nerve supply. Do not infer these from mere spatial
// proximity to a vessel/nerve. Sources: NCBI Bookshelf NBK532982 (femur),
// NBK526053 (tibia), NBK470591 (fibula), NBK534821 (humerus), NBK549847 (rib),
// NBK519534 (patella), NBK538319 (scapula), NBK545260 (forearm),
// NBK537024 (leg veins), NBK541141 (sternum), NBK525990 (clavicle),
// NBK551653 (sacrum), NBK519524 (hip), NBK532292 (mandible),
// NBK538527 (maxilla), NBK544319 (foot), NBK541086 (talus), NBK557549
// (frontal periosteum), NBK541093 (occipital), NBK482497 (temporal),
// NBK544257 (zygomatic), NBK578171 (foot sesamoids), PMC5640352 (parietal).
const boneSupply:Record<string,Relations>={
 'femur':{arterial:['Medial circumflex femoral artery','Lateral circumflex femoral artery','Deep femoral artery'],innervation:['Femoral nerve','Obturator nerve','Tibial nerve','Common fibular nerve']},
 'tibia':{arterial:['Posterior tibial artery','Anterior tibial artery'],venous:['Posterior tibial veins','Anterior tibial veins'],innervation:['Tibial nerve','Deep fibular nerve']},
 'fibula':{arterial:['Fibular artery','Anterior tibial artery'],venous:['Fibular veins']},
 'radius':{arterial:['Anterior interosseous artery','Radial artery']},
 'ulna':{arterial:['Ulnar artery','Anterior interosseous artery']},
 'humerus':{arterial:['Anterior circumflex humeral artery','Posterior circumflex humeral artery','Brachial artery','Deep brachial artery']},
 'scapula':{arterial:['Suprascapular artery','Circumflex scapular artery','Dorsal scapular artery']},
 'patella':{arterial:['Superior lateral genicular artery','Superior medial genicular artery','Inferior lateral genicular artery','Inferior medial genicular artery']},
 'clavicle':{arterial:['Suprascapular artery','Thoracoacromial artery','Internal thoracic artery']},
 'body of sternum':{arterial:['Internal thoracic artery'],venous:['Internal thoracic veins']},
 'manubrium of sternum':{arterial:['Internal thoracic artery'],venous:['Internal thoracic veins']},
 'xiphoid process':{arterial:['Internal thoracic artery'],venous:['Internal thoracic veins']},
 'sacrum':{arterial:['Median sacral artery']},
 'hip bone':{arterial:['Obturator artery','Iliolumbar artery','Superior gluteal artery']},
 'mandible':{arterial:['Inferior alveolar artery']},
 'maxilla':{arterial:['Maxillary artery','Infraorbital artery','Posterior superior alveolar artery']},
 'frontal bone':{arterial:['Supratrochlear artery']},
 'parietal bone':{arterial:['Superficial temporal artery']},
 'occipital bone':{arterial:['Occipital artery']},
 'temporal bone':{arterial:['Ascending pharyngeal artery']},
 'zygomatic bone':{arterial:['Maxillary artery','Facial artery']},
 'sesamoid bones of foot':{arterial:['Medial plantar artery','Plantar metatarsal arteries']},
 'talus':{arterial:['Posterior tibial artery','Fibular artery','Dorsalis pedis artery'],innervation:['Tibial nerve','Saphenous nerve','Deep fibular nerve']},
 'calcaneus':{arterial:['Calcaneal branches of posterior tibial artery','Calcaneal branches of fibular artery']},
 'navicular bone':{arterial:['Dorsalis pedis artery','Posterior tibial artery']},
 'cuboid bone':{arterial:['Lateral plantar artery']},
 'medial cuneiform bone':{arterial:['Medial plantar artery','Dorsalis pedis artery']},
 'intermediate cuneiform bone':{arterial:['Medial plantar artery','Dorsalis pedis artery']},
 'lateral cuneiform bone':{arterial:['Medial plantar artery','Dorsalis pedis artery']},
};
// Regional bone supply where the source describes a group of bones rather
// than an individual named bone. Do not extend this to cartilage or nerves.
// Sources: NCBI Bookshelf NBK525969, NBK459153 (vertebrae), NBK535382
// (carpals), NBK544319 (foot), NBK549872 (metatarsals), NBK538428 (hand).
function regionalBoneSupply(name:string):Relations|undefined{
 if(/^(?:atlas \(c1\)|axis \(c2\)|vertebra c\d+)$/.test(name))return {arterial:['Vertebral artery','Deep cervical artery']};
 if(/^vertebra t\d+$/.test(name))return {arterial:['Posterior intercostal arteries']};
 if(/^vertebra l\d+$/.test(name))return {arterial:['Lumbar arteries','Iliolumbar artery']};
 if(/^(?:capitate|hamate|lunate|pisiform|scaphoid|trapezium|trapezoid|triquetrum) bone$/.test(name))return {arterial:['Radial artery','Ulnar artery']};
 if(/^(?:first|second|third|fourth|fifth) metacarpal bone$/.test(name))return {arterial:['Palmar metacarpal arteries','Dorsal metacarpal arteries']};
 if(/^\w+ phalanx of \w+ finger of hand$/.test(name))return {arterial:['Proper palmar digital arteries','Dorsal digital arteries of hand']};
 if(/^(?:first|second|third|fourth|fifth) metatarsal bone$/.test(name))return {arterial:['Dorsal metatarsal arteries','Plantar metatarsal arteries']};
 if(/^\w+ phalanx of \w+ finger of foot$/.test(name))return {arterial:['Proper plantar digital arteries','Dorsal digital arteries of foot']};
 return undefined;
}
// Named muscular arteries are associated with modeled muscle leaves only.
// Sources: NCBI Bookshelf NBK537148, NBK537202, NBK513255, NBK538262,
// NBK537145, NBK539705, NBK539913.
const muscleArterialSupply:Record<string,string[]>={
 'deltoid muscle':['Posterior circumflex humeral artery','Thoracoacromial artery Deltoid br'],
 'supraspinatus muscle':['Suprascapular artery','Dorsal scapular artery'],
 'infraspinatus muscle':['Suprascapular artery','Circumflex scapular artery'],
 'transverse part of trapezius muscle':['Transverse cervical artery'],
 'ascending part of trapezius muscle':['Transverse cervical artery'],
 'descending part of trapezius muscle':['Transverse cervical artery'],
 'serratus anterior muscle':['Circumflex scapular artery'],
 'rectus femoris muscle':['Lateral circumflex femoral artery'],
 'vastus lateralis muscle':['Lateral circumflex femoral artery'],
 'vastus medialis muscle':['Femoral artery','Deep femoral artery'],
 'biceps brachii muscle':['Brachial artery'],
 'triceps brachii muscle':['Brachial artery','Deep brachial artery'],
 'gastrocnemius muscle':['Popliteal artery'],
 'soleus muscle':['Posterior tibial artery','Fibular artery'],
 'tibialis posterior muscle':['Posterior tibial artery'],
};
const ligamentAttachments:Record<string,string[]>={
 'anterior cruciate ligament':['Femur','Tibia'],
 'posterior cruciate ligament':['Femur','Tibia'],
 'fibular collateral ligament':['Femur','Fibula'],
 'tibial collateral ligament':['Femur','Tibia'],
 'patellar ligament':['Patella','Tibia'],
 'acromioclavicular ligament':['Clavicle','Scapula'],
 'coracohumeral ligament':['Scapula','Humerus'],
 'iliofemoral ligament':['Hip bone','Femur'],
 'ischiofemoral ligament':['Hip bone','Femur'],
 'pubofemoral ligament':['Hip bone','Femur'],
 'ligament of head of femur':['Hip bone','Femur'],
 'anterior talofibular ligament':['Talus','Fibula'],
 'calcaneofibular ligament':['Calcaneus','Fibula'],
};

// Whole-organ pelvic links, never inherited by unnamed wall/layer regions.
// NCBI Bookshelf NBK554601, NBK547660, NBK545187, NBK557575;
// PMC3312145 distinguishes upper vaginal autonomic and lower somatic supply.
// Targets must exist in the active model; these links do not certify geometry.
const uterineFunctional:Relations={innervation:['Uterovaginal plexus'],arterial:['Uterine artery','Ovarian artery'],venous:['Uterine vein','Uterine venous plexus']};
const cervicalFunctional:Relations={innervation:['Uterovaginal plexus'],arterial:['Uterine artery','Vaginal artery'],venous:['Uterine vein']};
const tubalFunctional:Relations={arterial:['Uterine artery','Ovarian artery'],venous:['Uterine vein','Ovarian vein']};
const femalePelvicFunctional:Record<string,Relations>={
 'uterus':uterineFunctional,
 'cervix':cervicalFunctional,
 'cervix of uterus':cervicalFunctional,
 'uterine cervix':cervicalFunctional,
 'vagina':{innervation:['Uterovaginal plexus','Pudendal nerve'],arterial:['Uterine artery','Vaginal artery','Internal pudendal artery'],venous:['Vaginal vein','Vaginal venous plexus']},
 'uterine tube':tubalFunctional,
 'fallopian tube':tubalFunctional,
 'ovary':{arterial:['Ovarian artery','Uterine artery'],venous:['Ovarian vein']},
};
const explicit:Record<string,Relations>={
 ...femalePelvicFunctional,
 'oesophagus':{innervation:['Vagus nerve (X)'],arterial:['Left gastric artery'],venous:['Left gastric vein']},
 'esophagus':{innervation:['Vagus nerve (X)'],arterial:['Left gastric artery'],venous:['Left gastric vein']},
 'stomach':{innervation:['Vagus nerve (X)','Vagus nerve'],arterial:['Left gastric artery','Right gastric artery','Left gastro-omental artery','Right gastro-omental artery','Left gastroepiploic artery','Right gastroepiploic artery','Gastroepiploic artery','Short gastric arteries','Splenic artery'],venous:['Left gastric vein','Right gastric vein','Left gastro-omental vein','Right gastro-omental vein','Left gastroepiploic vein','Right gastroepiploic vein','Splenic vein']},
 'duodenum':{innervation:['Vagus nerve (X)'],arterial:['Superior pancreaticoduodenal artery','Inferior pancreaticoduodenal artery'],venous:['Superior mesenteric vein','Hepatic portal vein']},
 'jejunum':{innervation:['Vagus nerve (X)'],arterial:['Superior mesenteric artery'],venous:['Superior mesenteric vein']},
 'ileum':{innervation:['Vagus nerve (X)'],arterial:['Superior mesenteric artery'],venous:['Superior mesenteric vein']},
 'ascending colon':{innervation:['Vagus nerve (X)'],arterial:['Ileocolic artery','Right colic artery','Superior mesenteric artery'],venous:['Right colic vein','Superior mesenteric vein']},
 'transverse colon':{arterial:['Middle colic artery'],venous:['Middle colic vein']},
 'descending colon':{arterial:['Left colic artery','Inferior mesenteric artery'],venous:['Inferior mesenteric vein']},
 'sigmoid colon':{arterial:['Sigmoid arteries','Inferior mesenteric artery'],venous:['Inferior mesenteric vein']},
 'rectum':{arterial:['Superior rectal artery','Superior anorectal artery'],venous:['Superior rectal vein','Inferior mesenteric vein']},
 'trachea':{innervation:['Vagus nerve (X)','Recurrent laryngeal nerve'],arterial:['Inferior thyroid artery','Bronchial artery'],venous:['Inferior thyroid vein','Bronchial vein']},
 'right main bronchus':{innervation:['Vagus nerve (X)'],arterial:['Bronchial artery'],venous:['Bronchial vein']},
 'left main bronchus':{innervation:['Vagus nerve (X)'],arterial:['Bronchial artery'],venous:['Bronchial vein']},
 'heart':{innervation:['Vagus nerve (X)','Sympathetic nerves'],arterial:['Left coronary artery','Right coronary artery'],venous:['Coronary sinus','Great cardiac vein']},
 'spleen':{arterial:['Splenic artery'],venous:['Splenic vein']},
 'liver':{innervation:['Vagus nerve (X)'],arterial:['Proper hepatic artery','Hepatic artery'],venous:['Hepatic veins']},
 'pancreas':{innervation:['Vagus nerve (X)'],arterial:['Splenic artery','Superior mesenteric artery'],venous:['Splenic vein','Superior mesenteric vein']},
 'kidney':{arterial:['{side} renal artery'],venous:['{side} renal vein']},
 'suprarenal gland':{arterial:['Superior suprarenal artery','Middle suprarenal artery','Inferior suprarenal artery'],venous:['Suprarenal vein']},
 'transverse part of trapezius muscle':{innervation:['Accessory nerve (XI)']},
 'ascending part of trapezius muscle':{innervation:['Accessory nerve (XI)']},
 'descending part of trapezius muscle':{innervation:['Accessory nerve (XI)']},
};

const clean=(name:string)=>structureName(name).replace(/\.[eo]\d*[lr]$/i,'').replace(/^\((.+)\)(\s*\((?:left|right)\))$/i,'$1$2').replace(/^\((.+)\)$/,'$1').replace(/\s+/g,' ').trim();
const key=(name:string)=>clean(name).toLowerCase();
const sideOf=(name:string):'left'|'right'|null=>{
 if(/^(?:left|right) (?:atrium|ventricle)$/i.test(clean(name)))return null;
 const suffix=structureName(name).match(/\.[eo]\d*([lr])$/i);
 if(suffix)return suffix[1].toLowerCase()==='l'?'left':'right';
 const m=clean(name).match(/\((left|right)\)$/i)??clean(name).match(/^(left|right)\s/i);
 if(m)return m[1].toLowerCase() as 'left'|'right';
 const embedded=clean(name).match(/\bof (left|right)\s/i);
 if(embedded)return embedded[1].toLowerCase() as 'left'|'right';
 const abbreviated=clean(name).match(/\s([LR])$/i);
 return abbreviated?abbreviated[1].toLowerCase()==='l'?'left':'right':null;
};
const partSide=(part:Part):'left'|'right'|null=>{
 const named=sideOf(part.name);
 if(named||!part.id.startsWith('HRA:'))return named;
 const source=part.sourceId??part.id;
 const written=source.match(/(?:^|_)(left|right)(?:_|$)/i)?.[1];
 if(written)return written.toLowerCase() as 'left'|'right';
 const suffix=source.match(/_([LR])(?:_|$)/i)?.[1];
 return suffix?suffix.toLowerCase()==='l'?'left':'right':null;
};
const withoutSide=(name:string)=>key(name).replace(/\s*\((?:left|right)\)$/,'').replace(/\s+[lr]$/,'');
const intestinalNames:Record<string,string>={'intestine duodenum':'duodenum','small intestine jejunum':'jejunum','small intestine ileum':'ileum','small intestine illium':'ileum','large intestine cecum':'cecum','large intestine descending colon':'descending colon','large intestine rectum':'rectum'};
const baseName=(name:string)=>{const value=/^(?:left|right) (?:atrium|ventricle)$/.test(withoutSide(name))?withoutSide(name):withoutSide(name).replace(/^(?:left|right)\s+/,'');return intestinalNames[value]??value;};
const belongsTo=(name:string,base:string)=>name===base||name.startsWith(`proximal part of ${base}`)||name.startsWith(`middle part of ${base}`)||name.startsWith(`distal part of ${base}`);
const vascularIdentity=(name:string)=>{
 const side=sideOf(name);
 return `${withoutSide(name).replace(/^(?:left|right)\s+/,'')}|${side??''}`;
};
const statedVascularParent=(name:string)=>clean(name).match(/\b(?:branch(?:es)?|tributar(?:y|ies)) of (.+)$/i)?.[1]??null;

const gastricRegionCodes=new Set(['A05.5.01.007','A05.5.01.009','A05.5.01.012','A05.5.01.014']);
const gastricAssemblyCache=new WeakMap<Part[],Part[]>();
function modeledGastricAssembly(parts:Part[]):Part[]{
 const cached=gastricAssemblyCache.get(parts);if(cached)return cached;
 const assemblies=new Map<string,Part[]>();
 for(const part of parts){
  const code=taExactEntityForConcept(part.conceptId);
  if(part.suppressed||part.system!=='digestive'||!part.sectionAssembly||!code||!gastricRegionCodes.has(code))continue;
  const group=assemblies.get(part.sectionAssembly)??[];group.push(part);assemblies.set(part.sectionAssembly,group);
 }
 const regions=[...assemblies.values()].find(group=>new Set(group.map(part=>taExactEntityForConcept(part.conceptId))).size===4)??[];
 gastricAssemblyCache.set(parts,regions);return regions;
}

const targetIndexes=new WeakMap<Part[],{exact:Map<string,Part[]>;base:Map<string,Part[]>;concept:Map<string,Part[]>}>();
function targets(query:string,source:Part,parts:Part[],accept:(part:Part)=>boolean=()=>true,allBonePieces=false,allAliases=false,allConceptPieces=false):Part[]{
 const side=partSide(source);
 const q=query.replace('{side}',side??'left');
 const alternatives=query.includes('{side}')&&side===null?[q,query.replace('{side}','right')]:[q];
 let index=targetIndexes.get(parts);
 if(!index){index={exact:new Map(),base:new Map(),concept:new Map()};for(const part of parts){const exact=key(part.name),base=baseName(part.name);const exactGroup=index.exact.get(exact)??[];exactGroup.push(part);index.exact.set(exact,exactGroup);const baseGroup=index.base.get(base)??[];baseGroup.push(part);index.base.set(base,baseGroup);const conceptGroup=index.concept.get(part.conceptId)??[];conceptGroup.push(part);index.concept.set(part.conceptId,conceptGroup);}targetIndexes.set(parts,index);}
 const matches=(part:Part,alt:string)=>accept(part)&&(!sideOf(alt)||!partSide(part)||partSide(part)===sideOf(alt))&&(!side||!partSide(part)||partSide(part)===side);
 const exact=alternatives.flatMap(alt=>(index.exact.get(alt.toLowerCase())??[]).filter(part=>matches(part,alt)));
 const matched=exact.length&&!allAliases?exact:[...exact,...alternatives.flatMap(alt=>(index.base.get(baseName(alt))??[]).filter(part=>matches(part,alt)))];
 if(!matched.length&&baseName(q)==='stomach')matched.push(...modeledGastricAssembly(parts).filter(part=>matches(part,q)));
 const byConcept=new Map<string,Part>();
 for(const part of matched.filter(part=>part.id!==source.id)){
  const previous=byConcept.get(part.conceptId);
  if(!previous||part.vertexCount>previous.vertexCount)byConcept.set(part.conceptId,part);
 }
 const resolved=[...byConcept.values()];
 if(allConceptPieces)return resolved.flatMap(part=>(index.concept.get(part.conceptId)??[]).filter(candidate=>candidate.id!==source.id&&accept(candidate)&&sameSide(source,candidate)));
 return allBonePieces?resolved.flatMap(part=>isBone(part)?(index.concept.get(part.conceptId)??[]).filter(candidate=>candidate.id!==source.id&&isBone(candidate)&&accept(candidate)&&sameSide(source,candidate)):[part]):resolved;
}

function digestiveStages(parts:Part[]):string[][]{
 const hasWhole=(name:string)=>parts.some(part=>part.system==='digestive'&&baseName(part.name)===name.toLowerCase());
 const collapsed=(name:'jejunum'|'ileum')=>{
  const subdivisions=[`Proximal part of ${name}`,`Middle part of ${name}`,`Distal part of ${name}`];
  return hasWhole(name)||!subdivisions.some(label=>parts.some(part=>part.system==='digestive'&&baseName(part.name)===label.toLowerCase()))
   ?[[name[0].toUpperCase()+name.slice(1),...subdivisions]]:subdivisions.map(label=>[label]);
 };
 return [
  ...digestiveRoute.slice(0,5),
  ...collapsed('jejunum'),
  ...collapsed('ileum'),
  ...digestiveRoute.slice(11),
 ];
}

function routeRelations(source:Part,parts:Part[],routes:string[][][]):ResolvedRelation[]{
 const sourceKey=baseName(source.name);
 const intestinalPart=sourceKey.match(/^(?:proximal|middle|distal) part of (jejunum|ileum)$/);
 if(intestinalPart&&parts.some(part=>part.system==='digestive'&&baseName(part.name)===intestinalPart[1]))return [];
 const side=partSide(source);
 const result:ResolvedRelation[]=[];
 const sameSide=(group:string[])=>group.filter(label=>!side||!sideOf(label)||sideOf(label)===side);
 for(const route of routes){
  const index=route.findIndex(group=>group.some(label=>baseName(label)===sourceKey&&(!side||!sideOf(label)||sideOf(label)===side)));
  if(index<0)continue;
  for(const [kind,direction] of [['before',-1],['after',1]] as const){
   const via:string[]=[];
   for(let step=index+direction;step>=0&&step<route.length;step+=direction){
    const group=sameSide(route[step]);
    const wholeIntestine=['jejunum','ileum'].includes(baseName(group[0]??''));
    const labels=wholeIntestine&&targets(group[0],source,parts).length?[group[0]]:group;
    const matches=labels.flatMap(label=>targets(label,source,parts,undefined,false,true,true));
    if(matches.length){
     for(const target of matches)result.push({kind,target,...(via.length?{via:direction===1?via:[...via].reverse()}:{})});
     break;
    }
    if(group.length)via.push(group[0]);
   }
  }
 }
 return result;
}

function downstreamShortcuts(source:Part,parts:Part[]):ResolvedRelation[]{
 const shortcuts:[string,string,string][]=[
  ['Kidney','Renal pelvis','Ureter'],
  ['Testis','Epididymis','Ductus deferens'],
 ];
 const result:ResolvedRelation[]=[];
 for(const [first,middle,last] of shortcuts){
  const sourceKey=baseName(source.name);
  if(sourceKey!==first.toLowerCase()&&sourceKey!==last.toLowerCase())continue;
  const intermediary=targets(middle,source,parts);
  if(!intermediary.length)continue;
  const kind=sourceKey===first.toLowerCase()?'after':'before';
  for(const target of targets(kind==='after'?last:first,source,parts,undefined,false,false,true)){
   result.push({kind,target,via:[middle],viaModeled:true});
  }
 }
 return result;
}

function branchRelations(source:Part,parts:Part[],branches:[string[],string[][]][]):ResolvedRelation[]{
 const sourceKey=baseName(source.name),side=partSide(source);
 const result:ResolvedRelation[]=[];
 const matching=(group:string[])=>group.some(label=>baseName(label)===sourceKey&&(!side||!sideOf(label)||sideOf(label)===side));
 const sameSide=(group:string[])=>group.filter(label=>!side||!sideOf(label)||sideOf(label)===side);
 const overlaps=(a:string[],b:string[])=>sameSide(a).some(left=>sameSide(b).some(right=>baseName(left)===baseName(right)&&(!sideOf(left)||!sideOf(right)||sideOf(left)===sideOf(right))));
 const walk=(group:string[],kind:'before'|'after',via:string[],visited:Set<string>):void=>{
  const labels=sameSide(group);
  if(!labels.length)return;
  const identity=labels.map(label=>`${baseName(label)}:${sideOf(label)??''}`).join('|');
  if(visited.has(identity))return;
  const nextVisited=new Set(visited);nextVisited.add(identity);
  const matches=labels.flatMap(label=>targets(label,source,parts,undefined,false,true,true));
  if(matches.length){
   for(const target of matches)result.push({kind,target,...(via.length?{via:kind==='before'?[...via].reverse():via}:{})});
   return;
  }
  const nextVia=[...via,labels[0]];
  for(const [parent,children] of branches){
   if(kind==='after'&&overlaps(parent,labels))for(const child of children)walk(child,kind,nextVia,nextVisited);
   if(kind==='before')for(const child of children)if(overlaps(child,labels))walk(parent,kind,nextVia,nextVisited);
  }
 };
 for(const [parent,children] of branches){
  if(matching(parent))for(const child of children)walk(child,'after',[],new Set());
  for(const child of children)if(matching(child))walk(parent,'before',[],new Set());
 }
 return result;
}

function namedVascularRelations(source:Part,parts:Part[]):Relations{
 if(source.system!=='arterial'&&source.system!=='venous')return {};
 const result:Relations={},parentName=statedVascularParent(source.name);
 const childKind=source.system==='arterial'?'before':'after';
 const parentKind=source.system==='arterial'?'after':'before';
 if(parentName)result[childKind]=parts.filter(part=>part.system===source.system&&vascularIdentity(part.name)===vascularIdentity(parentName)).map(part=>part.name);
 const identity=vascularIdentity(source.name);
 result[parentKind]=parts.filter(part=>part.system===source.system&&part.id!==source.id&&
  statedVascularParent(part.name)&&vascularIdentity(statedVascularParent(part.name)!)===identity).map(part=>part.name);
 return result;
}

const sameSide=(a:Part,b:Part)=>!partSide(a)||!partSide(b)||partSide(a)===partSide(b);
const matchesBone=(bone:Part,label:string)=>isBone(bone)&&baseName(bone.name)===baseName(label);
const jointGroups=(part:Part)=>Object.keys(jointBones).filter(joint=>(part.groups??[]).some(group=>group.toLowerCase()===joint));

function physicalRelations(source:Part,parts:Part[]):Relations{
 const result:Relations={},sourceBase=baseName(source.name);
 // OpenStax A&P 2e 27.2: the ovarian ligament attaches the ovary to the uterus.
 // Keep this named reproductive support separate from skeletal ligaments;
 // targets() enforces the source's stated side and provides reciprocal links.
 if(source.system==='reproductive'){
  if(sourceBase==='ovarian ligament')result.connects=['Ovary','Uterus'];
  else if(sourceBase==='ovary'||sourceBase==='uterus')result.connectedBy=['Ovarian ligament'];
 }
 const joints=jointGroups(source);
 if(source.system==='connective'&&/ligament|capsule/i.test(source.name)){
  const specific=ligamentAttachments[sourceBase];
  result.connects=specific??[];
  result.joint=joints.flatMap(joint=>jointBones[joint]).filter(bone=>!result.connects?.includes(bone));
 }
 if(isBone(source)){
  const relatedJoints=Object.entries(jointBones).filter(([,bones])=>bones.some(bone=>matchesBone(source,bone))).map(([joint])=>joint);
  result.articulates=relatedJoints.flatMap(joint=>jointBones[joint].filter(bone=>!matchesBone(source,bone)));
  const jointStructures=parts.filter(part=>part.system==='connective'&&/ligament|capsule/i.test(part.name)&&sameSide(source,part)&&jointGroups(part).some(joint=>relatedJoints.includes(joint)));
  result.connectedBy=jointStructures.filter(part=>{
   const specific=ligamentAttachments[baseName(part.name)];
   if(specific)return specific.some(bone=>matchesBone(source,bone));
   return false;
  }).map(part=>part.name);
  result.joint=jointStructures.filter(part=>!result.connectedBy?.includes(part.name)).map(part=>part.name);
  for(const [tendon,bonesAndMuscles] of Object.entries(namedTendons)){
   if(bonesAndMuscles.some(label=>matchesBone(source,label)))result.connectedBy=[...(result.connectedBy??[]),...parts.filter(part=>baseName(part.name)===tendon&&sameSide(source,part)).map(part=>part.name)];
  }
  for(const [ligament,bones] of Object.entries(ligamentAttachments)){
   if(bones.some(label=>matchesBone(source,label)))result.connectedBy=[...(result.connectedBy??[]),...parts.filter(part=>baseName(part.name)===ligament&&sameSide(source,part)).map(part=>part.name)];
  }
 }
 const tendon=/\btendon\b/i.test(source.name)&&!/\btendon sheath\b/i.test(source.name);
 if(tendon){
  result.connects=[...(result.connects??[]),...(namedTendons[sourceBase]??[])];
  const of=clean(source.name).match(/\btendon of (.+?)(?: \((?:left|right)\))?$/i)?.[1];
  if(of)result.connects.push(of,`${of} muscle`);
  if(/·\s*Tendon$/i.test(source.name))result.connects.push(...parts.filter(part=>part.conceptId===source.conceptId&&part.id!==source.id&&!/·\s*Tendon$/i.test(part.name)).map(part=>part.name));
 }
 if(source.system==='fascia'||source.system==='connective'&&/fascia|aponeurosis/i.test(source.name)){
  result.continuous=fascialContinuity.filter(pair=>pair.some(name=>name.toLowerCase()===sourceBase)).flatMap(pair=>pair.filter(name=>name.toLowerCase()!==sourceBase));
  const covered=fascialCoverings[sourceBase]??[];
  result.covers=parts.filter(part=>part.system==='muscular'&&sameSide(source,part)&&covered.some(muscle=>baseName(part.name)===muscle.toLowerCase()||baseName(part.name).endsWith(` of ${muscle.toLowerCase()}`))).map(part=>part.name);
  if(sourceBase==='fascia lata')result.continuous.push('Iliotibial tract');
  if(sourceBase==='tensor fasciae latae'&&!isTendon(source))result.continuous.push('Iliotibial tract');
 }
 if(source.system==='muscular'&&!isTendon(source)){
  result.coveredBy=Object.entries(fascialCoverings).filter(([,muscles])=>muscles.some(muscle=>sourceBase===muscle.toLowerCase()||sourceBase.endsWith(` of ${muscle.toLowerCase()}`))).map(([fascia])=>fascia);
  result.connectedBy=parts.filter(part=>isTendon(part)&&part.id!==source.id&&sameSide(source,part)&&(part.conceptId===source.conceptId||baseName(part.name).replace(/^tendon of (?:left|right) /,'tendon of ')===`tendon of ${sourceBase}`)).map(part=>part.name);
  if(sourceBase==='tensor fasciae latae')result.continuous=['Iliotibial tract'];
 }
 if(sourceBase==='iliotibial tract')result.continuous=['Fascia lata','Tensor fasciae latae'];
 return result;
}

const networkCache=new WeakMap<Part,Map<string,ResolvedRelation[]>>();
const isCartilage=(part:Part)=>/cartilage/i.test(part.name)&&!part.name.toLowerCase().includes('perichondular');
const isTendon=(part:Part)=>/\btendon\b/i.test(part.name)&&!/\btendon sheath\b/i.test(part.name);
export const isAnatomicalBone=(part:Part)=>part.system==='skeletal'&&!isCartilage(part)&&!/perichondular|gingiva|\bteeth\b|\btooth\b|bone tissue|bone marrow|·\s*(?:suture|ligament)\b/i.test(part.name);
const isBone=isAnatomicalBone;
const canCarryFunctionalLinks=(part:Part)=>part.system!=='attachments'&&part.system!=='regions'&&!/·\s*(?:ligament|suture)\b/i.test(part.name)&&(part.system!=='muscular'||!isTendon(part))&&(part.system!=='skeletal'||isBone(part));
function explicitRelations(part:Part):Relations{
 if(!canCarryFunctionalLinks(part))return {};
 const name=withoutSide(part.name);
 if(femalePelvicFunctional[name]&&part.system!=='reproductive')return {};
 return explicit[name]??
  (['jejunum','ileum'].find(base=>belongsTo(name,base))?explicit[name.includes('jejunum')?'jejunum':'ileum']:undefined)??{};
}
const structuralCache=new WeakMap<Part,Map<string,ResolvedRelation[]>>();
const genericGroups=/^\d+:|^(?:head|neck|trunk|abdomen|thorax|pelvis|upper limb|lower limb|skin|muscles|nerves|arteries|veins|skeletal system|digestive system|respiratory system|central nervous system|peripheral nervous system|endocrine glands|visceral systems)$/i;
const sideNeutralSourceId=(part:Part)=>part.id.startsWith('HRA:')?(part.sourceId??part.id).replace(/_([LR])(?=_|$)/i,'_{side}'):'';
function structuralNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=structuralCache.get(parts[0]);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),byName=new Map<string,Part[]>(),byMuscleCore=new Map<string,Part[]>();
 for(const part of parts){const name=baseName(part.name),list=byName.get(name)??[];list.push(part);byName.set(name,list);if(part.system==='muscular'){const core=name.replace(/ muscle$/,'');const peers=byMuscleCore.get(core)??[];peers.push(part);byMuscleCore.set(core,peers);}}
 const add=(from:Part,kind:RelationKind,to:Part)=>{const list=graph.get(from.id)??[];if(!list.some(relation=>relation.kind===kind&&relation.target.id===to.id))list.push({kind,target:to});graph.set(from.id,list);};
 const link=(source:Part,parent:Part)=>{
  if(source.id===parent.id||source.conceptId===parent.conceptId||source.suppressed||parent.suppressed||!sameSide(source,parent))return;
  add(source,'partOf',parent);add(parent,'contains',source);
 };
 for(const source of parts){
  const side=partSide(source);
  if(side){
   const role=source.name.match(/\s*·\s*([^·]+)$/)?.[1]??'';
   const site=source.system==='attachments'?source.name.match(/\.([eo]\d*)[lr]$/i)?.[1]??'':'';
   const sourcePairKey=sideNeutralSourceId(source);
   const peers=(byName.get(baseName(source.name))??[]).filter(part=>part.system===source.system&&partSide(part)&&partSide(part)!==side&&(part.name.match(/\s*·\s*([^·]+)$/)?.[1]??'')===role&&(part.system!=='attachments'||(part.name.match(/\.([eo]\d*)[lr]$/i)?.[1]??'')===site)&&(!sourcePairKey||sideNeutralSourceId(part)===sourcePairKey));
   const opposite=peers.sort((a,b)=>b.vertexCount-a.vertexCount)[0];
   if(opposite){add(source,'counterpart',opposite);add(opposite,'counterpart',source);}
  }
  if(source.system==='attachments'){
   const site=source.name.match(/\.([oe])\d*[lr]$/i)?.[1]?.toLowerCase();
   if(site)for(const muscle of (byMuscleCore.get(baseName(source.name).replace(/ muscle$/,''))??[]).filter(part=>!/\btendon\b/i.test(part.name)&&sameSide(source,part)).sort((a,b)=>b.vertexCount-a.vertexCount).slice(0,1)){
    add(source,site==='o'?'originFor':'insertionFor',muscle);
    add(muscle,site==='o'?'originSites':'insertionSites',source);
   }
   continue;
  }
  const name=clean(source.name);
  const of=name.match(/\bof (.+)$/i)?.[1];
  const candidates=[...(of?[of,of.replace(/\s*\([IVX]+\)$/i,'')]:[]),...(source.groups??[]).filter(group=>!genericGroups.test(group)&&baseName(group)!==baseName(source.name))];
  const parents=new Map<string,Part>();
  for(const candidate of candidates){
   for(const parent of byName.get(baseName(candidate))??[]){
    if(parent.system!==source.system||!sameSide(source,parent)||parent.conceptId===source.conceptId)continue;
    const previous=parents.get(parent.conceptId);
    if(!previous||parent.vertexCount>previous.vertexCount)parents.set(parent.conceptId,parent);
   }
  }
  // Keep the closest named structures, not every broad source chapter.
  for(const parent of [...parents.values()].sort((a,b)=>clean(b.name).length-clean(a.name).length).slice(0,2))link(source,parent);
 }
 if(parts[0])structuralCache.set(parts[0],graph);
 return graph;
}
function skeletalNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=networkCache.get(parts[0]);
 if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(from:Part,kind:RelationKind,target:Part)=>{
  if(from.id===target.id)return;
  const links=graph.get(from.id)??[];
  if(!links.some(link=>link.kind===kind&&link.target.id===target.id))links.push({kind,target});
  graph.set(from.id,links);
 };
 const bones=parts.filter(isBone),cartilages=parts.filter(isCartilage);
 const muscleBodies=parts.filter(part=>part.system==='muscular'&&!isTendon(part)&&!/\b(?:bursa|sheath|ligament|retinaculum|fascia|aponeurosis|tendinous arch)\b/i.test(part.name));
 const muscleBodyIds=new Set(muscleBodies.map(part=>part.id));
 const tendons=parts.filter(isTendon);
 const fasciae=parts.filter(part=>part.system==='fascia'||baseName(part.name)==='iliotibial tract');
 const bonesByName=new Map<string,Part[]>();
 for(const bone of bones){const name=baseName(bone.name),group=bonesByName.get(name)??[];group.push(bone);bonesByName.set(name,group);}
 const muscleNameCache=new Map<string,Part[]>();
 const byConcept=new Map<string,Part[]>();
 for(const part of parts){
  const group=byConcept.get(part.conceptId)??[];
  group.push(part);byConcept.set(part.conceptId,group);
 }
 const canonical=(candidates:Part[])=>{
  const found=new Map<string,Part>();
  for(const part of candidates){
   const prior=found.get(part.conceptId);
   if(!prior||part.vertexCount>prior.vertexCount)found.set(part.conceptId,part);
  }
  return [...found.values()];
 };
 const boneNamed=(name:string,side:'left'|'right'|null)=>canonical((bonesByName.get(baseName(name))??[]).filter(part=>!side||!partSide(part)||partSide(part)===side));
 const muscleNamed=(name:string,side:'left'|'right'|null)=>{
  const wanted=baseName(name).replace(/ muscle$/,'');
  let candidates=muscleNameCache.get(wanted);
  if(!candidates){candidates=muscleBodies.filter(part=>{
   const key=baseName(part.name).replace(/ muscle$/,'');
   return key===wanted||key.endsWith(` of ${wanted}`);
  });muscleNameCache.set(wanted,candidates);}
  return canonical(candidates.filter(part=>!side||!partSide(part)||partSide(part)===side));
 };
 for(const muscle of muscleBodies)for(const tendon of (byConcept.get(muscle.conceptId)??[]).filter(isTendon)){
  add(muscle,'tendons',tendon);add(tendon,'muscles',muscle);
 }

 // Material-split bone surfaces share a concept ID; a named costal cartilage
 // instead joins the matching rib. Neither rule links cartilage to vessels.
 for(const cartilage of cartilages){
  const peers=(byConcept.get(cartilage.conceptId)??[]).filter(isBone);
  let attached=canonical(peers);
  if(!attached.length){
   const name=baseName(cartilage.name);
   const rib=name.match(/(?:costal cartilage of (first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|eleventh|twelfth) rib|^(first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|eleventh|twelfth) costal cartilage$)/);
   if(rib)attached=boneNamed(`${rib[1]??rib[2]} rib`,partSide(cartilage));
   else if(name.startsWith('articular cartilage of')&&name.includes('femur'))attached=boneNamed('Femur',partSide(cartilage));
   else if(name==='triradiate cartilage')attached=boneNamed('Hip bone',partSide(cartilage));
  }
  for(const bone of attached){add(bone,'cartilages',cartilage);add(cartilage,'bones',bone);}
 }

 // The source model has separate origin and insertion markers. Their measured
 // contact with bone is generated from the actual imported mesh coordinates.
 for(const site of skeletalAttachments.entries){
  const side=site.side as 'left'|'right';
  const muscles=muscleNamed(site.muscle,side);
  if(!muscles.length)continue;
  for(const muscle of muscles){
   const tendonSurfaces=(byConcept.get(muscle.conceptId)??[]).filter(isTendon);
   for(const boneName of site.bones)for(const bone of boneNamed(boneName,side)){
    add(muscle,site.site==='origin'?'origin':'insertion',bone);
    add(bone,'muscles',muscle);
    if(site.site==='insertion')for(const tendon of tendonSurfaces){add(bone,'tendons',tendon);add(tendon,'bones',bone);}
   }
  }
 }
 // Preserve the few named tendons whose attachment is expressed in the name
 // rather than by a separate Z-Anatomy marker.
 for(const [name,connections] of Object.entries(namedTendons))for(const tendon of tendons.filter(part=>baseName(part.name)===name)){
  const attachedBones=canonical(connections.flatMap(connection=>boneNamed(connection,partSide(tendon))));
  const attachedMuscles=canonical([
   ...(byConcept.get(tendon.conceptId)??[]).filter(part=>muscleBodyIds.has(part.id)),
   ...connections.flatMap(connection=>muscleNamed(connection,partSide(tendon))),
  ]);
  for(const bone of attachedBones){add(tendon,'bones',bone);add(bone,'tendons',tendon);}
  for(const muscle of attachedMuscles){add(tendon,'muscles',muscle);add(muscle,'tendons',tendon);}
  for(const bone of attachedBones)for(const muscle of attachedMuscles){add(muscle,'insertion',bone);add(bone,'muscles',muscle);}
 }
 for(const [name,connections] of Object.entries(tendonFascialAttachments))for(const tendon of tendons.filter(part=>baseName(part.name)===name)){
  for(const connection of connections)for(const fascia of targets(connection,tendon,fasciae)){
   add(tendon,'fascia',fascia);add(fascia,'tendons',tendon);
  }
 }
 for(const fascia of fasciae)for(const boneName of fascialBoneAttachments[baseName(fascia.name)]??[]){
  for(const bone of boneNamed(boneName,partSide(fascia))){add(fascia,'bones',bone);add(bone,'fascia',fascia);}
 }
 const arteries=parts.filter(part=>part.system==='arterial');
 const veins=parts.filter(part=>part.system==='venous');
 const nerves=parts.filter(part=>part.system==='nervous');
 for(const bone of canonical(bones)){
  const name=baseName(bone.name);
  const supply=boneSupply[name]??regionalBoneSupply(name)??(/^(?:first|second|third|fourth|fifth|sixth|seventh|eighth|ninth|tenth|eleventh|twelfth) rib$/.test(name)
   ?{arterial:['Posterior intercostal arteries'],innervation:['Intercostal nerves']}:undefined);
  if(!supply)continue;
  for(const query of supply.arterial??[])for(const artery of targets(query,bone,arteries)){add(bone,'arterial',artery);add(artery,'supplies',bone);}
  for(const query of supply.venous??[])for(const vein of targets(query,bone,veins)){add(bone,'venous',vein);add(vein,'drains',bone);}
  for(const query of supply.innervation??[])for(const nerve of targets(query,bone,nerves)){add(bone,'innervation',nerve);add(nerve,'innervates',bone);}
 }
 for(const muscle of canonical(muscleBodies)){
  const muscleName=baseName(muscle.name).replace(/ muscle$/,'');
  for(const [name,supply] of Object.entries(muscleArterialSupply)){
   const targetName=name.replace(/ muscle$/,'');
   if(muscleName!==targetName&&!muscleName.endsWith(` of ${targetName}`))continue;
   for(const query of supply)for(const artery of targets(query,muscle,arteries)){
    add(muscle,'arterial',artery);add(artery,'supplies',muscle);
   }
  }
 }
 // A single anatomical bone or muscle can have several imported material
 // pieces. Selecting any piece must show the same concept-level relationships.
 for(const group of byConcept.values())for(const members of [group.filter(isBone),group.filter(part=>muscleBodyIds.has(part.id))]){
  if(members.length<2)continue;
  const union=members.flatMap(part=>graph.get(part.id)??[]);
  const memberIds=new Set(members.map(part=>part.id));
  for(const member of members)for(const relation of union){
   const reverse=(graph.get(relation.target.id)??[]).filter(link=>memberIds.has(link.target.id));
   add(member,relation.kind,relation.target);
   for(const link of reverse)add(relation.target,link.kind,member);
  }
 }
 if(parts[0])networkCache.set(parts[0],graph);
 return graph;
}

function nerveParent(source:Part,identities:Set<string>):string|null{
 if(source.system!=='nervous'||!/\bnerve\b/i.test(source.name))return null;
 const sourceSide=sideOf(source.name);
 const exists=(label:string)=>identities.has(`${baseName(label)}|${sourceSide??''}`)||!sourceSide&&[...identities].some(identity=>identity.startsWith(`${baseName(label)}|`));
 for(const [parent,children] of nerveBranches){
  if(children.some(group=>group.some(label=>baseName(label)===baseName(source.name))))return parent.find(exists)??null;
 }
 const lexical=clean(source.name).match(/\b(?:branch(?:es)?|division) of (.+? nerve(?: \([IVX]+\))?)(?: \((?:left|right)\))?$/i)?.[1];
 return lexical&&exists(lexical)?lexical:null;
}

function nerveRelations(source:Part,parts:Part[]):Relations{
 if(source.system!=='nervous')return {};
 const candidates=parts.filter(part=>part.system==='nervous');
 const identities=new Set(candidates.map(part=>`${baseName(part.name)}|${sideOf(part.name)??''}`));
 const result:Relations={},parent=nerveParent(source,identities);
 if(parent)result.before=[parent];
 result.after=candidates.filter(part=>part.id!==source.id&&sameSide(source,part)&&
  baseName(nerveParent(part,identities)??'')===baseName(source.name)).map(part=>part.name);
 const muscleName=clean(source.name).match(/^Nerve to (.+?)(?: \((?:left|right)\))?$/i)?.[1];
 if(muscleName)result.innervates=[muscleName,`${muscleName} muscle`];
 return result;
}

function muscleNerves(source:Part):string[]{
 if(source.system!=='muscular'||isTendon(source))return [];
 const groups=source.groups??[];
 const upper=groups.some(group=>/upper limb|muscles of hand/i.test(group));
 const lower=groups.some(group=>/lower limb|muscles of foot/i.test(group));
 const upperNerves=/^(?:Ulnar|Median|Radial|Musculocutaneous|Axillary|Lateral pectoral|Long thoracic|Suprascapular|Dorsal scapular|Subclavian) nerve$/i;
 const lowerNerves=/^(?:Sciatic|Tibial|Common fibular|Deep fibular|Superficial fibular|Medial plantar|Lateral plantar|Femoral|Obturator) nerve$/i;
 const groupNerves=source.system==='muscular'&&!/\btendon\b/i.test(source.name)?groups.filter(group=>
  / nerve(?: \([IVX]+\))?$/i.test(group)&&!(lower&&upperNerves.test(group))&&!(upper&&lowerNerves.test(group))):[];
 let direct=groupNerves.length>1?groupNerves.filter(group=>group!=='Sciatic nerve'):groupNerves;
 const name=baseName(source.name);
 if(direct.length===1&&direct[0]==='Sciatic nerve'&&groups.some(group=>/muscles of foot|posterior compartment of leg/i.test(group)))direct=[];
 if(name==='popliteus muscle')direct=['Tibial nerve'];
 if(name==='flexor digitorum brevis')direct=['Medial plantar nerve'];
 if(name==='long head of biceps femoris')direct=direct.filter(group=>group!=='Common fibular nerve');
 if(name==='short head of biceps femoris')direct=direct.filter(group=>group!=='Tibial nerve');
 return direct;
}

function functionalNote(source:Part,kind:RelationKind,target:Part):string|undefined{
 if(source.system!=='reproductive'||withoutSide(source.name)!=='vagina'||kind!=='innervation')return;
 if(withoutSide(target.name)==='pudendal nerve')return 'Lower vagina · somatic supply';
 if(withoutSide(target.name)==='uterovaginal plexus')return 'Upper vagina · autonomic supply';
}

function reverseFunctionalRelations(source:Part,parts:Part[]):ResolvedRelation[]{
 if(!['nervous','arterial','venous'].includes(source.system))return [];
 const result:ResolvedRelation[]=[];
 if(source.system==='nervous')for(const part of parts.filter(part=>part.system==='muscular'&&sameSide(source,part)&&muscleNerves(part).some(nerve=>baseName(nerve)===baseName(source.name))))result.push({kind:'innervates',target:part});
 const kind=source.system==='nervous'?'innervation':source.system==='arterial'?'arterial':'venous';
 const reverse=source.system==='nervous'?'innervates':source.system==='arterial'?'supplies':'drains';
 for(const part of parts){
  const lookup=explicitRelations(part);
  if(!lookup?.[kind]||!sameSide(source,part))continue;
  if(lookup[kind].some(query=>targets(query,part,parts).some(target=>target.conceptId===source.conceptId)))result.push({kind:reverse,target:part,note:functionalNote(part,kind,source)});
 }
 return result;
}

function accepts(kind:RelationKind,source:Part,target:Part):boolean{
 if(kind==='innervation')return target.system==='nervous';
 if(kind==='arterial')return target.system==='arterial';
 if(kind==='venous')return target.system==='venous';
 if(kind==='innervates')return target.system==='muscular'&&!/·\s*Tendon$/i.test(target.name)||['digestive','respiratory','cardiac','urinary','reproductive','skeletal'].includes(target.system);
 if(kind==='supplies'||kind==='drains')return !['arterial','venous','nervous','attachments','regions'].includes(target.system);
 if(kind==='articulates')return isBone(target);
 if(kind==='connects')return isBone(target)||target.system==='muscular'||source.system==='reproductive'&&baseName(source.name)==='ovarian ligament'&&target.system==='reproductive'&&['ovary','uterus'].includes(baseName(target.name));
 if(kind==='connectedBy')return target.system==='connective'||target.system==='muscular'&&/\btendon\b/i.test(target.name)||source.system==='reproductive'&&['ovary','uterus'].includes(baseName(source.name))&&target.system==='reproductive'&&baseName(target.name)==='ovarian ligament';
 if(kind==='joint')return isBone(source)?target.system==='connective':isBone(target);
 if(kind==='continuous')return target.system==='fascia'||target.system==='connective'||baseName(target.name)==='iliotibial tract';
 if(kind==='covers')return target.system==='muscular';
 if(kind==='coveredBy')return target.system==='fascia';
 return true;
}

export function anatomicalRelations(source:Part,allParts:Iterable<Part>):ResolvedRelation[]{
 const parts=Array.isArray(allParts)?allParts:[...allParts];
 const lookup=explicitRelations(source);
 const namedMuscleNerve=source.system==='muscular'&&!isTendon(source)?[`Nerve to ${baseName(source.name).replace(/ muscle$/,'')} muscle`]:[];
 const mapped:Relations={...lookup,innervation:[...(lookup.innervation??[]),...muscleNerves(source),...namedMuscleNerve]};
 const route=routeRelations(source,parts,[digestiveStages(parts),airwayRoute,leftAirwayRoute,arterialRoute,legArterialRoute,legVenousRoute,armArterialRoute,armVenousRoute,portalVenousRoute,portalMesentericRoute,greatSaphenousRoute,smallSaphenousRoute,rightHeartRoute,inferiorCavaRoute,leftHeartRoute,urinaryRoute,spermRoute]);
 const branches=branchRelations(source,parts,[...airwayBranches,...arterialBranches,...nerveBranches]);
 const namedVessels=namedVascularRelations(source,parts);
 const nerves=nerveRelations(source,parts);
 const physical=physicalRelations(source,parts);
 const reverse=reverseFunctionalRelations(source,parts);
 const result:ResolvedRelation[]=[];
 const seen=new Set<string>();
 const needsSkeletalNetwork=['skeletal','muscular','fascia','connective','arterial','venous','nervous'].includes(source.system);
 for(const relation of [...route,...branches,...downstreamShortcuts(source,parts),...(structuralNetwork(parts).get(source.id)??[]),...(needsSkeletalNetwork?skeletalNetwork(parts).get(source.id)??[]:[]),...reverse]){
  if(!accepts(relation.kind,source,relation.target))continue;
  const token=`${relation.kind}:${relation.target.id}`;
  if(!seen.has(token)){seen.add(token);result.push(relation);}
 }
 for(const kind of ['before','after','innervation','arterial','venous','innervates','supplies','drains','articulates','connects','connectedBy','joint','continuous','covers','coveredBy'] as const){
  for(const query of [...(namedVessels[kind]??[]),...(nerves[kind]??[]),...(physical[kind]??[]),...(mapped[kind]??[])]){
   for(const target of targets(query,source,parts,part=>accepts(kind,source,part),kind==='articulates'||kind==='connects'||kind==='joint',false,['before','after','innervation','arterial','venous'].includes(kind))){
    const token=`${kind}:${target.id}`;
    if(!seen.has(token)){seen.add(token);result.push({kind,target,note:functionalNote(source,kind,target)});}
   }
  }
 }
 return result;
}

/** Whole-concept functional links for a complete, explicitly named assembly.
 * This is a semantic query against the existing registry, never a geometry part
 * or an assertion that every region has the whole organ's vascular supply.
 */
export function anatomicalConceptFunctionalRelations(concept:Concept,selected:Part[],parts:Part[]):ResolvedRelation[]{
 if(selected.length<2||concept.id==='selection-set'||concept.id.startsWith('hierarchy:')||concept.id.startsWith('guest:'))return [];
 const elements=new Set(concept.elements),ids=new Set(selected.map(part=>part.id));
 if(elements.size!==ids.size||[...elements].some(id=>!ids.has(id)))return [];
 if(new Set(selected.map(part=>part.system)).size!==1)return [];
 const lookup=explicit[withoutSide(concept.name)];if(!lookup)return [];
 const source={...selected[0],id:`functional-query:${concept.id}`,conceptId:concept.id,name:concept.name};
 if(femalePelvicFunctional[withoutSide(concept.name)]&&source.system!=='reproductive')return [];
 if(!canCarryFunctionalLinks(source))return [];
 const result:ResolvedRelation[]=source.system==='digestive'?routeRelations(source,parts,[digestiveStages(parts)]).filter(relation=>!ids.has(relation.target.id)):[],seen=new Set<string>(result.map(relation=>`${relation.kind}:${relation.target.id}`));
 for(const kind of ['innervation','arterial','venous'] as const)for(const query of lookup[kind]??[]){
  for(const target of targets(query,source,parts,part=>accepts(kind,source,part),false,false,true)){
   if(ids.has(target.id))continue;
   const key=`${kind}:${target.id}`;if(seen.has(key))continue;seen.add(key);result.push({kind,target,note:functionalNote(source,kind,target)});
  }
 }
 return result;
}
