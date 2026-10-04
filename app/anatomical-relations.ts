import {endocrineRelations} from './endocrine-relations';
import {digestiveDuctRelations} from './digestive-duct-relations';
import {hepaticRelations,isNativeWholeLiver} from './hepatic-relations';
import {cranialRelations} from './cranial-relations';
import {cardiacRelations} from './cardiac-relations';
import {pulmonaryRelations,isReviewedPulmonaryCovering} from './pulmonary-relations';
import {coronaryRelations,isCoronaryArtery} from './coronary-relations';
import {containmentExceptions,isAorticBranchGrouping,nonComponentRelationship} from './anatomical-containment';
import {nativePartFrame} from './ta98-modeled-parent';
import {fibularRelations} from './fibular-relations';
import {handVascularRelations} from './hand-vascular-relations';
import {hasNativeUrinaryRoute,nativeUrinaryRelations,urinaryPartFrame} from './urinary-relations';
import {structureName,type Part,type Concept} from './anatomy';
import {taExactEntityForConcept} from './anatomical-terminology';
import skeletalAttachments from './generated-skeletal-attachments.json';

export type RelationKind='before'|'after'|'innervation'|'arterial'|'venous'|'innervates'|'supplies'|'drains'|'articulates'|'connects'|'connectedBy'|'joint'|'continuous'|'covers'|'coveredBy'|'cartilages'|'bones'|'tendons'|'fascia'|'muscles'|'origin'|'insertion'|'partOf'|'contains'|'originFor'|'insertionFor'|'originSites'|'insertionSites'|'counterpart'|'adjacent'|'lymphaticDrainage'|'lymphaticTributaries';
export interface ResolvedRelation {kind:RelationKind;target:Part;via?:string[];viaModeled?:boolean;viaRoutes?:{via:string[];modeled:boolean}[];note?:string}
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
 [['Common palmar digital nerves of median nerve'],[['Proper palmar digital nerves of median nerve']]],
 [['Superficial branch of ulnar nerve'],[['Common palmar digital nerves of ulnar nerve']]],
 // The female source is a distal digital network, not a palmar cutaneous nerve.
 [['Superficial branch of ulnar nerve'],[['Palmar branch ulnar nerve']]],
 [['Common palmar digital nerves of ulnar nerve'],[['Proper palmar digital nerves of ulnar nerve']]],
 [['Deep branch of radial nerve'],[['Posterior interosseous nerve of forearm']]],
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
 // Whole-organ links only. NCBI Bookshelf NBK519558 (diaphragm),
 // NBK470452 (thyroid), NBK459205 and
 // PMID 3376104 (appendicular artery/vein). Missing vessels stay unresolved.
 'diaphragm':{innervation:['Phrenic nerve'],arterial:['Musculophrenic artery','Superior phrenic arteries','Inferior phrenic artery','Inferior phrenic arteries']},
 'thyroid gland':{arterial:['Superior thyroid artery','Inferior thyroid artery'],venous:['Superior thyroid vein','Middle thyroid vein','Inferior thyroid vein']},
 'vermiform appendix':{arterial:['Appendicular artery'],venous:['Appendicular vein']},
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
 'adrenal gland':{arterial:['Superior suprarenal artery','Middle suprarenal artery','Inferior suprarenal artery'],venous:['Suprarenal vein']},
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
// Naming and branching reference: NCBI Bookshelf NBK532968. These synonyms
// resolve the existing native nerve records without inferring mesh junctions.
const anatomicalNameAliases:Record<string,string>={'common peroneal nerve':'common fibular nerve','deep peroneal nerve':'deep fibular nerve','superficial peroneal nerve':'superficial fibular nerve','peroneal artery':'fibular artery','peroneal veins':'fibular veins','intestine duodenum':'duodenum','small intestine jejunum':'jejunum','small intestine ileum':'ileum','small intestine illium':'ileum','large intestine cecum':'cecum','large intestine descending colon':'descending colon','large intestine rectum':'rectum','appendix':'vermiform appendix','bladder':'urinary bladder'};
// Qualified source labels retain the parent nerve and its sensory territory.
Object.assign(anatomicalNameAliases,{
 'palmar branch ulnar nerve':'unpartitioned palmar digital network of ulnar nerve',
 'common palmar digital branches median nerve':'common palmar digital nerves of median nerve',
 'common palmar digital branches of median nerve':'common palmar digital nerves of median nerve',
 'proper palmar digital branches median nerve':'proper palmar digital nerves of median nerve',
 'proper palmar digital branches of median nerve':'proper palmar digital nerves of median nerve',
 'common palmar digital branches of ulnar nerve':'common palmar digital nerves of ulnar nerve',
 'proper palmar digital branches of ulnar nerve':'proper palmar digital nerves of ulnar nerve',
 'deep branch radial nerve':'deep branch of radial nerve',
 'superficial branch radial nerve':'superficial branch of radial nerve',
});
const baseName=(name:string)=>{const value=/^(?:left|right) (?:atrium|ventricle)$/.test(withoutSide(name))?withoutSide(name):withoutSide(name).replace(/^(?:left|right)\s+/,'');return anatomicalNameAliases[value]??value;};
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

// A whole pharynx remains a whole structure; this route fallback does not
// rename it as an oropharynx or manufacture its missing regional meshes.
const digestiveStageCache=new WeakMap<Part[],string[][]>();
const airwayStageCache=new WeakMap<Part[],string[][]>();
function pharyngealStages(parts:Part[]):string[][]{
 return parts.some(part=>!part.suppressed&&['oropharynx','laryngopharynx'].includes(baseName(part.name)))
  ?digestiveRoute.slice(0,2):[['Pharynx']];
}
function airwayStages(parts:Part[]):string[][]{
 const cached=airwayStageCache.get(parts);if(cached)return cached;
 const stages=[airwayRoute[0],...pharyngealStages(parts),...airwayRoute.slice(3)];
 airwayStageCache.set(parts,stages);return stages;
}
function digestiveStages(parts:Part[]):string[][]{
 const cached=digestiveStageCache.get(parts);if(cached)return cached;
 const hasWhole=(name:string)=>parts.some(part=>part.system==='digestive'&&baseName(part.name)===name.toLowerCase());
 const collapsed=(name:'jejunum'|'ileum')=>{
  const subdivisions=[`Proximal part of ${name}`,`Middle part of ${name}`,`Distal part of ${name}`];
  return hasWhole(name)||!subdivisions.some(label=>parts.some(part=>part.system==='digestive'&&baseName(part.name)===label.toLowerCase()))
   ?[[name[0].toUpperCase()+name.slice(1),...subdivisions]]:subdivisions.map(label=>[label]);
 };
 const stages=[
  ['Mouth'],['Oral cavity'],...pharyngealStages(parts),...digestiveRoute.slice(2,5),
  ...collapsed('jejunum'),
  ...collapsed('ileum'),
  ...digestiveRoute.slice(11),
 ];
 digestiveStageCache.set(parts,stages);return stages;
}

const routeConceptNames=new WeakMap<Part[],Map<string,string[]>>();
function routeRelations(source:Part,parts:Part[],routes:string[][][]):ResolvedRelation[]{
 const sourceKey=baseName(source.name);
 let conceptNames=routeConceptNames.get(parts);
 if(!conceptNames){conceptNames=new Map();for(const part of parts){const names=conceptNames.get(part.conceptId)??[];const name=baseName(part.name);if(!names.includes(name))names.push(name);conceptNames.set(part.conceptId,names);}routeConceptNames.set(parts,conceptNames);}
 const wholeGastricRegion=source.system==='digestive'&&modeledGastricAssembly(parts).some(part=>part.id===source.id);
 const sourceKeys=wholeGastricRegion?['stomach']:[sourceKey,...(conceptNames.get(source.conceptId)??[])];
 const intestinalPart=sourceKey.match(/^(?:proximal|middle|distal) part of (jejunum|ileum)$/);
 if(intestinalPart&&parts.some(part=>part.system==='digestive'&&baseName(part.name)===intestinalPart[1]))return [];
 const side=partSide(source);
 const result:ResolvedRelation[]=[];
 const sameSide=(group:string[])=>group.filter(label=>!side||!sideOf(label)||sideOf(label)===side);
 for(const route of routes){
  if(route===urinaryRoute&&hasNativeUrinaryRoute(source))continue;
  const directIndex=route.findIndex(group=>group.some(label=>baseName(label)===sourceKey&&(!side||!sideOf(label)||sideOf(label)===side)));
  const index=directIndex>=0?directIndex:route.findIndex(group=>group.some(label=>sourceKeys.includes(baseName(label))&&(!side||!sideOf(label)||sideOf(label)===side)));
  if(index<0)continue;
  for(const [kind,direction] of [['before',-1],['after',1]] as const){
   const via:string[]=[];
   for(let step=index+direction;step>=0&&step<route.length;step+=direction){
    const group=sameSide(route[step]);
    const wholeIntestine=['jejunum','ileum'].includes(baseName(group[0]??''));
    const labels=wholeIntestine&&targets(group[0],source,parts).length?[group[0]]:group;
    const matches=labels.flatMap(label=>targets(label,source,parts,undefined,false,true,true)).filter(target=>{
     if(route!==arterialRoute)return true;
     const frame=nativePartFrame(source);return !!frame&&!target.suppressed&&nativePartFrame(target)===frame;
    });
    if(matches.length){
     for(const target of matches){const wholePharynx=sourceKey==='pharynx'||baseName(target.name)==='pharynx';const wholeRoute=wholeGastricRegion||!route[index].some(label=>baseName(label)===sourceKey)||modeledGastricAssembly(parts).some(part=>part.id===target.id);result.push({kind,target,...(via.length?{via:direction===1?via:[...via].reverse()}:{}),...(wholePharynx?{note:'Route through the pharynx’s oral/laryngeal regions; the named source surface has unverified regional extent and tissue continuity.'}:wholeRoute?{note:'Route through the whole named structure; this surface piece does not establish direct local continuity.'}:{})});}
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
  if(first==='Kidney'&&hasNativeUrinaryRoute(source))continue;
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
 if(source.system==='muscular'){
  const independentTissue=!isTendon(source)||!parts.some(part=>part.id!==source.id&&part.conceptId===source.conceptId&&!isTendon(part));
  if(independentTissue)result.coveredBy=Object.entries(fascialCoverings).filter(([,muscles])=>muscles.some(muscle=>sourceBase===muscle.toLowerCase()||sourceBase.endsWith(` of ${muscle.toLowerCase()}`))).map(([fascia])=>fascia);
  if(!isTendon(source))result.connectedBy=parts.filter(part=>isTendon(part)&&part.id!==source.id&&sameSide(source,part)&&(part.conceptId===source.conceptId||baseName(part.name).replace(/^tendon of (?:left|right) /,'tendon of ')===`tendon of ${sourceBase}`)).map(part=>part.name);
  if(sourceBase==='tensor fasciae latae')result.continuous=['Iliotibial tract'];
 }
 if(sourceBase==='iliotibial tract')result.continuous=['Fascia lata','Tensor fasciae latae'];
 return result;
}

const networkCache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const isCartilage=(part:Part)=>/cartilage/i.test(part.name)&&!part.name.toLowerCase().includes('perichondular');
const isTendon=(part:Part)=>/\btendon\b/i.test(part.name)&&!/\btendon sheath\b/i.test(part.name);
export const isAnatomicalBone=(part:Part)=>part.system==='skeletal'&&!isCartilage(part)&&!/perichondular|gingiva|\bteeth\b|\btooth\b|bone tissue|bone marrow|·\s*(?:suture|ligament)\b/i.test(part.name);
const isBone=isAnatomicalBone;
const canCarryFunctionalLinks=(part:Part)=>part.system!=='attachments'&&part.system!=='regions'&&!/·\s*(?:ligament|suture)\b/i.test(part.name)&&(part.system!=='muscular'||!isTendon(part))&&(part.system!=='skeletal'||isBone(part));
const isAdrenal=(part:Part)=>part.system==='endocrine'&&['adrenal gland','suprarenal gland'].includes(baseName(part.name));
const adrenalFunctionalFrame=(source:Part,target:Part)=>{
 if(!isAdrenal(source)&&!isAdrenal(target))return true;
 const frame=nativePartFrame(source);
 return !source.suppressed&&!target.suppressed&&!!frame&&nativePartFrame(target)===frame&&sameSide(source,target);
};
function explicitRelations(part:Part):Relations{
 if(!canCarryFunctionalLinks(part))return {};
 const name=withoutSide(part.name);
 if(femalePelvicFunctional[name]&&part.system!=='reproductive')return {};
 return explicit[name]??explicit[baseName(part.name)]??
  (['jejunum','ileum'].find(base=>belongsTo(name,base))?explicit[name.includes('jejunum')?'jejunum':'ileum']:undefined)??{};
}
const structuralCache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const genericGroups=/^\d+:|^(?:head|neck|trunk|abdomen|thorax|pelvis|upper limb|lower limb|skin|muscles|nerves|arteries|veins|skeletal system|digestive system|respiratory system|central nervous system|peripheral nervous system|endocrine glands|visceral systems)$/i;
const sideNeutralSourceId=(part:Part)=>part.id.startsWith('HRA:')?(part.sourceId??part.id).replace(/_([LR])(?=_|$)/i,'_{side}'):'';
function structuralNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=structuralCache.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),byName=new Map<string,Part[]>(),byMuscleCore=new Map<string,Part[]>();
 for(const part of parts){const name=baseName(part.name),list=byName.get(name)??[];list.push(part);byName.set(name,list);if(part.system==='muscular'){const core=name.replace(/ muscle$/,'');const peers=byMuscleCore.get(core)??[];peers.push(part);byMuscleCore.set(core,peers);}}
 const add=(from:Part,kind:RelationKind,to:Part,note?:string)=>{const list=graph.get(from.id)??[];if(!list.some(relation=>relation.kind===kind&&relation.target.id===to.id))list.push({kind,target:to,note});graph.set(from.id,list);};
 const link=(source:Part,parent:Part)=>{
  if(source.id===parent.id||source.conceptId===parent.conceptId||source.suppressed||parent.suppressed||!sameSide(source,parent))return;
  if(isAorticBranchGrouping(source.name,parent.name))return;
  if(baseName(source.name)==='proximal segment of subclavian artery'&&baseName(parent.name)==='subclavian artery'){
   const frame=nativePartFrame(source);if(!frame||nativePartFrame(parent)!==frame)return;
  }
  const note=nonComponentRelationship(source.name,parent.name);
  if(note){add(source,'adjacent',parent,note);add(parent,'adjacent',source,note);return;}
  add(source,'partOf',parent);add(parent,'contains',source);
 };
 for(const source of parts){
  const side=partSide(source);
  // A left-labeled brachiocephalic source may be a variant or misidentified
  // segment; it is not a bilateral counterpart of the conventional trunk.
  if(side&&baseName(source.name)!=='brachiocephalic trunk'){
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
 for(const rule of containmentExceptions){
  for(const parent of byName.get(rule.parent)??[])for(const child of byName.get(rule.child)??[]){
   const frame=nativePartFrame(parent);
   if(parent.system!=='nervous'||child.system!=='nervous'||!frame||nativePartFrame(child)!==frame||!sameSide(parent,child))continue;
   add(parent,'adjacent',child,rule.note);add(child,'adjacent',parent,rule.note);
  }
 }
 structuralCache.set(parts,graph);
 return graph;
}
function skeletalNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=networkCache.get(parts);
 if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>();
 const add=(from:Part,kind:RelationKind,target:Part,note?:string)=>{
  if(from.id===target.id)return;
  const links=graph.get(from.id)??[];
  if(!links.some(link=>link.kind===kind&&link.target.id===target.id))links.push({kind,target,note});
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
 for(const disc of parts.filter(part=>!part.suppressed&&part.system==='connective'&&baseName(part.name)==='articular disc of distal radio-ulnar joint')){
  const frame=nativePartFrame(disc);if(!frame)continue;
  for(const name of ['Radius','Ulna'])for(const bone of boneNamed(name,partSide(disc)).filter(bone=>!bone.suppressed&&nativePartFrame(bone)===frame)){
   const note='Structures at the distal radioulnar joint; the source disc is a coarse surface and its full fibrocartilage attachments and interfaces are unverified.';
   add(disc,'joint',bone,note);add(bone,'joint',disc,note);
  }
 }
 const upperCervicalBones:Record<string,string[]>={
  'anterior atlanto-occipital membrane':['Atlas vertebra c1','Atlas (C1)'],
  'posterior atlanto-occipital membrane':['Atlas vertebra c1','Atlas (C1)'],
  'tectorial membrane of atlanto-axial joint':['Axis vertebra c2','Axis (C2)'],
  'alar ligaments':['Axis vertebra c2','Axis (C2)'],
 };
 for(const tissue of parts.filter(part=>!part.suppressed&&part.system==='connective')){
  const names=upperCervicalBones[baseName(tissue.name)],frame=nativePartFrame(tissue);
  if(!names||!frame)continue;
  for(const name of names)for(const bone of boneNamed(name,null).filter(bone=>!bone.suppressed&&nativePartFrame(bone)===frame)){
   const note='Named craniocervical bone association. Exact tissue extent and entheses remain unverified; the occipital attachment has no separate native bone mesh in this source. This link does not certify a shared mesh interface.';
   add(tissue,'joint',bone,note);add(bone,'joint',tissue,note);
  }
 }
 for(const muscle of muscleBodies)for(const tendon of (byConcept.get(muscle.conceptId)??[]).filter(isTendon)){
  add(muscle,'tendons',tendon);add(tendon,'muscles',muscle);
 }
 // Resolve the same named attachment queries used by the tendon inspector.
 // Common/intermediate tendon names and alternate source muscle names need
 // not share a concept ID or the literal "tendon of" prefix.
 for(const tendon of tendons)for(const query of physicalRelations(tendon,parts).connects??[]){
  for(const muscle of targets(query,tendon,parts,part=>muscleBodyIds.has(part.id),false,false,false)){
   add(tendon,'connects',muscle);add(muscle,'connectedBy',tendon);
  }
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
 networkCache.set(parts,graph);
 return graph;
}

function nerveParent(source:Part,identities:Set<string>):string|null{
 if(source.system!=='nervous'||!/\bnerves?\b/i.test(source.name))return null;
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
 if(muscleName){const core=muscleName.replace(/ muscle$/i,'');result.innervates=[core,`${core} muscle`];}
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

/** Existing curated supply labels; no inferred spinal level or territory inheritance. */
export function anatomicalInnervationNames(source:Part):string[]{
 return [...new Set([...(explicitRelations(source).innervation??[]),...muscleNerves(source)])];
}

function functionalNote(source:Part,kind:RelationKind,target:Part):string|undefined{
 if(['before','after'].includes(kind)&&source.system==='nervous'&&target.system==='nervous'){
  const names=[baseName(source.name),baseName(target.name)];
  if(names.includes('unpartitioned palmar digital network of ulnar nerve'))return 'Provisional distal digital network through the superficial ulnar branch. Common/proper digital subdivisions and the physical junction are unverified; the palmar cutaneous branch is a distinct structure.';
  if(names.some(name=>/^(?:common|proper) palmar digital nerves of (?:median|ulnar) nerve$/.test(name)))return 'Named digital branching relationship, preserving median versus ulnar territory. Source fascicles, individual sensory territories and physical junctions remain unverified; this is not a one-direction axonal-flow claim.';
  if(names.includes('deep branch of radial nerve'))return 'Deep radial motor/proprioceptive branch. Its source extent and posterior interosseous transition remain unverified; no cutaneous sensory territory is inferred.';
  if(names.includes('superficial branch of radial nerve'))return 'Superficial radial cutaneous branch. Individual skin territories, variants and physical source junctions remain unverified.';
 }
 if(isBone(source)&&baseName(source.name)==='fibula'&&kind==='arterial'){
  if(baseName(target.name)==='fibular artery')return 'Fibular (peroneal) arterial nutrient and periosteal branches supply the fibular shaft. Individual source branches, nutrient foramina and mesh contacts remain unverified.';
  if(baseName(target.name)==='anterior tibial artery')return 'Proximal fibular head/epiphyseal arterial branches; the named artery is not a reconstructed local nutrient branch.';
 }
 if(isAdrenal(source)&&['arterial','venous'].includes(kind))return 'Named adrenal vascular reference. Individual capsular branches, source ostia, continuous lumens and variants remain unverified; missing native vessels are not replaced by donor meshes.';
 if(source.system!=='reproductive'||withoutSide(source.name)!=='vagina'||kind!=='innervation')return;
 if(withoutSide(target.name)==='pudendal nerve')return 'Lower vagina · somatic supply';
 if(withoutSide(target.name)==='uterovaginal plexus')return 'Upper vagina · autonomic supply';
}

// Atlas arrays are immutable snapshots. Build inverse supply links once per
// snapshot instead of scanning every organ again for each vessel or nerve.
const reverseFunctionalCache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
function reverseFunctionalRelations(source:Part,parts:Part[]):ResolvedRelation[]{
 if(!['nervous','arterial','venous'].includes(source.system))return [];
 let graph=reverseFunctionalCache.get(parts);
 if(!graph){
  graph=new Map();
  const byConcept=new Map<string,Part[]>(),nervesByName=new Map<string,Part[]>();
  for(const part of parts){
   const peers=byConcept.get(part.conceptId)??[];peers.push(part);byConcept.set(part.conceptId,peers);
   if(part.system==='nervous'){const name=baseName(part.name),nerves=nervesByName.get(name)??[];nerves.push(part);nervesByName.set(name,nerves);}
  }
  const add=(from:Part,kind:RelationKind,target:Part,note?:string)=>{
   const links=graph!.get(from.id)??[];
   if(!links.some(link=>link.kind===kind&&link.target.id===target.id))links.push({kind,target,...(note?{note}:{})});
   graph!.set(from.id,links);
  };
  for(const part of parts){
   if(part.system==='muscular')for(const query of muscleNerves(part))for(const nerve of nervesByName.get(baseName(query))??[]){
    if(sameSide(nerve,part))add(nerve,'innervates',part);
   }
   const lookup=explicitRelations(part);
   for(const [kind,system,reverse] of [['innervation','nervous','innervates'],['arterial','arterial','supplies'],['venous','venous','drains']] as const){
    const concepts=new Set((lookup[kind]??[]).flatMap(query=>targets(query,part,parts,target=>adrenalFunctionalFrame(part,target)).map(target=>target.conceptId)));
    for(const concept of concepts)for(const target of byConcept.get(concept)??[]){
     if(target.system===system&&sameSide(target,part)&&adrenalFunctionalFrame(part,target))add(target,reverse,part,functionalNote(part,kind,target));
    }
   }
  }
  reverseFunctionalCache.set(parts,graph);
 }
 return graph.get(source.id)??[];
}

const isRetina=(part:Part)=>['nervous','sensory'].includes(part.system)&&baseName(part.name)==='retina';
const isOpticNerve=(part:Part)=>['nervous','sensory'].includes(part.system)&&['optic nerve','optic nerve (ii)'].includes(baseName(part.name));
// OpenStax A&P 14.1: retinal ganglion axons leave through the optic nerve.
// PMC11890264: central retinal artery supplies inner retina; outer retina has
// a separate choroidal supply. Links describe anatomy, not mesh registration.
function retinalRelations(source:Part,parts:Part[]):ResolvedRelation[]{
 if(source.suppressed)return [];
 const retina=isRetina(source),optic=isOpticNerve(source);
 const artery=source.system==='arterial'&&baseName(source.name)==='central retinal artery';
 const vein=source.system==='venous'&&baseName(source.name)==='central retinal vein';
 if(!retina&&!optic&&!artery&&!vein)return [];
 const sourceFamily=(part:Part)=>part.id.startsWith('ZA:')?'ZA':part.id.startsWith('HRA:')?'HRA':part.id.match(/^LOCAL:(male|female):/)?.[0]??part.id;
 const side=partSide(source);
 if(!side)return [];
 const links:ResolvedRelation[]=[];
 for(const target of parts){
  if(target.suppressed||target.id===source.id||partSide(target)!==side||sourceFamily(target)!==sourceFamily(source))continue;
  let kind:RelationKind|undefined,note:string|undefined;
  if(retina&&isOpticNerve(target)||optic&&isRetina(target)){
   kind='continuous';note='Retinal ganglion-cell axons collect at the optic disc and continue into the optic nerve.';
  }else if(retina&&target.system==='arterial'&&baseName(target.name)==='central retinal artery'||artery&&isRetina(target)){
   kind=retina?'arterial':'supplies';note='The central retinal artery supplies the inner retina; outer retinal supply is choroidal.';
  }else if(retina&&target.system==='venous'&&baseName(target.name)==='central retinal vein'||vein&&isRetina(target)){
   kind=retina?'venous':'drains';note='Retinal venous blood drains through the central retinal vein.';
  }
  if(kind)links.push({kind,target,note:source.id.startsWith('HRA:')?`${note} These donor surfaces still require registration and packing review.`:note});
 }
 return links;
}

function reviewedCranialConnection(source:Part,target:Part):boolean{
 const frame=nativePartFrame(source),names=[baseName(source.name),baseName(target.name)];
 return source.system==='nervous'&&target.system==='nervous'&&!!frame&&nativePartFrame(target)===frame&&!!partSide(source)&&partSide(target)===partSide(source)&&names.includes('chorda tympani')&&names.some(name=>['lingual nerve','facial nerve','facial nerve (vii)'].includes(name));
}
function reviewedEpiglotticAttachment(ligament:Part,tissue:Part):boolean{
 const frame=nativePartFrame(ligament),name=baseName(ligament.name).replace(/-/g,''),target=baseName(tissue.name);
 if(!frame||nativePartFrame(tissue)!==frame||ligament.suppressed||tissue.suppressed)return false;
 return name==='thyroepiglottic ligament'&&['epiglottis','thyroid cartilage'].includes(target)
  ||name==='hyoepiglottic ligament'&&['epiglottis','hyoid bone','body of hyoid bone'].includes(target);
}
function reviewedPelvicConnection(source:Part,target:Part):boolean{
 if(source.suppressed||target.suppressed||source.system!=='reproductive'||target.system!=='reproductive')return false;
 const frame=nativePartFrame(source);if(!frame||nativePartFrame(target)!==frame)return false;
 const cervix=(part:Part)=>['cervix','uterine cervix','cervix uteri'].includes(baseName(part.name));
 return cervix(source)&&baseName(target.name)==='vagina'||cervix(target)&&baseName(source.name)==='vagina';
}

function accepts(kind:RelationKind,source:Part,target:Part):boolean{
 if(['before','after','arterial','supplies','venous','drains','innervation','innervates'].includes(kind)&&[source,target].some(p=>isNativeWholeLiver(p)||['Hepatic portal vein','Hepatic veins'].includes(p.name))){
  const frame=nativePartFrame(source),other=nativePartFrame(target);
  if((frame||other)&&(source.suppressed||target.suppressed||!frame||frame!==other))return false;
 }
 if(['before','after'].includes(kind)&&[source,target].some(p=>p.system==='respiratory'&&/trachea|bronchus|bronchi/.test(baseName(p.name)))){
  const frame=nativePartFrame(source),other=nativePartFrame(target);
  if((frame||other)&&(source.suppressed||target.suppressed||!frame||frame!==other))return false;
 }
 if(['before','after','adjacent'].includes(kind)&&[source,target].some(isCoronaryArtery)){
  const frame=nativePartFrame(source),other=nativePartFrame(target);
  if(frame||other)if(source.suppressed||target.suppressed||kind!=='adjacent'&&(source.system!=='arterial'||target.system!=='arterial')||!frame||frame!==other)return false;
 }
 if(['before','after'].includes(kind)&&[source,target].some(p=>p.system==='urinary')){
  const frame=urinaryPartFrame(source),other=urinaryPartFrame(target);
  if((frame||other)&&(source.suppressed||target.suppressed||!frame||frame!==other))return false;
 }
 // UAMS upper-limb nerve table: these are named branching associations,
 // not source-tissue containment or certified mesh junctions.
 if(['before','after'].includes(kind)&&[source,target].some(part=>/^(?:(?:median|ulnar|radial) nerve|(?:common|proper) palmar digital nerves of (?:median|ulnar) nerve|(?:deep|superficial) branch of (?:radial|ulnar) nerve|unpartitioned palmar digital network of ulnar nerve)$/.test(baseName(part.name)))){
  const frame=nativePartFrame(source);
  if(source.suppressed||target.suppressed||source.system!=='nervous'||target.system!=='nervous'||!frame||nativePartFrame(target)!==frame)return false;
 }

 if([source,target].some(part=>part.id.startsWith('LOCAL:female:ta98-fibula:'))&&(source.suppressed||target.suppressed||nativePartFrame(source)!==nativePartFrame(target)))return false;
 // Peroneal is the source synonym of fibular, not an additional nerve.
 // Resolving that alias must not connect independently registered donors.
 if(['before','after','innervation','innervates'].includes(kind)&&[source,target].some(part=>part.system==='nervous'&&/^(?:common|deep|superficial) fibular nerve$/.test(baseName(part.name)))){
  if(source.suppressed||target.suppressed)return false;
  if(['before','after'].includes(kind)&&(source.system!=='nervous'||target.system!=='nervous'))return false;
  const frame=nativePartFrame(source),other=nativePartFrame(target);
  if((frame||other)&&frame!==other)return false;
 }
 if(['arterial','supplies'].includes(kind)){
  const fibula=[source,target].find(p=>isBone(p)&&/^(?:fibula|.+ of fibula)$/.test(baseName(p.name)));
  if(fibula){const frame=nativePartFrame(fibula);if(!frame||nativePartFrame(source)!==frame||nativePartFrame(target)!==frame)return false;}
 }
 if(['arterial','venous','innervation','supplies','drains','innervates'].includes(kind)&&!adrenalFunctionalFrame(source,target))return false;
 if(kind==='innervation')return target.system==='nervous';
 if(kind==='arterial')return target.system==='arterial';
 if(kind==='venous')return target.system==='venous';
 if(kind==='innervates')return target.system==='muscular'&&!/·\s*Tendon$/i.test(target.name)||['digestive','respiratory','cardiac','urinary','reproductive','skeletal'].includes(target.system);
 if(kind==='supplies'||kind==='drains')return isRetina(target)||!['arterial','venous','nervous','attachments','regions'].includes(target.system);
 if(kind==='articulates')return isBone(target);
 if(kind==='connects')return reviewedEpiglotticAttachment(source,target)||isBone(target)||target.system==='muscular'||source.system==='reproductive'&&baseName(source.name)==='ovarian ligament'&&target.system==='reproductive'&&['ovary','uterus'].includes(baseName(target.name));
 if(kind==='connectedBy')return reviewedEpiglotticAttachment(target,source)||target.system==='connective'||target.system==='muscular'&&/\btendon\b/i.test(target.name)||source.system==='reproductive'&&['ovary','uterus'].includes(baseName(source.name))&&target.system==='reproductive'&&baseName(target.name)==='ovarian ligament';
 if(kind==='joint')return isBone(source)?target.system==='connective':isBone(target);
 if(kind==='continuous')return reviewedPelvicConnection(source,target)||reviewedCranialConnection(source,target)||isRetina(source)&&isOpticNerve(target)||isOpticNerve(source)&&isRetina(target)||target.system==='fascia'||target.system==='connective'||baseName(target.name)==='iliotibial tract';
 if(kind==='covers')return isReviewedPulmonaryCovering(source,target)||target.system==='muscular';
 if(kind==='coveredBy')return isReviewedPulmonaryCovering(target,source)||target.system==='fascia';
 return true;
}

// Anatomical flow relationships in a known native source frame. They do not
// certify mesh junctions, ostia or patient-specific terminal anatomy.
// References: PMID17415746 (cisterna investigation), PMID39124550 (terminal
// variation). An unspecified 'Lymphatic duct' is not assumed to be right-sided.
const lymphaticNetworks=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const renalVenousNetworks=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const aorticBranchNetworks=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const adrenalNetworks=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
function adrenalNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=adrenalNetworks.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),kidneys=parts.filter(part=>!part.suppressed&&part.system==='urinary'&&baseName(part.name)==='kidney');
 // PMID31501707 distinguishes normally separate adrenal/renal capsules;
 // PMID18246295 reviews their sectional regional anatomy. This is adjacency,
 // never a shared-tissue interface, containment or source packing acceptance.
 for(const gland of parts.filter(part=>!part.suppressed&&isAdrenal(part))){
  const frame=nativePartFrame(gland),side=partSide(gland);if(!frame||!side)continue;
  for(const kidney of kidneys.filter(part=>nativePartFrame(part)===frame&&partSide(part)===side)){
   const note='Distinct ipsilateral organs in the perinephric region. Their capsules are normally separate; this link does not imply renal containment, adrenal-renal fusion or verified source packing.';
   for(const [from,to] of [[gland,kidney],[kidney,gland]]){
    const list=graph.get(from.id)??[];list.push({kind:'adjacent',target:to,note});graph.set(from.id,list);
   }
  }
 }
 adrenalNetworks.set(parts,graph);return graph;
}
function aorticBranchNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=aorticBranchNetworks.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),arteries=parts.filter(part=>!part.suppressed&&part.system==='arterial');
 const add=(source:Part,target:Part,note='Conventional aortic arch branching reference. Source ostia, continuous lumens and anatomical variants are unverified; this relationship does not certify a mesh junction.')=>{
  for(const [from,to,kind] of [[source,target,'after'],[target,source,'before']] as const){
   const links=graph.get(from.id)??[];if(!links.some(link=>link.kind===kind&&link.target.id===to.id))links.push({kind,target:to,note});graph.set(from.id,links);
  }
 };
 for(const source of arteries){
  const frame=nativePartFrame(source);if(!frame)continue;
  const name=baseName(source.name);
  if(name==='common carotid artery'){
   const side=partSide(source);if(!side)continue;
   for(const target of arteries.filter(part=>nativePartFrame(part)===frame&&partSide(part)===side&&['internal carotid artery','external carotid artery'].includes(baseName(part.name))))add(source,target,'Common carotid bifurcation into the ipsilateral internal and external carotid arteries. Source bifurcation surfaces, lumen continuity and anatomical variants remain unverified.');
   continue;
  }
  if(name==='proximal segment of subclavian artery'&&partSide(source)==='left'){
   for(const target of arteries.filter(part=>nativePartFrame(part)===frame&&baseName(part.name)==='subclavian artery'&&partSide(part)==='left'))add(source,target,'Adjoining native source partitions share 12 distal boundary coordinates. This establishes the source surface interface, not biological lumen continuity or complete vessel walls.');
   continue;
  }
  const branch=name==='aortic arch'||name==='arch of aorta'?'arch':name==='brachiocephalic trunk'&&partSide(source)!=='left'?'trunk':undefined;
  if(!branch)continue;
  for(const target of arteries){
   if(nativePartFrame(target)!==frame)continue;
   const key=baseName(target.name),side=partSide(target);
   if(branch==='arch'&&side==='left'){
    if(key==='proximal segment of subclavian artery'){add(source,target);continue;}
    if(key==='subclavian artery'&&arteries.some(part=>nativePartFrame(part)===frame&&baseName(part.name)==='proximal segment of subclavian artery'&&partSide(part)==='left'))continue;
   }
   if(branch==='arch'?(key==='brachiocephalic trunk'&&side!=='left'||['common carotid artery','subclavian artery'].includes(key)&&side==='left'):['common carotid artery','subclavian artery'].includes(key)&&side==='right')add(source,target);
  }
 }
 aorticBranchNetworks.set(parts,graph);return graph;
}
function renalVenousNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=renalVenousNetworks.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>();
 const veins=parts.filter(part=>!part.suppressed&&part.system==='venous');
 for(const source of veins){
  const side=/^intrarenal veins of (left|right) kidney$/i.exec(source.name)?.[1]?.toLowerCase();
  const frame=nativePartFrame(source);if(!side||!frame)continue;
  for(const target of veins.filter(part=>baseName(part.name)==='renal vein'&&partSide(part)===side&&nativePartFrame(part)===frame)){
   const note='Intrarenal venous tributaries drain into the ipsilateral renal vein. This named source network does not establish every tributary, lumen or a continuous mesh junction.';
   for(const [from,to,kind] of [[source,target,'after'],[target,source,'before']] as const){
    const list=graph.get(from.id)??[];list.push({kind,target:to,note});graph.set(from.id,list);
   }
  }
 }
 renalVenousNetworks.set(parts,graph);return graph;
}
function lymphaticNetwork(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=lymphaticNetworks.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),named=new Map<string,Part[]>();
 for(const part of parts){if(part.suppressed||!['lymphatic','venous'].includes(part.system))continue;const key=baseName(part.name),list=named.get(key)??[];list.push(part);named.set(key,list);}
 const add=(source:Part,target:Part,note:string,via?:string[])=>{
  for(const [from,to,kind] of [[source,target,'lymphaticDrainage'],[target,source,'lymphaticTributaries']] as const){
   const list=graph.get(from.id)??[];if(!list.some(link=>link.kind===kind&&link.target.id===to.id))list.push({kind,target:to,note,via,viaModeled:via?false:undefined});graph.set(from.id,list);
  }
 };
 const sameFrame=(source:Part,name:string,system:string,side?:'left'|'right')=>{
  const frame=nativePartFrame(source);if(!frame)return [];
  const candidates=(named.get(name)??[]).filter(part=>part.system===system&&nativePartFrame(part)===frame&&(!side||partSide(part)===side));
  const largest=new Map<string,Part>();for(const part of candidates){const prior=largest.get(part.conceptId);if(!prior||part.vertexCount>prior.vertexCount)largest.set(part.conceptId,part);}return [...largest.values()];
 };
 for(const spleen of named.get('spleen')??[]){
  if(spleen.system!=='lymphatic')continue;
  for(const name of ['splenic nodes','splenic lymph nodes'])for(const nodes of sameFrame(spleen,name,'lymphatic'))add(spleen,nodes,'Regional splenic lymph drainage reference at the capsule and hilum. This source group does not resolve each collector, hilar node or efferent junction; source interfaces remain unverified. Blood-borne splenic cell traffic is not represented by this link.');
 }
 for(const cisterna of named.get('cisterna chyli')??[]){if(cisterna.system!=='lymphatic')continue;for(const duct of sameFrame(cisterna,'thoracic duct','lymphatic'))add(cisterna,duct,'Lymph passes from the modeled cisterna chyli into the thoracic duct. Cisterna anatomy varies; surface continuity and lumen junction are unverified.');}
 for(const duct of named.get('thoracic duct')??[]){
  if(duct.system!=='lymphatic')continue;
  for(const name of ['internal jugular vein','subclavian vein'])for(const vein of sameFrame(duct,name,'venous','left'))add(duct,vein,'Typical left jugulosubclavian termination, represented here by the two adjoining venous segments. These are references to one venous-angle region, not two independent outlets. Terminal anatomy varies; source mesh junctions are unverified.',['Left venous angle']);
 }
 for(const duct of named.get('lymphatic duct')??[]){
  if(duct.system!=='lymphatic'||partSide(duct)!=='right')continue;
  for(const name of ['internal jugular vein','subclavian vein'])for(const vein of sameFrame(duct,name,'venous','right'))add(duct,vein,'Typical right jugulosubclavian termination, represented by its adjoining venous segments. The terminal junction and its anatomical variants are unverified.',['Right venous angle']);
 }
 lymphaticNetworks.set(parts,graph);return graph;
}

const reviewedRouteCache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
function reviewedNativeRoutes(parts:Part[]):Map<string,ResolvedRelation[]>{
 const cached=reviewedRouteCache.get(parts);if(cached)return cached;
 const graph=new Map<string,ResolvedRelation[]>(),available=parts.filter(p=>!p.suppressed);
 const add=(source:Part,target:Part,kind:RelationKind,inverse:RelationKind,note:string)=>{
  const frame=nativePartFrame(source);if(!frame||nativePartFrame(target)!==frame||source.id===target.id)return;
  for(const [from,to,k] of [[source,target,kind],[target,source,inverse]] as const){const list=graph.get(from.id)??[];if(!list.some(r=>r.kind===k&&r.target.id===to.id))list.push({kind:k,target:to,note});graph.set(from.id,list);}
 };
 const named=(source:Part,names:string[],system:string,side?:'left'|'right'|null)=>available.filter(p=>p.system===system&&nativePartFrame(p)===nativePartFrame(source)&&names.includes(baseName(p.name))&&(!side||partSide(p)===side));
 for(const source of available){
  const name=baseName(source.name),side=partSide(source);if(!nativePartFrame(source))continue;
  // IARC cervical anatomy, NCBI Bookshelf NBK568392: the vaginal
  // portion of the cervix enters the vaginal vault, with the fornices around it.
  if(source.system==='reproductive'&&['cervix','uterine cervix','cervix uteri'].includes(name)){
   for(const target of named(source,['vagina'],'reproductive'))add(source,target,'continuous','continuous','Cervix–vagina anatomical junction: the vaginal portion of the cervix enters the vaginal vault, surrounded by the fornices; the cervical canal opens into the vagina at the external os. Source wall attachment, shared coordinates and lumen continuity remain unverified.');
  }
  // NCBI Bookshelf NBK538202: named epiglottic ligament attachments.
  // Resolve native peers only; independently registered donor ligaments must
  // not appear attached to the main body's epiglottis merely by name.
  const epiglotticLigament=name.replace(/-/g,'');
  if(['hyoepiglottic ligament','thyroepiglottic ligament'].includes(epiglotticLigament)){
   const note='Named epiglottic ligament attachment reference. Source tissue extent, entheses and shared mesh interfaces remain unverified.';
   for(const target of named(source,['epiglottis'],'respiratory'))add(source,target,'connects','connectedBy',note);
   if(epiglotticLigament==='thyroepiglottic ligament'){
    for(const target of available.filter(p=>nativePartFrame(p)===nativePartFrame(source)&&baseName(p.name)==='thyroid cartilage'))add(source,target,'connects','connectedBy',note);
   }else{
    for(const target of available.filter(p=>nativePartFrame(p)===nativePartFrame(source)&&['hyoid bone','body of hyoid bone'].includes(baseName(p.name))))add(source,target,'connects','connectedBy',note);
   }
  }
  if(name==='chorda tympani'&&source.system==='nervous'&&side){
   for(const target of named(source,['facial nerve','facial nerve (vii)','lingual nerve'],'nervous',side))add(source,target,'continuous','continuous','Chorda tympani branches from the facial/intermediate nerve and joins the ipsilateral lingual nerve. It carries taste afferents and parasympathetic efferents; this connection is not a one-direction axonal-flow claim. Source fascicles and physical junctions remain unverified.');
  }
  if(name==='ascending lumbar vein'&&source.system==='venous'&&side){
   const outlets=side==='right'?['azygos vein']:['hemi-azygos vein','hemiazygos vein'];
   for(const target of named(source,outlets,'venous'))add(source,target,'after','before','Named ascending-lumbar collateral drainage into the ipsilateral azygos/hemiazygos system. Caudal venous roots and subcostal communications vary; a direct end-to-end mesh continuation or universal junction is not asserted. Source lumen junctions remain unverified.');
  }
  if(name==='lumbar veins'&&source.system==='venous'&&side){
   for(const target of named(source,['inferior vena cava','inferior vena cava (abdominal part)'],'venous'))add(source,target,'after','before','Lumbar venous drainage to the inferior vena cava; tributary levels and variants in this grouped source are unverified.');
   for(const target of named(source,['ascending lumbar vein'],'venous',side))add(source,target,'after','before','Segmental lumbar veins communicate with the ipsilateral ascending lumbar collateral pathway. This grouped surface does not verify each segmental tributary or lumen junction.');
  }
  if(name==='dorsal digital arteries of foot'&&source.system==='arterial'&&side){
   for(const target of named(source,['dorsal metatarsal arteries'],'arterial',side))add(target,source,'after','before','Dorsal metatarsal branches continue into dorsal digital arteries of the toes. Individual branch coverage, source wall layers and lumen junctions remain unverified.');
  }
 }
 reviewedRouteCache.set(parts,graph);return graph;
}
export function anatomicalRelations(source:Part,allParts:Iterable<Part>):ResolvedRelation[]{
 const parts=Array.isArray(allParts)?allParts:[...allParts];
 const lookup=explicitRelations(source);
 const namedMuscleNerve=source.system==='muscular'&&!isTendon(source)?[`Nerve to ${baseName(source.name).replace(/ muscle$/,'')} muscle`]:[];
 const mapped:Relations={...lookup,innervation:[...(lookup.innervation??[]),...muscleNerves(source),...namedMuscleNerve]};
 const route=routeRelations(source,parts,[digestiveStages(parts),airwayStages(parts),leftAirwayRoute,arterialRoute,legArterialRoute,legVenousRoute,armArterialRoute,armVenousRoute,portalVenousRoute,portalMesentericRoute,greatSaphenousRoute,smallSaphenousRoute,rightHeartRoute,inferiorCavaRoute,leftHeartRoute,urinaryRoute,spermRoute]);
 const branches=branchRelations(source,parts,[...airwayBranches,...arterialBranches,...nerveBranches]);
 const namedVessels=namedVascularRelations(source,parts);
 const nerves=nerveRelations(source,parts);
 const physical=physicalRelations(source,parts);
 const reverse=reverseFunctionalRelations(source,parts);
 const result:ResolvedRelation[]=[];
 const seen=new Set<string>();
 const needsSkeletalNetwork=['skeletal','muscular','fascia','connective','arterial','venous','nervous'].includes(source.system);
 for(const relation of [...(endocrineRelations(parts).get(source.id)??[]),...(digestiveDuctRelations(parts).get(source.id)??[]),...(hepaticRelations(parts).get(source.id)??[]),...route,...branches,...downstreamShortcuts(source,parts),...(structuralNetwork(parts).get(source.id)??[]),...(needsSkeletalNetwork?skeletalNetwork(parts).get(source.id)??[]:[]),...reverse,...retinalRelations(source,parts),...(lymphaticNetwork(parts).get(source.id)??[]),...(renalVenousNetwork(parts).get(source.id)??[]),...(aorticBranchNetwork(parts).get(source.id)??[]),...(adrenalNetwork(parts).get(source.id)??[]),...(reviewedNativeRoutes(parts).get(source.id)??[]),...(cranialRelations(parts).get(source.id)??[]),...(cardiacRelations(parts).get(source.id)??[]),...(pulmonaryRelations(parts).get(source.id)??[]),...(coronaryRelations(parts).get(source.id)??[]),...(fibularRelations(parts).get(source.id)??[]),...(handVascularRelations(parts).get(source.id)??[]),...(nativeUrinaryRelations(parts).get(source.id)??[])]){
  if(!accepts(relation.kind,source,relation.target))continue;
  const token=`${relation.kind}:${relation.target.id}`;
  if(!seen.has(token)){seen.add(token);const note=relation.note??functionalNote(source,relation.kind,relation.target);result.push(note?{...relation,note}:relation);}
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
 // The whole native fibula can be assembled from regional and surface leaves.
 // Its supply belongs to the stored complete assembly, not each named leaf.
 const fibularWhole=baseName(concept.name)==='fibula'&&/^(?:ZA:Fibula\.[lr]|LOCAL:(?:male|female):[lr]_fibula)$/.test(concept.id)
  &&selected.every(part=>part.system==='skeletal'&&part.sectionAssembly===concept.id&&!part.suppressed)
  &&parts.filter(part=>!part.suppressed&&part.sectionAssembly===concept.id).every(part=>ids.has(part.id));
 const lookup=explicit[withoutSide(concept.name)]??(fibularWhole?boneSupply.fibula:undefined);if(!lookup)return [];
 const source={...selected[0],id:fibularWhole?concept.id:`functional-query:${concept.id}`,conceptId:concept.id,name:concept.name};
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
