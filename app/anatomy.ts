export type SystemId = 'skeletal'|'muscular'|'arterial'|'venous'|'nervous'|'digestive'|'respiratory'|'urinary'|'reproductive'|'lymphatic'|'endocrine'|'integumentary'|'connective'|'sensory'|'cardiac'|'pregnancy'|'attachments'|'regions'|'fascia'|'schematic'|'cell-boundary'|'cell-nucleus'|'cell-endomembrane'|'cell-mitochondria'|'cell-cytoskeleton'|'cell-vesicles'|'cell-protein-synthesis'|'cell-division';
export type AnatomySex='male'|'female';
export const SYSTEMS: {id:SystemId;name:string;color:string;description:string}[] = [
 {id:'skeletal',name:'Skeleton',color:'#e2d9ba',description:'Bones form the supporting framework of the body, protect organs, and provide attachment points for muscles. Their internal tissue also stores minerals and produces blood cells.'},
 {id:'muscular',name:'Muscles',color:'#a85b50',description:'Skeletal muscles generate movement by pulling on their attachments. Together with tendons, they move joints, stabilize posture, and produce heat.'},
 {id:'cardiac',name:'Heart',color:'#b96760',description:'The heart is a muscular pump with four chambers. Its valves direct blood forward through the pulmonary and systemic circuits.'},
 {id:'sensory',name:'Sensory organs',color:'#b0c8ce',description:'These structures contribute to special senses, including sight, hearing, and balance. Their specialized tissues detect stimuli and work with the nervous system to convey information.'},
 {id:'arterial',name:'Arteries',color:'#c05245',description:'The heart drives blood through the circulation. Arteries carry blood away from the heart to supply tissues or, in the pulmonary circuit, to the lungs.'},
 {id:'venous',name:'Veins',color:'#527c9f',description:'Veins return blood toward the heart. Superficial and deep networks collect blood from the tissues; the pulmonary veins bring oxygenated blood back from the lungs.'},
 {id:'nervous',name:'Nervous system',color:'#d8b565',description:'The brain, spinal cord, and peripheral nerves carry and process signals. They support sensation, movement, coordination, and automatic regulation of body functions.'},
 {id:'respiratory',name:'Respiratory',color:'#b98991',description:'The airways conduct air to the lungs, where oxygen and carbon dioxide move between air and blood. Breathing depends on pressure changes produced by respiratory muscles.'},
 {id:'digestive',name:'Digestive',color:'#b8916b',description:'The digestive tract breaks down food, absorbs nutrients and water, and moves waste onward. Accessory organs contribute bile and digestive enzymes.'},
 {id:'urinary',name:'Urinary',color:'#b47961',description:'The kidneys filter blood and regulate fluid, electrolyte, and acid–base balance. Urine travels through the ureters to the bladder and exits through the urethra.'},
 {id:'lymphatic',name:'Lymphatic',color:'#879f7c',description:'Lymphatic vessels return excess tissue fluid to the circulation. Lymph nodes and other lymphoid organs support immune surveillance and responses.'},
 {id:'endocrine',name:'Endocrine',color:'#c5a09a',description:'Endocrine organs release hormones into the blood to coordinate processes such as metabolism, growth, stress responses, and reproduction.'},
 {id:'reproductive',name:'Reproductive',color:'#bda098',description:'The male reproductive structures represented here contribute to sperm production, maturation, transport, and the production of sex hormones.'},
 {id:'integumentary',name:'Body surface',color:'#ba9b7d',description:'The body surface provides an outer anatomical reference. The integumentary system forms a protective barrier and contributes to sensation and temperature regulation.'},
 {id:'pregnancy',name:'Pregnancy reference',color:'#b88380',description:'The placenta and umbilical cord support exchange between maternal and fetal circulations during pregnancy. These reference structures are shown separately from the default adult anatomy.'},
 {id:'fascia',name:'Fascia',color:'#aec3bb',description:'Connective-tissue sheets surround muscles and other structures. Hide this layer to inspect the anatomy beneath it.'},
 {id:'schematic',name:'Schematic details',color:'#d8b565',description:'Source-designated schematic anatomy for isolated study. Its position and scale do not represent placement within the assembled body.'},
 {id:'attachments',name:'Muscle attachments',color:'#bb896b',description:'Source-marked regions where muscles originate or insert. These are attachment annotations, rather than additional muscles.'},
 {id:'regions',name:'Surface regions',color:'#aaa69c',description:'Named regions on the body surface provide topographical context. They are reference areas rather than separate organs.'},
 {id:'connective',name:'Connective tissue',color:'#aec3bb',description:'Cartilage, ligaments, and other connective tissues support, connect, and separate structures. Their roles include stabilizing joints and distributing mechanical loads.'},
 {id:'cell-boundary',name:'Cell boundary',color:'#e6a8cf',description:'The plasma membrane separates the cell interior from its surroundings and controls exchange and signaling.'},
 {id:'cell-nucleus',name:'Nucleus',color:'#b77bba',description:'The nucleus contains most of the cell’s DNA. Its envelope and pores regulate exchange with the cytoplasm.'},
 {id:'cell-endomembrane',name:'Endomembrane system',color:'#78beb6',description:'The endoplasmic reticulum and Golgi apparatus help make, modify, and transport cell products.'},
 {id:'cell-mitochondria',name:'Mitochondria',color:'#e49a54',description:'Mitochondria convert energy from nutrients into ATP through cellular respiration.'},
 {id:'cell-cytoskeleton',name:'Cytoskeleton',color:'#8bc5a8',description:'Microfilaments help maintain cell shape and support movement and intracellular organization.'},
 {id:'cell-vesicles',name:'Vesicles and lysosomes',color:'#9b9dd6',description:'Membrane-bound compartments sort, transport, and break down material inside the cell.'},
 {id:'cell-protein-synthesis',name:'Protein synthesis',color:'#d1a056',description:'Ribosomes read RNA to assemble proteins. The model also shows representative RNA and protein structures.'},
 {id:'cell-division',name:'Cell division',color:'#d9c45f',description:'Centrioles help organize microtubules in the centrosome during cell division.'},
];
export interface Part {suppressed?:boolean;tissue?:string;regions?:string[];depth?:number;sourceId?:string;sourceOffset?:number[];groups?:string[];provenance?:{label:string;url:string;detail:string};material?:string;uvs?:number;colors?:number;id:string;name:string;conceptId:string;system:SystemId;chunk:number;positions:number;normals:number;indices:number;vertexCount:number;indexCount:number;bounds:[number[],number[]]}
export interface Concept {id:string;name:string;elements:string[]}
export interface Atlas {quality?:'full'|'optimized';version:string;sex?:'male'|'female';source?:string;scope?:string;parts:Part[];concepts:Concept[];chunks:{url:string;bytes:number;gzip?:string;gzipBytes?:number}[];materials?:Record<string,{color:number[];map?:string;normalMap?:string;normalScale?:number;roughness?:number;metalness?:number;opacity?:number;vertexColors?:boolean}>;triangles:number}
export type View = 'three-quarter'|'front'|'back'|'side'|'right'|'superior'|'inferior';
export interface SceneState {focus?:number;contextOpacity?:number;region?:string;peel?:number;depthHidden?:string[];hidden?:string[];labels?:boolean;camera?:number[];section?:import('./section-plane').SectionState;sections?:import('./section-plane').SectionState[];activeSection?:number;skinOpacity?:number;inspectorOpen?:boolean;explode:number;visible:SystemId[];selected:string[];isolate:boolean;view:View;rotate:boolean;reset:number}
export const DEFAULT_VISIBLE:SystemId[] = ['cardiac','sensory','skeletal','muscular','arterial','venous','nervous','respiratory','digestive','urinary','lymphatic','endocrine','reproductive','connective'];
export const CELL_VISIBLE:SystemId[] = ['cell-boundary','cell-nucleus','cell-endomembrane','cell-mitochondria','cell-cytoskeleton','cell-vesicles','cell-protein-synthesis','cell-division'];
export const isSurfaceSystem=(system:SystemId)=>system==='integumentary'||system==='cell-boundary';
// Source meshes can be split by rendering material (for example Bone-7 or
// Brain-Inner). Those labels identify surfaces, not anatomical structures.
export const structureName=(name:string)=>name.replace(/\s*·\s*[^·]+$/u,'').trim();
export const surfaceRole=(name:string)=>{
 const role=/\s*·\s*([^·]+)$/u.exec(name)?.[1];
 return role&&!/^(?:Bone(?:-\d+)?|Brain(?:-Inner)?|Skin-.*|Lung-\d+|Suture-\d+)$/i.test(role)?role:null;
};
export const EXPLANATIONS:Record<string,string> = {
 'nuclear envelope':'A double membrane surrounding the nucleus. Its pores regulate the movement of molecules between nucleus and cytoplasm.',
 'nuclear pores':'Protein complexes in the nuclear envelope that control transport of RNA and proteins between nucleus and cytoplasm.',
 'nucleolus':'A region inside the nucleus where ribosomal RNA is made and ribosome subunits begin to assemble.',
 'heterochromatin':'Densely packed chromatin in the nucleus. It contains DNA and associated proteins and is generally less accessible for transcription.',
 'microfilaments':'Thin actin filaments of the cytoskeleton that help shape the cell and support movement.',
 'microfilaments (set 2)':'An additional modeled set of actin microfilaments that supports cell shape and movement.',
 'golgi apparatus':'Stacked membrane compartments that modify, sort, and package proteins and lipids for delivery.',
 'rough endoplasmic reticulum':'A membrane network studded with ribosomes where many proteins are made and begin folding.',
 'ribosomes':'Molecular machines that translate messenger RNA into proteins.',
 'proteins':'Representations of proteins, molecules built from amino acids that carry out many cellular functions.',
 'mitochondrial outer membranes':'The outer boundaries of the modeled mitochondria, organelles that help generate cellular ATP.',
 'mitochondrial interior':'Internal mitochondrial structures involved in energy conversion and ATP production.',
 'endosome contents':'Material carried inside an endosome, a compartment that sorts molecules taken into the cell.',
 'endosome outer membrane':'The membrane enclosing a modeled endosome.',
 'endosome inner membrane':'An internal membrane shown within the modeled endosome.',
 'centrioles':'Cylindrical microtubule structures within the centrosome that help organize the cell division machinery.',
 'plasma membrane':'The cell’s outer boundary, made mainly of a lipid bilayer with embedded proteins.',
 'lysosomal membrane':'The boundary of a lysosome, which contains enzymes that break down cellular material.',
 'lysosomal contents':'Material inside a lysosome, a compartment used for breakdown and recycling.',
 'rna':'A nucleic acid that helps carry and use genetic information, including instructions for making proteins.',
 'sternum':'The sternum, or breastbone, is the flat bone at the front of the chest. Its manubrium, body, and xiphoid process form the front of the rib cage and connect to the ribs through costal cartilage.',
 'body of sternum':'The body is the long middle portion of the sternum, between the manubrium above and the xiphoid process below. Its sides meet the costal cartilages of the ribs.',
 'manubrium of sternum':'The manubrium is the broad upper portion of the sternum. It meets the clavicles and the body of the sternum at the sternal angle.',
 'xiphoid process':'The xiphoid process is the small lower portion of the sternum, below its body. It provides attachment for the diaphragm and abdominal muscles.',
 'ureter':'A muscular tube that carries urine from a kidney to the urinary bladder. One ureter descends on each side of the abdomen and pelvis.',
 'ureters':'The two ureters are muscular tubes that carry urine from the kidneys to the urinary bladder.',
 'femoral nerve':'The femoral nerve arises from the lumbar plexus and passes into the front of the thigh. It supplies muscles that extend the knee and sensation to parts of the thigh and leg.',
 'sciatic nerve':'The sciatic nerve leaves the pelvis and runs down the back of the thigh. Its branches provide motor and sensory innervation to much of the lower limb.',
 'femoral artery':'The femoral artery continues from the external iliac artery into the thigh. It supplies the lower limb and gives rise to the deep femoral artery.',
 'femoral vein':'The femoral vein drains blood from the lower limb. It accompanies the femoral artery in the thigh and continues as the external iliac vein above the inguinal ligament.',
 'uterus':'A hollow muscular organ in the pelvis. Its lining changes through the menstrual cycle and can support implantation and development during pregnancy.',
 'vagina':'A muscular canal connecting the cervix of the uterus to the outside of the body. It provides a passage for menstrual flow and forms part of the birth canal.',
 'ovary':'An organ that contains developing oocytes and produces hormones including estrogen and progesterone.',
 'heart':'A muscular pump in the chest. Its right side sends blood to the lungs; its left side sends blood through the systemic circulation.',
 'liver':'A large organ beneath the right side of the diaphragm. It processes absorbed nutrients, produces bile, and synthesizes many proteins carried in the blood.',
 'brain':'The central organ of the nervous system. Its interconnected regions support perception, movement, memory, language, and the regulation of bodily functions.',
 'stomach':'A muscular chamber between the esophagus and small intestine. It stores and mixes food with acid and enzymes before releasing it into the duodenum.',
 'spleen':'A lymphoid organ in the upper left abdomen. It filters blood, removes aging blood cells, and participates in immune responses.',
 'pancreas':'An abdominal organ with digestive and endocrine roles. It supplies enzymes to the small intestine and releases hormones including insulin and glucagon.',
 'urinary bladder':'A muscular reservoir in the pelvis that stores urine arriving from the kidneys through the ureters.',
 'trachea':'The main airway connecting the larynx to the bronchi. Its cartilage supports keep the airway open during breathing.',
 'diaphragm':'A broad muscle separating the chest and abdomen. When it contracts, it increases chest volume and helps draw air into the lungs.',
};
export function explanation(name:string){return EXPLANATIONS[structureName(name).replace(/\s*\((?:left|right)\)$/i,'').toLowerCase()] ?? '';}
