import type {Part} from './anatomy';
import type {ResolvedRelation,RelationKind} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';
// Native regional exteriors; references: NCBI Bookshelf NBK470591/NBK532968.
const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const name=(p:Part)=>p.name.toLowerCase().replace(/\s*\((?:left|right)\)$/,'').replace(/1$/,'');
export function fibularRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const previous=cache.get(parts);if(previous)return previous;
 const graph=new Map<string,ResolvedRelation[]>(),assemblies=new Map<string,Map<string,Part>>();
 for(const p of parts){
  const match=/^LOCAL:female:ta98-fibula:(A02\.5\.07\.(?:002|005|006|014)):(left|right)$/.exec(p.id);
  if(!match||p.suppressed||p.system!=='skeletal'||p.sectionAssembly!==`LOCAL:female:${match[2]==='left'?'l':'r'}_fibula`)continue;
  const group=assemblies.get(p.sectionAssembly)??new Map();group.set(match[1],p);assemblies.set(p.sectionAssembly,group);
 }
 const add=(from:Part,kind:RelationKind,to:Part,note:string)=>{const list=graph.get(from.id)??[];if(!list.some(r=>r.kind===kind&&r.target.id===to.id))list.push({kind,target:to,note});graph.set(from.id,list);};
 const pair=(a:Part,kind:RelationKind,b:Part,inverse:RelationKind,note:string)=>{add(a,kind,b,note);add(b,inverse,a,note);};
 for(const group of assemblies.values()){
  if(group.size!==4)continue;
  const head=group.get('A02.5.07.002')!,neck=group.get('A02.5.07.005')!,shaft=group.get('A02.5.07.006')!,malleolus=group.get('A02.5.07.014')!,side=head.id.endsWith(':left')?'left':'right',frame=nativePartFrame(head);
  const find=(labels:string[],systems:string[])=>parts.filter(p=>!p.suppressed&&nativePartFrame(p)===frame&&p.name.endsWith(`(${side})`)&&systems.includes(p.system)&&labels.includes(name(p)));
  const connect=(region:Part,labels:string[],systems:string[],kind:RelationKind,inverse:RelationKind,note:string)=>{for(const p of find(labels,systems))pair(region,kind,p,inverse,note);};
  connect(head,['tibia'],['skeletal'],'articulates','articulates','Proximal tibiofibular articulation at the head. Native articular footprints, cartilage thickness and source packing remain unverified.');
  connect(malleolus,['talus'],['skeletal'],'articulates','articulates','Lateral malleolar articulation with the talus. The exterior does not independently establish the talar facet or cartilage layer.');
  connect(head,['fibular collateral ligament'],['connective'],'connectedBy','connects','Fibular collateral ligament attachment at the head; source enthesis and packing remain unaccepted.');
  connect(head,['biceps femoris longus','biceps femoris short'],['muscular'],'muscles','insertion','Biceps femoris distal fibular-head attachment. Source tendon continuity and full enthesis footprint remain unverified.');
  connect(shaft,['extensor digitorum longus','extensor hallucis longus'],['muscular'],'muscles','origin','Fibular shaft origin of the named anterior-compartment muscle. Exact native origin footprint and compartment packing remain unverified.');
  connect(malleolus,['anterior talofibular ligament','posterior talofibular ligament','anterior tibiofibular ligament','posterior tibiofibular ligament'],['connective'],'connectedBy','connects','Named ankle ligament attachment at the distal fibula. Source attachment footprints and packing remain unaccepted.');
  connect(shaft,['peroneal artery','fibular artery'],['arterial'],'arterial','supplies','Fibular (peroneal) arterial nutrient and periosteal branches supply the shaft. Individual source branches and nutrient foramina remain unverified.');
  connect(head,['anterior tibial artery'],['arterial'],'arterial','supplies','Proximal fibular head/epiphyseal supply through anterior tibial branches; the named artery is not a reconstructed local branch.');
  connect(neck,['common peroneal nerve','common fibular nerve'],['nervous'],'adjacent','adjacent','Common fibular nerve course around the neck. This is a regional course relation, not fibular-bone innervation or certified mesh contact.');
  for(const [a,b] of [[head,neck],[neck,shaft],[shaft,malleolus]])pair(a,'adjacent',b,'adjacent','Adjoining regions of the native exterior partition with shared seams. Anatomical limits and independent internal volumes remain provisional.');
 }
 cache.set(parts,graph);return graph;
}
