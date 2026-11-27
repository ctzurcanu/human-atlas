import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';
import {organInstanceScale,organOrientation,organOrigin} from './organ-frame';

const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const nameOf=(part:Part)=>{
 const region=/^TA98HEART:(male|female):myocardium:(right-atrium|right-ventricle|left-atrium|left-ventricle)$/.exec(part.id);
 return region&&part.tissueWallAssembly===`heart:${region[1]}`?region[2].replaceAll('-',' '):part.name.toLowerCase().trim().replace(/^root of (?=(?:left|right) (?:inferior|superior) pulmonary vein$)/,'');
};
const chambers={
 leftAtrium:['left atrium','atrium (left)','cardiac atrium (left)'],
 rightAtrium:['right atrium','atrium (right)','cardiac atrium (right)'],
 leftVentricle:['left ventricle','heart left ventricle','ventricle (left)'],
 rightVentricle:['right ventricle','heart right ventricle','ventricle (right)'],
};
const valves=[
 {name:'Mitral',pattern:/^(?:mitral valve|bicuspid atrioventricular valve(?: \(left\))?|left atrioventricular valve|(?:anterior|posterior) leaflet of (?:mitral|left atrioventricular) valve)$/,up:chambers.leftAtrium,down:chambers.leftVentricle},
 {name:'Tricuspid',pattern:/^(?:tricuspid(?: atrioventricular)? valve|right atrioventricular valve|(?:anterior|posterior|inferior|septal) leaflet of (?:tricuspid|right atrioventricular) valve)$/,up:chambers.rightAtrium,down:chambers.rightVentricle},
 {name:'Aortic',pattern:/^(?:aortic(?: semilunar)? valve|(?:right coronary|left coronary|non-coronary) leaflet|(?:anterior|left posterior|right posterior) cusp of aortic valve)$/,up:chambers.leftVentricle,down:['ascending aorta']},
 {name:'Pulmonary',pattern:/^(?:pulmonary(?: semilunar)? valve|(?:anterior|right|left) semilunar leaflet of pulmonary valve|(?:right anterior|posterior|left anterior) cusp of pulmonary valve)$/,up:chambers.rightVentricle,down:['pulmonary trunk']},
];
type ChordalPartition={schema:string;role:string;leafletPartId:string;chordalPartId:string;papillaryPartId:string;sharedLeafletInterface:boolean;sharedPapillaryInterface?:boolean;papillaryInterfaceVerified:boolean};
const chordalPartition=(part:Part)=>(part as Part&{cardiacChordalPartition?:ChordalPartition}).cardiacChordalPartition;
const apparatusOwnership=(part:Part)=>(part as Part&{cardiacApparatusOwnership?:{schema:string;bankSha256:string}}).cardiacApparatusOwnership;
type CoupledAttachment={schema:string;bankSha256:string;interfaces:{peer:string;kind:string;continuous:boolean}[]};
const coupledAttachment=(part:Part)=>(part as Part&{cardiacCoupledAttachment?:CoupledAttachment}).cardiacCoupledAttachment;
// Reciprocal charts in one canonical bank remain coincident only while their
// actual placements match. A later heart-only fit can separate a fixed cuff.
const sameCoupledPose=(first:Part,second:Part)=>organOrigin(first).distanceTo(organOrigin(second))<=1e-9&&Math.abs(organInstanceScale(first)-organInstanceScale(second))<=1e-9&&1-Math.abs(organOrientation(first).dot(organOrientation(second)))<=1e-12;
export function cardiacValveDescription(part:Part):{text:string;url:string}|undefined{
 if(part.system!=='cardiac')return;
 const valve=valves.find(v=>v.pattern.test(nameOf(part)));if(!valve)return;
 const subject=/\b(?:leaflet|cusp)\b/.test(nameOf(part))?`This is a leaflet of the ${valve.name.toLowerCase()} valve.`:`The ${valve.name.toLowerCase()} valve controls flow from the ${valve.up[0]} to the ${valve.down[0]}.`;
 const detail=/\b(?:leaflet|cusp)\b/.test(nameOf(part))?` The valve controls flow from the ${valve.up[0]} to the ${valve.down[0]}.`:'';
 return {text:`${subject}${detail} Its leaflets open and close to prevent backward blood flow.`,url:'https://www.nhlbi.nih.gov/health/heart/blood-flow'};
}

export function cardiacSourceLimitation(part:Part):string|undefined{
 const chordal=chordalPartition(part);
 if(chordal?.schema==='native-chordal-partition/v1')return chordal.papillaryInterfaceVerified&&chordal.sharedPapillaryInterface?'Native leaflet and branching chordal source tissue with reciprocal finite interfaces at the leaflet neck and papillary insertion. Computational neck/basal boundaries and complete fine chordal distribution remain unverified.':'Native leaflet and branching chordal source tissue with a shared finite neck interface. The strand-neck ownership boundary is illustrative; the papillary source interfaces remain unjoined and unverified.';
 if(part.conceptId==='LOCAL:female:ventricles')return 'Combined ventricular envelope. Separate left/right chamber meshes, myocardial/endocardial layers and complete lumen anatomy are not established by this source.';
 const identity=(part as Part&{cardiacIdentityEvidence?:{named_source?:string}}).cardiacIdentityEvidence;
 const resolvedPosterior=(part.ta98Term==='A12.1.02.005'&&identity?.named_source==='FJ2433')||(part.ta98Term==='A12.1.02.018'&&identity?.named_source==='FJ2430');
 if(!resolvedPosterior&&['ZA:Inferior leaflet of right atrioventricular valve','ZA:Inferior papillary muscle of left ventricle','ZA:Inferior papillary muscle of right ventricle'].includes(part.conceptId))return 'Named inferior cardiac source tissue. Its correspondence to the TA98 posterior subdivision remains unresolved; hierarchy placement identifies the broader valve or papillary group only. Extent and physical attachments remain unverified.';
 if(part.system==='cardiac'&&nativePartFrame(part)&&valves.some(v=>v.pattern.test(nameOf(part))))return 'Named native valve tissue. Separate wall layers, nodules, lunules, commissures, annular attachments and lumen continuity remain unverified; the source label does not certify those finer features.';
}

/** Named valvular relationships, never an assertion of certified mesh contact.
 * NIH NHLBI Heart valves; OpenStax A&P 2e 19.1. Donor HRA organs do not
 * acquire a native frame just because their labels describe the same heart. */
export function cardiacRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),frames=new Map<string,Part[]>();
 for(const part of parts){const frame=nativePartFrame(part);if(!frame||part.suppressed||!['cardiac','arterial','venous'].includes(part.system))continue;const list=frames.get(frame)??[];list.push(part);frames.set(frame,list);}
 const add=(source:Part,target:Part,note:string,kind:'adjacent'|'continuous'='adjacent')=>{
  if(source.id===target.id)return;
  for(const [from,to] of [[source,target],[target,source]]){const list=graph.get(from.id)??[];if(!list.some(r=>r.target.id===to.id&&r.kind===kind))list.push({kind,target:to,note});graph.set(from.id,list);}
 };
 for(const peers of frames.values()){
  for(const source of peers){
   const p=chordalPartition(source);if(p?.schema!=='native-chordal-partition/v1'||p.role!=='chordae'||!source.pairedTissue||!source.sectionSolid||!p.sharedLeafletInterface)continue;
   const leaflet=peers.find(q=>q.id===p.leafletPartId),pm=peers.find(q=>q.id===p.papillaryPartId);
   if(leaflet?.pairedTissue&&leaflet.sectionSolid&&source.tissueDrawGroups?.some(g=>g.neighbor===leaflet.id)&&leaflet.tissueDrawGroups?.some(g=>g.neighbor===source.id))add(source,leaflet,'Native chordal strands and leaflet tissue share the same finite neck surface, with opposed ownership. The computational tissue boundary is illustrative.','continuous');
   if(pm&&/papillary/i.test(pm.name)){
    const joined=p.papillaryInterfaceVerified&&p.sharedPapillaryInterface&&pm.pairedTissue&&pm.sectionSolid&&pm.normalSigns!==undefined&&source.tissueDrawGroups?.some(g=>g.neighbor===pm.id)&&pm.tissueDrawGroups?.some(g=>g.neighbor===source.id);
    add(source,pm,joined?'The native source insertion has a reciprocal finite chordal–papillary interface. Basal ownership and complete finer chordal anatomy remain provisional.':'The native chordal fan meets this papillary source region. Its shared tissue interface has not been established.',joined?'continuous':'adjacent');
   }
  }
  for(const source of peers){
   const ownership=apparatusOwnership(source);if(ownership?.schema!=='native-apparatus-ownership/v1'||!source.pairedTissue||!source.sectionSolid)continue;
   for(const g of source.tissueDrawGroups??[]){
    const target=peers.find(p=>p.id===g.neighbor);if(!target?.pairedTissue||!target.sectionSolid||apparatusOwnership(target)?.bankSha256!==ownership.bankSha256||!target.tissueDrawGroups?.some(h=>h.neighbor===source.id)||graph.get(source.id)?.some(r=>r.kind==='continuous'&&r.target.id===target.id))continue;
    const leaflets=/leaflet|cusp/i.test(source.name)&&/leaflet|cusp/i.test(target.name),papillary=/papillary/i.test(source.name)&&/papillary/i.test(target.name);
    add(source,target,leaflets?'Native leaflet coaptation/sliding contact. The common mesh interface does not fuse the leaflets.':papillary?'Common native basal muscle tissue is partitioned between these source regions; the ownership boundary remains illustrative.':'These native source tissues have a shared computational contact boundary. It does not establish a fused anatomical attachment.');
   }
  }
  for(const source of peers){
   const attachment=coupledAttachment(source);if(attachment?.schema!=='coupled-cardiac-attachment/v1'||!source.pairedTissue||!source.sectionSolid||source.normalSigns===undefined)continue;
   for(const contract of attachment.interfaces??[]){
    if(!contract.continuous||!['papillary-wall-insertion','leaflet-host-attachment','vessel-host-attachment','annular-leaflet-insertion','annular-myocardial-insertion'].includes(contract.kind))continue;
    const target=peers.find(p=>p.id===contract.peer),other=target&&coupledAttachment(target);
    if(!target?.pairedTissue||!target.sectionSolid||target.normalSigns===undefined||other?.bankSha256!==attachment.bankSha256||!other.interfaces?.some(r=>r.peer===source.id&&r.kind===contract.kind&&r.continuous)||!source.tissueDrawGroups?.some(g=>g.neighbor===target.id)||!target.tissueDrawGroups?.some(g=>g.neighbor===source.id))continue;
    const note=contract.kind.startsWith('annular-')?'Finite annular connective tissue shares this exact boundary with the leaflet hinge or myocardial host. The support region replaces local myocardial volume; its histological extent and discontinuities remain provisional.':contract.kind==='papillary-wall-insertion'?'Native papillary base and ventricular wall have a reciprocal finite tissue interface. Regional muscle ownership remains illustrative.':contract.kind==='vessel-host-attachment'?'The finite caval wall and atrial tissue have reciprocal contact charts. The exported entry course and eight offset tracks are unobstructed; wall thickness and the distal recipient join remain provisional.':'Native leaflet tissue attaches through a reciprocal finite host interface. The hinge footprint is retained; finer fibrous-ring histology remains provisional.';
    const coincident=sameCoupledPose(source,target);
    add(source,target,coincident?note:'The canonical source charts were reciprocal, but their current placements differ. Their shared interface requires reconnection and is unverified.',coincident?'continuous':'adjacent');
   }
  }
  const named=(names:string[])=>peers.filter(p=>names.includes(nameOf(p))&&(p.system==='cardiac'||names.includes('ascending aorta')||names.includes('pulmonary trunk')));
  for(const valve of valves)for(const source of peers.filter(p=>p.system==='cardiac'&&valve.pattern.test(nameOf(p)))){
   const note=`${valve.name} valve region separates ${valve.up[0]} from ${valve.down[0]}. Named chamber/outflow adjacency; source leaflet extent, annular contact and lumen continuity remain unverified.`;
   for(const target of [...named(valve.up),...named(valve.down)])add(source,target,note);
   if(valve.name==='Aortic'&&!named(['ascending aorta']).length)for(const target of peers.filter(p=>nameOf(p)==='aortic arch'))add(source,target,`${note} The recipient retains its native aortic arch; the intervening root join remains provisional.`);
   // A single ventricular envelope may be shown as context, but never renamed
   // or counted as an independently represented left/right ventricular chamber.
   if([...valve.up,...valve.down].some(n=>n.includes('ventricle')))for(const target of named(['ventricles']))add(source,target,`${note} This source is the combined ventricular envelope; separate left/right ventricular representation is missing.`);
  }
  for(const source of peers.filter(p=>p.system==='cardiac'&&/papillary.*(?:left|right) ventricle/.test(nameOf(p)))){
   const side=nameOf(source).includes('left ventricle')?'left':'right';
   for(const target of named(side==='left'?chambers.leftVentricle:chambers.rightVentricle))add(source,target,`Named ${side} ventricular papillary muscle projects from its ventricular wall. Chordal attachment and shared source tissue interfaces remain unverified; semilunar valves have no papillary/chordal apparatus.`);
  }
  const inflows=[
   {names:['superior vena cava','inferior vena cava','inferior vena cava (thoracic part)'],chambers:chambers.rightAtrium},
   {names:['left inferior pulmonary vein','left superior pulmonary vein','right inferior pulmonary vein','right superior pulmonary vein'],chambers:chambers.leftAtrium},
  ];
  for(const inflow of inflows)for(const source of peers.filter(p=>['venous','cardiac'].includes(p.system)&&inflow.names.includes(nameOf(p)))){
   for(const target of named(inflow.chambers)){
    const note=`Named venous return to the ${inflow.chambers[0]}. Source openings, wall/lumen continuity and physical venous–atrial contact remain unverified.`;
    for(const [from,to,kind] of [[source,target,'after'],[target,source,'before']] as const){
     const list=graph.get(from.id)??[];if(!list.some(r=>r.kind===kind&&r.target.id===to.id))list.push({kind,target:to,note});graph.set(from.id,list);
    }
   }
  }
 }
 cache.set(parts,graph);return graph;
}
