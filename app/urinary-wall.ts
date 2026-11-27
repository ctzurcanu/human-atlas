import type {Atlas,Part} from './anatomy';
import type {HierarchyChoice} from './hierarchy-choice';
import type {WallPeelStep} from './bowel-wall-peeling';

type Chart={name:string;partIds:string[];bankSha256:string};
type OpeningBoundary={name:string;representation:string;partIds:string[];points:number[][];samplePartIds:string[];faceIds:number[];barycentric:number[][];sourcePlanSha256:string;maximumSourceOffsetM:number;sourceCurveControls:number;samplesPerSourceEdge:number};
type PointReference={representation:'point-landmark';partId:string;point:[number,number,number];face:number;barycentric:[number,number,number];sourceBankSha256:string;sourcePlanSha256:string};
type ApexTransfer={schema:'bladder-apex-similarity/v1';sourcePartId:string;sourceCatalogueSha256:string;sourceBankSha256:string;sourcePoint:[number,number,number];matrix:number[][];transferPlanSha256:string;maximumSourceOffsetM:number;targetPartId:string};
type OpeningConfiguration={state:'voiding';openingTerm:'A08.3.01.027M'|'A08.3.01.028F';sourcePlanSha256:string;sourceBankSha256:string};
type TermGroup={term:string;conceptId:string;partIds:string[];openingBoundaries?:OpeningBoundary[];pointReference?:PointReference;openingConfiguration?:OpeningConfiguration};
type Interface={first:string;second:string;triangles:number;coordinate_normal_signature:string;actual_decoded_positions_and_opposed_words_match:boolean};
type ReferenceContact={first:string;second:string;firstFace:number;secondFace:number;triangles:1;representation:'open-reference-sheet-contact';coordinate_normal_signature:string;actual_decoded_positions_and_same_words_match:boolean};
type ApicalAttachment={partId:string;wallPartId:string;sourcePlanSha256:string;referenceContact:ReferenceContact};
type NativeContact={first:string;second:string;firstFaces:number[];secondFaces:number[];triangles:1;representation:'fixed-native-surface-contact';coordinate_normal_signature:string;actual_decoded_positions_and_opposed_words_match:true};
type MuscleExtensions={sourcePlanSha256:string;posteriorPartId:string;pubicAttachments:{partId:string;chart:string;sourceContact:NativeContact}[]};
type Marker={schema:string;sex:string;partIds:string[];charts:Chart[];termGroups:TermGroup[];externalInterfaces:Interface[];openingPlanSha256?:string;apicalAttachment?:ApicalAttachment;apexTransfer?:ApexTransfer;muscleExtensions?:MuscleExtensions};
type UrinaryPart=Part&{urinaryOrgan?:'bladder'|'urethra'|'ureter';urinaryRegion?:string;urinarySourcePartId?:string};
type Registry={marker:Marker;parts:Part[];byId:Map<string,Part>;chartFor:Map<string,string>};
const cache=new WeakMap<Atlas,Registry|null>();
const same=(a:readonly string[],b:readonly string[])=>a.length===b.length&&new Set(a).size===a.length&&a.every(id=>b.includes(id));
const sha=(s:unknown)=>typeof s==='string'&&/^[a-f0-9]{64}$/.test(s);

/** Cross-bank joins are permitted only by the independently decoded seam
 * ledger, never by proximity, a shared anatomical parent or a part's name. */
export function urinaryWallRegistry(atlas:Atlas):Registry|null{
 if(cache.has(atlas))return cache.get(atlas)!;
 const reject=()=>{cache.set(atlas,null);return null;};
 const m=(atlas as Atlas&{finiteUrinaryIntegration?:Marker}).finiteUrinaryIntegration;
 if(!m||!['urinary-local-coats-v1','urinary-local-coats-v2','urinary-local-coats-v3','urinary-local-coats-v4'].includes(m.schema)||m.sex!==atlas.sex||!['male','female'].includes(m.sex)
  ||!Array.isArray(m.partIds)||!Array.isArray(m.charts)||!Array.isArray(m.termGroups)||!Array.isArray(m.externalInterfaces))return reject();
 if(m.partIds.some(id=>typeof id!=='string')||m.charts.some(g=>!g||typeof g.name!=='string'||!Array.isArray(g.partIds)))return reject();
 const featureRoles:Record<string,string[]>={'A08.3.01.025':['mucosa-interureteric-crest'],'A08.3.01.006':['mucosa-uvula'],
  'A08.3.01.026':['mucosa-ureteric-orifice-left','mucosa-ureteric-orifice-right'],
  'A08.3.01.027M':['mucosa-internal-urethral-orifice'],'A08.3.01.028F':['mucosa-internal-urethral-orifice']};
 const features=m.schema!=='urinary-local-coats-v1'?m.termGroups.filter(g=>g&&g.term in featureRoles):[];
 if(m.schema!=='urinary-local-coats-v1'&&(!features.some(g=>g.term==='A08.3.01.025')
  ||m.sex==='male'&&!features.some(g=>g.term==='A08.3.01.006')
  ||features.some(g=>!Array.isArray(g.partIds)||!g.partIds.length)))return reject();
 const hasOpenings=['urinary-local-coats-v3','urinary-local-coats-v4'].includes(m.schema),regionalUreters=m.schema==='urinary-local-coats-v4';
 if(hasOpenings&&(!sha(m.openingPlanSha256)
  ||!features.some(g=>g.term==='A08.3.01.026')||!features.some(g=>g.term===(m.sex==='male'?'A08.3.01.027M':'A08.3.01.028F'))))return reject();
 const extra=features.flatMap(g=>g.partIds);if(new Set(extra).size!==extra.length)return reject();
 const attachment=m.apicalAttachment;
 const apexTransfer=m.apexTransfer;
 if(apexTransfer){
  const t=apexTransfer,mat=t.matrix;
  if(m.sex!=='female'||t.schema!=='bladder-apex-similarity/v1'||t.sourcePartId!=='BLADDER:male:body:adventitia'||t.targetPartId!=='BLADDER:female:body:adventitia'
   ||!sha(t.sourceCatalogueSha256)||!sha(t.sourceBankSha256)||!sha(t.transferPlanSha256)
   ||!Array.isArray(t.sourcePoint)||t.sourcePoint.length!==3||!t.sourcePoint.every(Number.isFinite)
   ||!Array.isArray(mat)||mat.length!==4||mat.some(row=>!Array.isArray(row)||row.length!==4||!row.every(Number.isFinite))
   ||mat[3].some((v,i)=>v!==(i===3?1:0))||!Number.isFinite(t.maximumSourceOffsetM)||t.maximumSourceOffsetM<0||t.maximumSourceOffsetM>1e-7)return reject();
  const scale2=mat.slice(0,3).reduce((sum,row)=>sum+row[0]*row[0],0);
  if(scale2<=0)return reject();
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)if(Math.abs(mat.slice(0,3).reduce((sum,row)=>sum+row[i]*row[j],0)-(i===j?scale2:0))>scale2*1e-8)return reject();
  const determinant=mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])-mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])+mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]);
  if(determinant<=0)return reject();
 }
 const muscles=m.muscleExtensions;
 if(muscles&&(m.sex!=='male'||!sha(muscles.sourcePlanSha256)||muscles.posteriorPartId!=='BLADDER:male:neck:vesicoprostatic'
  ||!Array.isArray(muscles.pubicAttachments)||!same(muscles.pubicAttachments.map(p=>p?.chart),['pubovesical-left','pubovesical-right'])))return reject();
 if(attachment&&(!sha(attachment.sourcePlanSha256)||attachment.partId!==`BLADDER:${m.sex}:median-umbilical-ligament`
  ||attachment.wallPartId!==`BLADDER:${m.sex}:body:adventitia`))return reject();
 const expected:Record<string,number>={bladder:22+extra.length+(muscles?1:0),urethra:m.sex==='male'?9:3,left:regionalUreters?30:33,right:regionalUreters?30:33,...(m.sex==='male'?{prostate:1}:{}),...(attachment?{attachments:1}:{}),...(muscles?{'pubovesical-left':1,'pubovesical-right':1}:{}),
  ...(regionalUreters?{'ureter-left-extramural':6,'ureter-right-extramural':6,'ureter-left-intramural':3,'ureter-right-intramural':3}:{})};
 if(!same(m.charts.map(g=>g.name),Object.keys(expected))||!same(m.partIds,m.charts.flatMap(g=>g.partIds??[])))return reject();
 const byId=new Map(atlas.parts.map(p=>[p.id,p])),chartFor=new Map<string,string>(),parts:Part[]=[];
 for(const g of features)for(const id of g.partIds){const p=byId.get(id)as UrinaryPart|undefined;
  const physicalTerm=['A08.3.01.026','A08.3.01.027M','A08.3.01.028F'].includes(g.term)?'A08.3.01.023':g.term;
  if(!p||p.urinaryOrgan!=='bladder'||!featureRoles[g.term].includes(p.tissueWallLayer??'')||p.ta98Term!==physicalTerm)return reject();}
 for(const g of m.charts){
  if(!Array.isArray(g.partIds)||g.partIds.length!==expected[g.name]||!sha(g.bankSha256))return reject();
  const members=g.partIds.map(id=>byId.get(id)),first=members[0];
  if(!first)return reject();
  const bank=atlas.chunks[first.chunk]as {sha256?:string}|undefined;
  if(bank?.sha256!==g.bankSha256||members.some(p=>!p||p.suppressed||!p.pairedTissue||!p.sectionSolid||p.sourceOffset
   ||p.normalSigns===undefined||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.vertexCount!==first.vertexCount
   ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return reject();
  for(const p of members as Part[]){chartFor.set(p.id,g.name);parts.push(p);}
 }
 const ledgers=new Map<string,Interface>();
 const key=(a:string,b:string)=>[a,b].sort().join('\n');
 for(const i of m.externalInterfaces){
  if(!i||!chartFor.has(i.first)||!chartFor.has(i.second)||chartFor.get(i.first)===chartFor.get(i.second)
   ||!Number.isInteger(i.triangles)||i.triangles<=0||!sha(i.coordinate_normal_signature)
   ||i.actual_decoded_positions_and_opposed_words_match!==true||ledgers.has(key(i.first,i.second)))return reject();
  ledgers.set(key(i.first,i.second),i);
 }
 if(attachment){
  const p=byId.get(attachment.partId)as UrinaryPart|undefined,ref=attachment.referenceContact,peer=ref&&byId.get(ref.second);
  if(!p||chartFor.get(p.id)!=='attachments'||p.urinaryOrgan!=='bladder'||p.urinaryRegion!=='attachment'
   ||p.tissueWallLayer!=='median-umbilical-ligament'||p.ta98Term!=='A08.3.01.007'
   ||!ledgers.has(key(p.id,attachment.wallPartId))||!m.termGroups.some(g=>g.term==='A08.3.01.007'&&same(g.partIds,[p.id]))
   ||!ref||ref.first!==p.id||!peer||peer.suppressed||peer.pairedTissue||chartFor.has(peer.id)
   ||ref.representation!=='open-reference-sheet-contact'||ref.triangles!==1||!sha(ref.coordinate_normal_signature)
   ||ref.actual_decoded_positions_and_same_words_match!==true
   ||!Number.isInteger(ref.firstFace)||ref.firstFace<0||ref.firstFace>=p.indexCount/3
   ||!Number.isInteger(ref.secondFace)||ref.secondFace<0||ref.secondFace>=peer.indexCount/3
   ||!p.tissueDrawGroups?.some(g=>g.neighbor===ref.second&&g.count===3&&g.start===ref.firstFace*3))return reject();
 }
 const nativeContacts=new Map<string,NativeContact>();
 if(muscles){
  const posterior=byId.get(muscles.posteriorPartId)as UrinaryPart|undefined;
  if(!posterior||chartFor.get(posterior.id)!=='bladder'||posterior.tissueWallLayer!=='vesicoprostatic'||posterior.ta98Term!=='A08.3.01.021M'
   ||!m.termGroups.some(g=>g.term==='A08.3.01.021M'&&same(g.partIds,[posterior.id])))return reject();
  for(const entry of muscles.pubicAttachments){
   const p=byId.get(entry.partId)as UrinaryPart|undefined,ref=entry.sourceContact,peer=ref&&byId.get(ref.second);
   if(!p||p.id!==`BLADDER:male:neck:${entry.chart}`||chartFor.get(p.id)!==entry.chart||p.ta98Term!=='A08.3.01.020'
    ||p.urinaryOrgan!=='bladder'||p.urinaryRegion!=='attachment'||p.tissueWallLayer!==entry.chart
    ||!peer||peer.suppressed||chartFor.has(peer.id)||peer.system!=='skeletal'
    ||ref.first!==p.id||ref.second!==`ZA:Hip bone.${entry.chart.endsWith('left')?'l':'r'}:Bone-1`
    ||ref.representation!=='fixed-native-surface-contact'||ref.triangles!==1||!sha(ref.coordinate_normal_signature)
    ||ref.actual_decoded_positions_and_opposed_words_match!==true||!Array.isArray(ref.firstFaces)||!Array.isArray(ref.secondFaces)
    ||ref.firstFaces.length!==1||ref.secondFaces.length!==1||!Number.isInteger(ref.firstFaces[0])||ref.firstFaces[0]<0||ref.firstFaces[0]>=p.indexCount/3
    ||!Number.isInteger(ref.secondFaces[0])||ref.secondFaces[0]<0||ref.secondFaces[0]>=peer.indexCount/3
    ||!p.tissueDrawGroups?.some(g=>g.neighbor===ref.second&&g.count===3&&g.start===ref.firstFaces[0]*3))return reject();
   nativeContacts.set(key(p.id,ref.second),ref);
  }
  if(!m.termGroups.some(g=>g.term==='A08.3.01.020'&&same(g.partIds,muscles.pubicAttachments.map(p=>p.partId))))return reject();
 }
 const used=new Set<string>();
 for(const p of parts){
  let end=0;
  for(const g of p.tissueDrawGroups??[]){
   if(g.start!==end||!Number.isInteger(g.count)||g.count<=0||g.count%3)return reject();end+=g.count;
   if(!g.neighbor)continue;
   if(attachment&&p.id===attachment.partId&&g.neighbor===attachment.referenceContact.second){
    if(g.count!==3||g.start!==attachment.referenceContact.firstFace*3)return reject();
    continue;
   }
   const contact=nativeContacts.get(key(p.id,g.neighbor));
   if(contact){if(g.count!==3||g.start!==contact.firstFaces[0]*3)return reject();continue;}
   const peer=byId.get(g.neighbor);
   if(!peer||peer.id===p.id||!chartFor.has(peer.id)
    ||peer.tissueDrawGroups?.filter(q=>q.neighbor===p.id).reduce((n,q)=>n+q.count,0)!==g.count)return reject();
   if(chartFor.get(p.id)!==chartFor.get(peer.id)){
    const k=key(p.id,peer.id),ledger=ledgers.get(k);if(!ledger||ledger.triangles*3!==g.count)return reject();used.add(k);
   }
  }
  if(end!==p.indexCount)return reject();
 }
 if(used.size!==ledgers.size)return reject();
 const concepts=new Map(atlas.concepts.map(c=>[c.id,c]));
 for(const g of m.termGroups){
  if(!g||typeof g.term!=='string'||typeof g.conceptId!=='string')return reject();
  const c=concepts.get(g.conceptId);
  if(!/^A\d{2}\.\d\.\d{2}\.\d{3}[FM]?$/.test(g.term)||!Array.isArray(g.partIds)||!g.partIds.length
   ||!c||c.ta98Term!==g.term||!same(c.elements,g.partIds)||g.partIds.some(id=>!chartFor.has(id)))return reject();
  if(g.openingConfiguration||/^A09\.(?:4\.02|2\.03)\.00[12]$/.test(g.term)){
   const cfg=g.openingConfiguration,parent=cfg&&m.termGroups.find(q=>q.term===cfg.openingTerm);
   if(g.term!==(m.sex==='male'?'A09.4.02.002':'A09.2.03.002')||!cfg||cfg.state!=='voiding'||cfg.openingTerm!==(m.sex==='male'?'A08.3.01.027M':'A08.3.01.028F')
    ||!parent||c.partOf!==parent.conceptId||!same(g.partIds,parent.partIds)||g.openingBoundaries
    ||cfg.sourcePlanSha256!==m.openingPlanSha256||!sha(cfg.sourcePlanSha256)
    ||g.partIds.some(id=>cfg.sourceBankSha256!==m.charts.find(q=>q.name===chartFor.get(id))?.bankSha256)
    ||parent.openingBoundaries?.length!==1||parent.openingBoundaries[0].name!=='neck')return reject();
  }
  if(g.term==='A08.3.01.002'||g.pointReference){
   const ref=g.pointReference,p=ref&&byId.get(ref.partId);
   const transfer=apexTransfer&&ref?.partId===apexTransfer.targetPartId&&ref.sourcePlanSha256===apexTransfer.transferPlanSha256;
   const attached=attachment&&ref?.partId===attachment.wallPartId&&ref.sourcePlanSha256===attachment.sourcePlanSha256;
   if(g.term!=='A08.3.01.002'||!ref||ref.representation!=='point-landmark'||!attached&&!transfer
    ||!same(g.partIds,[ref.partId])||!p||ref.sourceBankSha256!==m.charts.find(c=>c.name===chartFor.get(p.id))?.bankSha256
    ||!sha(ref.sourcePlanSha256)
    ||!Array.isArray(ref.point)||ref.point.length!==3||ref.point.some((v,k)=>!Number.isFinite(v)||v<p.bounds[0][k]||v>p.bounds[1][k])
    ||!Number.isInteger(ref.face)||ref.face<0||ref.face>=p.indexCount/3
    ||!p.tissueDrawGroups?.some(q=>(attached?q.neighbor===attachment!.partId:!q.neighbor)&&ref.face>=q.start/3&&ref.face<(q.start+q.count)/3)
    ||!Array.isArray(ref.barycentric)||ref.barycentric.length!==3||ref.barycentric.some(v=>!Number.isFinite(v)||v<0||v>1)
    ||Math.abs(ref.barycentric.reduce((a,b)=>a+b,0)-1)>1e-12)return reject();
   if(transfer&&apexTransfer){const t=apexTransfer,mapped=t.matrix.slice(0,3).map(row=>row[3]+t.sourcePoint.reduce((sum,v,i)=>sum+v*row[i],0));if(Math.hypot(...mapped.map((v,i)=>v-ref.point[i]))>t.maximumSourceOffsetM+1e-12)return reject();}
  }
 }
 if(regionalUreters){
  const coats=['adventitia','muscularis','mucosa'],regions=['abdominal','pelvic','intramural'];
  for(const side of ['left','right'])for(const region of regions)for(const coat of coats){
   const id=`URETER:${m.sex}:${side}:${region}:${coat}`,p=byId.get(id)as UrinaryPart|undefined;
   if(!p||p.urinaryOrgan!=='ureter'||p.urinaryRegion!==`${side}:${region}`||p.tissueWallLayer!==coat||p.ta98Term!==`A08.2.01.00${coats.indexOf(coat)+5}`
    ||chartFor.get(id)!==`ureter-${side}-${region==='intramural'?'intramural':'extramural'}`)return reject();
  }
  const ureters=parts.filter(p=>(p as UrinaryPart).urinaryOrgan==='ureter');
  const groups:Record<string,Part[]>={'A08.2.01.001':ureters,
   ...Object.fromEntries(regions.map((region,i)=>[`A08.2.01.00${i+2}`,ureters.filter(p=>p.id.includes(`:${region}:`))])),
   ...Object.fromEntries(coats.map((coat,i)=>[`A08.2.01.00${i+5}`,ureters.filter(p=>p.tissueWallLayer===coat)]))};
  for(const [term,members]of Object.entries(groups))if(!m.termGroups.some(g=>g.term===term&&same(g.partIds,members.map(p=>p.id))))return reject();
 }
 if(hasOpenings)for(const term of ['A08.3.01.026',m.sex==='male'?'A08.3.01.027M':'A08.3.01.028F']){
  const g=m.termGroups.find(g=>g.term===term)!,names=term==='A08.3.01.026'?['left','right']:['neck'];
  if(!Array.isArray(g.openingBoundaries)||!same(g.openingBoundaries.map(b=>b.name),names))return reject();
  for(const b of g.openingBoundaries){
   const role=b.name==='neck'?'mucosa-internal-urethral-orifice':'mucosa-ureteric-orifice-'+b.name;
   if(b.representation!=='wall-rim'||!Array.isArray(b.partIds)||!same(b.partIds,g.partIds.filter(id=>byId.get(id)?.tissueWallLayer===role))
    ||!b.partIds.length||b.sourcePlanSha256!==m.openingPlanSha256||!Number.isFinite(b.maximumSourceOffsetM)||b.maximumSourceOffsetM<0||b.maximumSourceOffsetM>1e-7
    ||!Number.isInteger(b.sourceCurveControls)||b.sourceCurveControls<3||b.samplesPerSourceEdge!==4
    ||!Array.isArray(b.points)||b.points.length!==b.sourceCurveControls*4
    ||!Array.isArray(b.samplePartIds)||!Array.isArray(b.faceIds)||!Array.isArray(b.barycentric)
    ||[b.samplePartIds,b.faceIds,b.barycentric].some(a=>a.length!==b.points.length))return reject();
   for(let i=0;i<b.points.length;i++){
    const p=byId.get(b.samplePartIds[i]),q=b.points[i],weights=b.barycentric[i],face=b.faceIds[i];
    if(!p||!b.partIds.includes(p.id)||!Array.isArray(q)||q.length!==3||q.some((v,a)=>!Number.isFinite(v)||v<p.bounds[0][a]-1e-7||v>p.bounds[1][a]+1e-7)
     ||!Number.isInteger(face)||face<0||face>=p.indexCount/3
     ||!p.tissueDrawGroups?.some(g=>!g.neighbor&&face>=g.start/3&&face<(g.start+g.count)/3)
     ||!Array.isArray(weights)||weights.length!==3||weights.some(w=>!Number.isFinite(w)||w< -1e-6||w>1+1e-6)||Math.abs(weights.reduce((a,b)=>a+b,0)-1)>1e-6)return reject();
   }
  }
 }
 const value={marker:m,parts,byId,chartFor};cache.set(atlas,value);return value;
}

export function urinaryExternalPeer(atlas:Atlas,part:Part,peerId:string):Part|undefined{
 const r=urinaryWallRegistry(atlas);if(!r||!r.chartFor.has(part.id)||!r.chartFor.has(peerId))return;
 return r.marker.externalInterfaces.some(i=>i.first===part.id&&i.second===peerId||i.second===part.id&&i.first===peerId)?r.byId.get(peerId):undefined;
}

export function urinaryTermParts(atlas:Atlas):Map<string,Part[]>{
 const r=urinaryWallRegistry(atlas),result=new Map<string,Part[]>();
 if(r)for(const g of r.marker.termGroups)result.set(g.term,g.partIds.map(id=>r.byId.get(id)!));
 return result;
}

/** Immaterial landmarks use the ordinary label, attached to literal tissue.
 * The host tissue is context for the point, not a second anatomical solid. */
export function urinaryPointReference(atlas:Atlas,choice:HierarchyChoice):{partId:string;point:[number,number,number]}|undefined{
 const r=urinaryWallRegistry(atlas);if(!r)return;
 const term=choice.ta98Term??choice.terminology?.ta98;
 const group=r.marker.termGroups.find(g=>g.pointReference&&(g.conceptId===choice.id||g.term===term));
 if(!group?.pointReference||!same(choice.elements,group.partIds))return;
 return {partId:group.pointReference.partId,point:[...group.pointReference.point]};
}

/** A configuration reuses its real clinical rim and the ordinary callout.
 * The deeper histological tissue join is not the clinical opening. */
export function urinaryLabelReference(atlas:Atlas,choice:HierarchyChoice):{partId:string;point:[number,number,number]}|undefined{
 const point=urinaryPointReference(atlas,choice);if(point)return point;
 const r=urinaryWallRegistry(atlas);if(!r)return;
 const term=choice.ta98Term??choice.terminology?.ta98;
 const group=r.marker.termGroups.find(g=>g.openingConfiguration&&(g.conceptId===choice.id||g.term===term));
 if(!group?.openingConfiguration||!same(choice.elements,group.partIds))return;
 const boundary=r.marker.termGroups.find(g=>g.term===group.openingConfiguration!.openingTerm)?.openingBoundaries?.[0];
 if(!boundary?.points.length)return;
 const index=boundary.points.reduce((best,q,i)=>q[2]>boundary.points[best][2]?i:best,0);
 return {partId:boundary.samplePartIds[index],point:[...boundary.points[index]]as[number,number,number]};
}

/** Whole-organ supply belongs only to the complete retained replacement. */
export function urinaryWholeSource(selected:Part[],parts:Part[]):Part|undefined{
 if(!selected.length)return;
 const attachments=selected.filter(p=>['median-umbilical-ligament','pubovesical-left','pubovesical-right'].includes(p.tissueWallLayer??'')),wall=selected.filter(p=>!attachments.includes(p));
 const first=wall[0]as UrinaryPart|undefined;if(!first)return;
 const id=first.urinarySourcePartId,organ=first.urinaryOrgan;
 if(!id||!organ)return;
 const source=parts.find(p=>p.id===id);
 if(!source?.suppressed||source.system!=='urinary'||!same(source.replacementParts??[],selected.map(p=>p.id)))return;
 if(organ==='bladder'&&selected.length<22||organ==='urethra'&&![3,9].includes(selected.length))return;
 if(organ==='ureter'){
  const prefix=first.id.match(/^URETER:(male|female):(left|right):/u)?.[0];
  if(!prefix||selected.length!==9||!same(selected.map(p=>p.id),['abdominal','pelvic','intramural'].flatMap(region=>['adventitia','muscularis','mucosa'].map(coat=>`${prefix}${region}:${coat}`))))return;
 }
 if((selected as UrinaryPart[]).some(p=>p.suppressed||p.urinarySourcePartId!==id||p.urinaryOrgan!==organ
  ||p.system!=='urinary'||!p.pairedTissue||!p.sectionSolid||p.normalSigns===undefined||p.sourceOffset
  ||p.tissueWallAssembly!==first.tissueWallAssembly))return;
 if(new Set(attachments.map(p=>p.tissueWallLayer)).size!==attachments.length||attachments.some(p=>organ!=='bladder'
  ||p.ta98Term!==(p.tissueWallLayer==='median-umbilical-ligament'?'A08.3.01.007':'A08.3.01.020')
  ||(p as UrinaryPart).urinaryRegion!=='attachment'||!p.tissueDrawGroups?.some(g=>g.neighbor&&wall.some(w=>w.id===g.neighbor))))return;
 if(organ!=='ureter'&&wall.some(p=>p.chunk!==first.chunk||p.positions!==first.positions
  ||p.normals!==first.normals||p.vertexCount!==first.vertexCount||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)))return;
 return source;
}

export function urinaryTissuePeelSteps(atlas:Atlas,choice:HierarchyChoice):WallPeelStep[]|null{
 const r=urinaryWallRegistry(atlas);if(!r||choice.id==='selection-set'||!choice.elements.length)return null;
 const selected=choice.elements.map(id=>r.byId.get(id)as UrinaryPart|undefined);
 if(selected.some(p=>!p||!r.chartFor.has(p.id)||!p.urinaryOrgan))return null;
 const parts=selected as UrinaryPart[],organs=new Set(parts.map(p=>p.urinaryOrgan));if(organs.size!==1)return null;
 const ids=new Set(choice.elements),regions=new Set(parts.map(p=>p.urinaryRegion));
 for(const region of regions){
  const peers=r.parts.filter(p=>(p as UrinaryPart).urinaryOrgan===parts[0].urinaryOrgan&&(p as UrinaryPart).urinaryRegion===region);
  if(!peers.every(p=>ids.has(p.id)))return null;
 }
 const muscle=parts[0].urinaryOrgan==='bladder'?
  [{name:'Outer muscular tissue',test:(r:string)=>['detrusor-unstratified','neck-external-longitudinal','vesicoprostatic','pubovesical-left','pubovesical-right'].includes(r)},
   {name:'Circular neck and deep trigonal muscle',test:(r:string)=>['neck-circular','trigone-deep'].includes(r)},
   {name:'Inner longitudinal and superficial trigonal muscle',test:(r:string)=>['neck-internal-longitudinal','trigone-superficial'].includes(r)}]
  :[{name:'Muscular coat',test:(r:string)=>r==='muscularis'}];
 const steps=[{name:'Median umbilical ligament',test:(r:string)=>r==='median-umbilical-ligament'},
  {name:'External covering',test:(r:string)=>['serosa','adventitia'].includes(r)},
  {name:'Spongy layer',test:(r:string)=>r==='spongy'},
  {name:'Subserosa',test:(r:string)=>r==='subserosa'},...muscle,
  {name:'Submucosa',test:(r:string)=>r==='submucosa'},{name:'Mucosa',test:(r:string)=>r.startsWith('mucosa')}];
 return steps.map(s=>({name:s.name,ids:parts.filter(p=>s.test(p.tissueWallLayer??'')).map(p=>p.id)})).filter(s=>s.ids.length);
}
