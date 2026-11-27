import type {Atlas,Concept,Part} from './anatomy';

/** Stored native regions and open articular patches; no anatomy acceptance. */
export interface CricoidPartition {
 schema:'ta98-cricoid-source-partition-v1';
 sex:'female'|'male';
 wholeConcept:string;
 sourcePartId:string;
 partIds:string[];
 carrierSha256:string;
 constructionSha256:string;
 provisionalBoundaries:true;
 clinicalShapeAccepted:false;
 jointOrFacetAccepted:false;
 laminaRemainderIndependentSolid:false;
 laminaClosureParts:string[];
 sourceNeighborsMoved:false;
 airwayCapsEmitted:false;
 installed:boolean;
}

type CricoidPart=Part&{cricoidSourcePartId?:string;originalPartId?:string;cricoidRegion?:string;cricoidSurface?:string;defaultHidden?:boolean};
type PartialConcept=Concept&{ta98ParentTerm?:string};
export interface CricoidPartitionRegistry {
 marker:CricoidPartition;
 source:Part;
 whole:Concept;
 arch:Concept;
 lamina:Concept;
 arytenoid:Concept;
 thyroid:Concept;
 members:Part[];
 termParts:ReadonlyMap<string,Part[]>;
}

const suffixes=['arch','lamina-remainder','arytenoid-articular:left','arytenoid-articular:right','thyroid-articular:left','thyroid-articular:right'] as const;
const terms=['A06.2.03.001','A06.2.03.002','A06.2.03.003','A06.2.03.004','A06.2.03.005'] as const;
const registries=new WeakMap<Atlas,CricoidPartitionRegistry|null>();
const sameMembers=(value:unknown,expected:readonly string[]):value is string[]=>Array.isArray(value)&&value.length===expected.length
 &&new Set(value).size===value.length&&value.every(id=>typeof id==='string'&&expected.includes(id));
const offset=(value:unknown):value is number=>typeof value==='number'&&Number.isInteger(value)&&value>=0;
const digest=(value:unknown)=>typeof value==='string'&&/^[0-9a-f]{64}$/.test(value);

/** Validate literal source identity and complete stored memberships, not shape. */
export function cricoidPartitionRegistry(atlas:Atlas):CricoidPartitionRegistry|null {
 if(registries.has(atlas))return registries.get(atlas)!;
 const reject=()=>{registries.set(atlas,null);return null;},marker=atlas.cricoidPartition;
 if(!marker||(marker.sex!=='female'&&marker.sex!=='male')||atlas.sex!==marker.sex)return reject();
 const sex=marker.sex,sourceId=sex==='female'?'LOCAL:female:cricoid_cartilage:1072':'ZA:Cricoid cartilage';
 const wholeId=sex==='female'?'LOCAL:female:cricoid_cartilage':sourceId;
 const nativeSource=sex==='female'?'cricoid_cartilage':'Cricoid cartilage';
 const ids=suffixes.map(suffix=>`CRICOID:${sex}:${suffix}`),laminaIds=ids.slice(1);
 if(marker.schema!=='ta98-cricoid-source-partition-v1'
  ||marker.sourcePartId!==sourceId||marker.wholeConcept!==wholeId||!sameMembers(marker.partIds,ids)
  ||!sameMembers(marker.laminaClosureParts,laminaIds)||!digest(marker.carrierSha256)||!digest(marker.constructionSha256)
  ||marker.provisionalBoundaries!==true||marker.clinicalShapeAccepted!==false||marker.jointOrFacetAccepted!==false
  ||marker.laminaRemainderIndependentSolid!==false||marker.sourceNeighborsMoved!==false||marker.airwayCapsEmitted!==false
  ||typeof marker.installed!=='boolean')return reject();
 const uniquePart=(id:string)=>{const rows=atlas.parts.filter(p=>p.id===id);return rows.length===1?rows[0]:undefined;};
 const uniqueConcept=(id:string)=>{const rows=atlas.concepts.filter(c=>c.id===id);return rows.length===1?rows[0]:undefined;};
 const source=uniquePart(sourceId),members=ids.map(uniquePart),first=members[0];
 if(!source||source.suppressed!==true||source.sourceOffset
  ||(source.system!=='respiratory'&&!(sex==='male'&&source.system==='skeletal'))
  ||source.conceptId!==wholeId||source.sourceId!==nativeSource||!sameMembers(source.replacementParts,ids)
  ||!first||members.some(p=>!p))return reject();
 const bank=atlas.chunks[first.chunk];
 if(!bank||!offset(bank.bytes)||!first.geometryOrigin||first.geometryOrigin.length!==3||!first.geometryOrigin.every(Number.isFinite))return reject();
 for(let i=0;i<members.length;i++){
  const p=members[i] as CricoidPart,suffix=suffixes[i];
  if(p.suppressed||p.defaultHidden||atlas.defaultHidden?.includes(p.id)||p.sourceOffset||p.system!=='respiratory'
   ||p.conceptId!==p.id||p.sourceId!==nativeSource||p.cricoidSourcePartId!==sourceId||p.originalPartId!==sourceId
   ||p.cricoidRegion!==(i===0?'arch':'lamina')
   ||p.cricoidSurface!==(i===0?'none':i===1?'remainder':suffix.replace(':','-'))
   ||!p.material||atlas.materials?.[p.material]?.pairedNormalSigns!==true
   ||p.chunk!==first.chunk||p.positions!==first.positions||p.normals!==first.normals||p.uvs!==first.uvs||p.vertexCount!==first.vertexCount
   ||JSON.stringify(p.geometryOrigin)!==JSON.stringify(first.geometryOrigin)
   ||!offset(p.positions)||!offset(p.normals)||!offset(p.uvs)||!offset(p.indices)||!offset(p.normalSigns)
   ||!Number.isInteger(p.vertexCount)||p.vertexCount<=0||!Number.isInteger(p.indexCount)||p.indexCount<=0||p.indexCount%3!==0
   ||p.positions+12*p.vertexCount>bank.bytes||p.normals+6*p.vertexCount>bank.bytes||p.uvs+8*p.vertexCount>bank.bytes
   ||p.indices+4*p.indexCount>bank.bytes||p.normalSigns+p.vertexCount>bank.bytes
   ||p.sectionAssembly!==`TA98:cricoid:${sex}:${i===0?'arch':'lamina'}`||p.sectionSolid!==(i===0)
   ||!sameMembers(p.sectionClosureParts,i===0?[p.id]:laminaIds))return reject();
 }
 if(atlas.parts.some(p=>!p.suppressed&&(p as CricoidPart).cricoidSourcePartId===sourceId&&!ids.includes(p.id)))return reject();
 // Reciprocal cartilage-only cuts stay in the common bank. The four open
 // articular records contribute to the lamina closure, never their own cap.
 const cutCount=(p:Part,neighbor:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===neighbor).reduce((n,g)=>n+g.count,0);
 for(const p of members as Part[]){
  const groups=p.tissueDrawGroups;if(!groups?.length)return reject();
  let end=0;
  for(const group of groups){
   if(!offset(group.start)||!Number.isInteger(group.count)||group.count<=0||group.count%3!==0||group.start!==end)return reject();
   end+=group.count;
   if(group.neighbor){
    const expected=p.id===ids[0]?ids[1]:p.id===ids[1]?ids[0]:undefined;
    const peer=members.find(row=>row?.id===group.neighbor);
    if(group.neighbor!==expected||!peer||cutCount(p,peer.id)!==cutCount(peer,p.id))return reject();
   }
  }
  if(end!==p.indexCount)return reject();
 }
 if(!cutCount(first,ids[1]))return reject();
 const whole=uniqueConcept(wholeId),canonical=terms.map(term=>uniqueConcept(`${sex==='female'?'TA98FEMALE':'TA98MALE'}:${term}`));
 const conceptMatches=(c:Concept|undefined,term:string,parts:readonly string[],parent?:string)=>!!c
  &&c.ta98Term===term&&c.ta98Kind==='exact'&&sameMembers(c.elements,parts)&&(parent===undefined||c.partOf===parent);
 if(!conceptMatches(whole,terms[0],ids)||!conceptMatches(canonical[0],terms[0],ids)
  ||!conceptMatches(canonical[1],terms[1],[ids[0]],wholeId)||!conceptMatches(canonical[2],terms[2],laminaIds,wholeId)
  ||!conceptMatches(canonical[3],terms[3],ids.slice(2,4),canonical[2]?.id)
  ||!conceptMatches(canonical[4],terms[4],ids.slice(4,6),canonical[2]?.id))return reject();
 for(let i=0;i<ids.length;i++){
  const c=uniqueConcept(ids[i]) as PartialConcept|undefined,k=i===0?1:i===1?2:i<4?3:4;
  if(!c||!sameMembers(c.elements,[ids[i]])||c.partOf!==canonical[k]!.id)return reject();
  if(i===1){if(c.ta98Term!==undefined||c.ta98Kind!=='parent'||c.ta98ParentTerm!==terms[2])return reject();}
  else if(c.ta98Term!==terms[k]||c.ta98Kind!=='exact')return reject();
 }
 const termParts=new Map(terms.map((term,i)=>[term,canonical[i]!.elements.map(id=>members.find(p=>p!.id===id)!)]));
 const result={marker,source,whole:whole!,arch:canonical[1]!,lamina:canonical[2]!,arytenoid:canonical[3]!,thyroid:canonical[4]!,members:members as Part[],termParts};
 registries.set(atlas,result);return result;
}

/** Only these five TA98 terms replace legacy source-name matching. */
export function cricoidPartitionTermParts(atlas:Atlas):ReadonlyMap<string,Part[]>|null {
 return cricoidPartitionRegistry(atlas)?.termParts??null;
}
