import type {Atlas,Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import reverifiedBanks from './data/ocular-reverified-banks.json';

type Mouth={targetId:string;contractSha256:string;boundaryVertices:6;opening:true};
type OcularPart=Part&{ocularMouthConnections?:Mouth[]};
const prefix='BDM:male:human_20_male_cardiovascular_system-';
/** Banks whose vessels were moved after certification: every certified opening was re-measured and still coincides,
 * so a part in such a bank counts as the certified bank it came from. */
const reverified=reverifiedBanks as Record<string,Record<string,string>>;
const role=(id:string)=>id.startsWith(prefix)?/^(left|right)_(ophthalmic_artery|short_posterior_ciliary_arteries|long_posterior_ciliary_arteries)_ID$/.exec(id.slice(prefix.length)):null;
const venousRole=(id:string)=>id.startsWith(prefix)?/^(left|right)_(inferior_ophthalmic_vein|superior_ophthalmic_vein|vorticose_veins|pterygoid_plexus|angular_vein|cavernous_sinus|central_retinal_vein|infraorbital_vein|anterior_superior_alveolar_vein|middle_superior_alveolar_vein)_ID(?::region:A12\.3\.06\.112)?$/.exec(id.slice(prefix.length)):null;

/** Attach the verified relation only to the exact installed geometry bank. */
export function withVerifiedAngularIovMouth(atlas:Atlas):Atlas{
 const bank='d65fa6eeaca3bad2d54a2b8988a559a7b8e9343c0d07e99902fdc0530f492482';
 const evidence='c77aa0991d231a1130c5b9751d551ed94a495989eb27f9ff1d7716dca752edf3';
 const ids=[prefix+'right_angular_vein_ID',prefix+'right_inferior_ophthalmic_vein_ID'];
 const data=atlas as Atlas&{chunks?:{sha256?:string}[]};
 const matches=(id:string,hash:string)=>{const owner=data.parts.find(part=>part.id===id) as (OcularPart&{chunk?:number})|undefined;return !!owner&&!owner.suppressed&&owner.system==='venous'&&owner.chunk!==undefined&&(data.chunks?.[owner.chunk]?.sha256===hash||reverified[data.chunks?.[owner.chunk]?.sha256??'']?.[id]===hash);};
 const pairs:{a:string;b:string;proof:string}[]=[];
 if(ids.every(id=>matches(id,bank)))pairs.push({a:ids[0],b:ids[1],proof:evidence});
 const vortexBank='fb09c2a26b19d9620bbf3df52a997c3ac110917524e58af9c7a9224a6bd2fdc2';
 const superiorBank='7b611ae5a0cc96c1fa7fb180d1624b36e9068a42dcd464fd4b8d3c9bcb29f715';
 const proof='74098448fe4fc8e80b9d05bc178a636cd2d9e93a18ec720be6cd841516fa6028';
 for(const side of ['left','right']){
  const vortex=prefix+side+'_vorticose_veins_ID';
  for(const [role,hash] of [['superior_ophthalmic_vein',superiorBank],['inferior_ophthalmic_vein',side==='left'?'dc4eca43e153e033a0f391ce237aac04b9c60c6849f78ed9b99effb7569c00da':bank]]){
   const target=prefix+side+'_'+role+'_ID';
   if((matches(vortex,vortexBank)||matches(vortex,'7d6f29d24df2d882d7acd2d54d228550152bff2a503123db66744537fbabab2b'))&&matches(target,hash))pairs.push({a:vortex,b:target,proof});
  }
  const verified='f746803e0d0d2f98d00270e3eefd4c4abf3c3f68f2fa6ef1e6d818dd2a994e12';
  const superior=prefix+side+'_superior_ophthalmic_vein_ID',inferior=prefix+side+'_inferior_ophthalmic_vein_ID';
  const inferiorBank=side==='left'?'dc4eca43e153e033a0f391ce237aac04b9c60c6849f78ed9b99effb7569c00da':bank;
  const junctions=[
   [superior,superiorBank,prefix+side+'_cavernous_sinus_ID','e19e1c6bda8c761c167904ac285924f9f89e333cee354e08d5c77f2646bac753'],
   [inferior,inferiorBank,prefix+side+'_cavernous_sinus_ID','e19e1c6bda8c761c167904ac285924f9f89e333cee354e08d5c77f2646bac753'],
   [superior,superiorBank,prefix+side+'_central_retinal_vein_ID:region:A12.3.06.112','b219ae1c22b81f2c1ed471372bc90db9c5d541a1251d71c1f656f434e11a04b4'],
   [inferior,inferiorBank,prefix+side+'_pterygoid_plexus_ID','dc4eca43e153e033a0f391ce237aac04b9c60c6849f78ed9b99effb7569c00da'],
   [inferior,inferiorBank,prefix+side+'_angular_vein_ID',side==='left'?'111004a8dfe4d0bb9d8215563b65ed641b18cdcca7d7f7580b7d67f64d0ce519':bank],
  ];
  for(const [a,ah,b,bh] of junctions)if(matches(a,ah)&&matches(b,bh))pairs.push({a,b,proof:verified});
  const orbital=prefix+side+'_infraorbital_vein_ID';
  const facialBank='b219ae1c22b81f2c1ed471372bc90db9c5d541a1251d71c1f656f434e11a04b4';
  const orbitalProof='7b33098b06c665e52267c97a1b3c21f1034ef9d0514d39616b6618cb79bee0de';
  for(const [name,hash] of [['angular_vein',side==='left'?'111004a8dfe4d0bb9d8215563b65ed641b18cdcca7d7f7580b7d67f64d0ce519':bank],['anterior_superior_alveolar_vein',facialBank],['middle_superior_alveolar_vein',facialBank]]){
   const peer=prefix+side+'_'+name+'_ID';
   if(matches(orbital,facialBank)&&matches(peer,hash))pairs.push({a:orbital,b:peer,proof:orbitalProof});
  }
 }
 if(!pairs.length)return atlas;
 return {...atlas,parts:atlas.parts.map(part=>{const links=pairs.filter(pair=>pair.a===part.id||pair.b===part.id);if(!links.length)return part;let prior=(part as OcularPart).ocularMouthConnections??[];for(const link of links){const targetId=link.a===part.id?link.b:link.a;prior=[...prior.filter(item=>item.targetId!==targetId),{targetId,contractSha256:link.proof,boundaryVertices:6,opening:true}];}return {...part,ocularMouthConnections:prior};})};
}

/** Only reciprocal, explicitly named and source-verified mouths establish continuity. */
export function ocularMouthRelations(source:Part,parts:Part[]):ResolvedRelation[]{
 const sourceRole=role(source.id);
 const sourceVenousRole=venousRole(source.id);
 if(source.suppressed||!(sourceRole&&source.system==='arterial'||sourceVenousRole&&source.system==='venous'))return [];
 const result:ResolvedRelation[]=[];
 for(const mouth of (source as OcularPart).ocularMouthConnections??[]){
  if(mouth.opening!==true||mouth.boundaryVertices!==6||!/^[a-f0-9]{64}$/.test(mouth.contractSha256))continue;
  const target=parts.find(part=>part.id===mouth.targetId) as OcularPart|undefined;
  if(!target||target.suppressed||target.system!==source.system)continue;
  const targetRole=role(target.id),targetVenousRole=venousRole(target.id);
  const arterial=!!(sourceRole&&targetRole&&targetRole[1]===sourceRole[1]&&
   (sourceRole[2]==='ophthalmic_artery')!==(targetRole[2]==='ophthalmic_artery'));
  const roles=[sourceVenousRole?.[2],targetVenousRole?.[2]];
  const venous=!!(sourceVenousRole&&targetVenousRole&&targetVenousRole[1]===sourceVenousRole[1]&&sourceVenousRole[2]!==targetVenousRole[2]&&(roles.includes('infraorbital_vein')?(roles.includes('angular_vein')||roles.includes('anterior_superior_alveolar_vein')||roles.includes('middle_superior_alveolar_vein')):roles.includes('vorticose_veins')?(roles.includes('superior_ophthalmic_vein')||roles.includes('inferior_ophthalmic_vein')):roles.includes('cavernous_sinus')?(roles.includes('superior_ophthalmic_vein')||roles.includes('inferior_ophthalmic_vein')):roles.includes('central_retinal_vein')?roles.includes('superior_ophthalmic_vein'):(roles.includes('inferior_ophthalmic_vein')&&(roles.includes('angular_vein')||roles.includes('pterygoid_plexus')))));
  if(!arterial&&!venous)continue;
  if(!target.ocularMouthConnections?.some(other=>other.targetId===source.id&&other.opening===true&&other.boundaryVertices===6&&other.contractSha256===mouth.contractSha256))continue;
  if(!result.some(relation=>relation.target.id===target.id))result.push({kind:'continuous',target,note:venous?(roles.includes('infraorbital_vein')?'These veins have matching vessel openings.':roles.includes('vorticose_veins')?'Vorticose and ophthalmic veins have matching vessel openings.':roles.includes('angular_vein')?'Angular and inferior ophthalmic veins share a vessel opening.':roles.includes('cavernous_sinus')?'Ophthalmic vein and cavernous sinus have matching openings.':roles.includes('central_retinal_vein')?'Retinal and superior ophthalmic veins have matching openings.':'Open inferior ophthalmic–pterygoid venous connection.'):'Open ciliary–ophthalmic arterial connection.'});
 }
 return result;
}
