import type {Atlas,Part} from './anatomy';

type BodyRegion='head'|'upper'|'lower'|'core';
type SortKey={depth:number;distal:number;inferior:number};
const center=(part:Part)=>part.bounds[0].map((low,axis)=>(low+part.bounds[1][axis])/2);

function regionFor(part:Part):BodyRegion{
 const groups=(part.groups??[]).join(' ').toLowerCase();
 if(/\b(upper limb|upper extremity|forearm|arm|hand|wrist|shoulder)\b/.test(groups))return 'upper';
 if(/\b(lower limb|lower extremity|thigh|leg|foot|ankle)\b/.test(groups))return 'lower';
 if(/\b(head|face|skull|brain|cranium|cerebrum|cerebellum)\b/.test(groups))return 'head';
 const [x,y]=center(part);
 if(y>1.45)return 'head';
 if(Math.abs(x)>.18&&y>.65)return 'upper';
 if(y<.85)return 'lower';
 return 'core';
}

function distanceToBounds(point:number[],part:Part){
 const squared=point.reduce((sum,value,axis)=>{
  const gap=Math.max(part.bounds[0][axis]-value,0,value-part.bounds[1][axis]);
  return sum+gap*gap;
 },0);
 return Math.sqrt(squared);
}

/** Source depth wins when supplied; otherwise distance from nearby bone is a local depth estimate. */
export function createDepthOrder(atlas:Atlas){
 const bones=atlas.parts.filter(part=>part.system==='skeletal').map(part=>({part,region:regionFor(part),side:Math.sign(center(part)[0])}));
 const keys=new WeakMap<Part,SortKey>();
 const key=(part:Part):SortKey=>{
  const existing=keys.get(part);if(existing)return existing;
  const point=center(part),region=regionFor(part),side=Math.sign(point[0]);
  const nearby=bones.filter(bone=>bone.region===region&&(side===0||bone.side===0||bone.side===side)&&Math.abs(center(bone.part)[1]-point[1])<.22);
  const candidates=nearby.length?nearby:bones.filter(bone=>bone.region===region);
  const boneDistance=candidates.length?Math.min(...candidates.map(bone=>distanceToBounds(point,bone.part))):0;
  const thickness=region==='core'?.16:region==='head'?.09:region==='upper'?.075:.08;
  const name=part.name.toLowerCase();
  const depth=Number.isFinite(part.depth)?-part.depth!:Math.round(Math.min(2,boneDistance/thickness)*5)/5+(/\bsuperficial\b/.test(name)?.2:0)-(/\bdeep\b/.test(name)?.2:0);
  const anchorX=side*(region==='upper'?.17:.09),anchorY=region==='upper'?1.38:.87;
  const distal=region==='upper'||region==='lower'?Math.round(Math.hypot(point[0]-anchorX,point[1]-anchorY,point[2])*20)/20:0;
  const result={depth,distal,inferior:point[1]};keys.set(part,result);return result;
 };
 return (a:{name:string;parts:Part[]},b:{name:string;parts:Part[]})=>{
  const average=(parts:Part[])=>parts.reduce((sum,part)=>{const value=key(part);return {depth:sum.depth+value.depth,distal:sum.distal+value.distal,inferior:sum.inferior+value.inferior};},{depth:0,distal:0,inferior:0});
  const left=average(a.parts),right=average(b.parts),aCount=a.parts.length,bCount=b.parts.length;
  return right.depth/bCount-left.depth/aCount||right.distal/bCount-left.distal/aCount||left.inferior/aCount-right.inferior/bCount||a.name.localeCompare(b.name,undefined,{numeric:true,sensitivity:'base'});
 };
}
