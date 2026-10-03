import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const nameOf=(part:Part)=>part.name.toLowerCase().trim();
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
export function cardiacValveDescription(part:Part):{text:string;url:string}|undefined{
 if(part.system!=='cardiac')return;
 const valve=valves.find(v=>v.pattern.test(nameOf(part)));if(!valve)return;
 const subject=/\b(?:leaflet|cusp)\b/.test(nameOf(part))?`This is a leaflet of the ${valve.name.toLowerCase()} valve.`:`The ${valve.name.toLowerCase()} valve controls flow from the ${valve.up[0]} to the ${valve.down[0]}.`;
 const detail=/\b(?:leaflet|cusp)\b/.test(nameOf(part))?` The valve controls flow from the ${valve.up[0]} to the ${valve.down[0]}.`:'';
 return {text:`${subject}${detail} Its leaflets open and close to prevent backward blood flow.`,url:'https://www.nhlbi.nih.gov/health/heart/blood-flow'};
}

/** Named valvular relationships, never an assertion of certified mesh contact.
 * NIH NHLBI Heart valves; OpenStax A&P 2e 19.1. Donor HRA organs do not
 * acquire a native frame just because their labels describe the same heart. */
export function cardiacRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),frames=new Map<string,Part[]>();
 for(const part of parts){const frame=nativePartFrame(part);if(!frame||part.suppressed||!['cardiac','arterial'].includes(part.system))continue;const list=frames.get(frame)??[];list.push(part);frames.set(frame,list);}
 const add=(source:Part,target:Part,note:string)=>{
  if(source.id===target.id)return;
  for(const [from,to] of [[source,target],[target,source]]){const list=graph.get(from.id)??[];if(!list.some(r=>r.target.id===to.id))list.push({kind:'adjacent',target:to,note});graph.set(from.id,list);}
 };
 for(const peers of frames.values()){
  const named=(names:string[])=>peers.filter(p=>names.includes(nameOf(p))&&(p.system==='cardiac'||names.includes('ascending aorta')||names.includes('pulmonary trunk')));
  for(const valve of valves)for(const source of peers.filter(p=>p.system==='cardiac'&&valve.pattern.test(nameOf(p)))){
   const note=`${valve.name} valve region separates ${valve.up[0]} from ${valve.down[0]}. Named chamber/outflow adjacency; source leaflet extent, annular contact and lumen continuity remain unverified.`;
   for(const target of [...named(valve.up),...named(valve.down)])add(source,target,note);
   // A single ventricular envelope may be shown as context, but never renamed
   // or counted as an independently represented left/right ventricular chamber.
   if([...valve.up,...valve.down].some(n=>n.includes('ventricle')))for(const target of named(['ventricles']))add(source,target,`${note} This source is the combined ventricular envelope; separate left/right ventricular representation is missing.`);
  }
  for(const source of peers.filter(p=>p.system==='cardiac'&&/papillary.*(?:left|right) ventricle/.test(nameOf(p)))){
   const side=nameOf(source).includes('left ventricle')?'left':'right';
   for(const target of named(side==='left'?chambers.leftVentricle:chambers.rightVentricle))add(source,target,`Named ${side} ventricular papillary muscle projects from its ventricular wall. Chordal attachment and shared source tissue interfaces remain unverified; semilunar valves have no papillary/chordal apparatus.`);
  }
 }
 cache.set(parts,graph);return graph;
}
