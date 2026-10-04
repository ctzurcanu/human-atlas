import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

const cache=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
const handDigital=(p:Part)=>/^LOCAL:female:[lr]_dorsal_digital_veins$/.test(p.conceptId)||/^ZA:Dorsal digital veins of hand\.[lr]$/.test(p.conceptId);
const palmarDigital=(p:Part)=>/^LOCAL:female:[lr]_proper_palmar_digital_veins$/.test(p.conceptId)||/^ZA:Palmar digital veins\.[lr]$/.test(p.conceptId);
const palmarArch=(p:Part)=>/^LOCAL:female:[lr]_superficial_palmar_venous_arch$/.test(p.conceptId)||/^ZA:Superficial venous palmar arch\.[lr]$/.test(p.conceptId);
const side=(p:Part)=>p.name.match(/\((left|right)\)$/i)?.[1]?.toLowerCase();
const name=(p:Part)=>p.name.toLowerCase().replace(/\s*\((left|right)\)$/,'');
export function handVesselCoverageLimitation(p:Part):string|undefined{
 if(p.system==='arterial'&&/^LOCAL:female:[lr]_palmar_digital_arteries$/.test(p.conceptId))return 'This native source combines common and proper palmar digital arterial branches. It is grouped under common palmar digital arteries; separate common/proper meshes and digit-specific supply territories remain unverified. Known distal surface self-crossings remain; this geometry is provisional. Source walls, lumens, junctions and tissue packing remain unverified.';
 if(p.system==='arterial'&&/^LOCAL:female:[lr]_(?:superficial_palmar_branch_radial_artery|deep_palmar_branch_ulnar_artery)$/.test(p.conceptId))return 'Named native palmar arterial branch. Source junction coordinates are shared with the adjoining artery and arch, but this does not establish vessel walls, lumen continuity, complete supply territories or tissue packing.';
 if(p.system==='venous'&&handDigital(p))return 'These are native dorsal digital branches of the hand. They are grouped under the dorsal venous network of hand, not identified as foot veins or as its complete network. Individual tributaries, vessel walls, lumens and physical junctions remain unverified.';
 if(p.system==='venous'&&palmarDigital(p))return 'Native palmar digital venous branches of the hand. '+(p.conceptId.startsWith('LOCAL:female:')?'This source is a partial proper digital network, placed under palmar digital veins without claiming their complete representation. ':'')+'Source tributaries, walls, lumens, junctions and tissue packing remain unverified.';
}

/** Group-level drainage; not a fixed serial route for every individual vein.
 * Original anatomical studies: PMID 4020057 and PMID 28417223 document
 * digital venous connections and variable dorsal metacarpal arrangements.
 */
export function handVascularRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cache.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>();
 const note='Group-level dorsal hand drainage through digital/metacarpal tributaries and the dorsal hand network toward the cephalic and basilic outlets. Individual tributaries and anatomical variants are not fixed by this link; source walls, lumens, physical junctions and packing remain unverified.';
 const pair=(a:Part,b:Part,via:string[]=[],relationNote=note)=>{
  for(const [from,to,kind] of [[a,b,'after'],[b,a,'before']] as const){
   const links=graph.get(from.id)??[];
   if(!links.some(r=>r.kind===kind&&r.target.id===to.id))links.push({kind,target:to,note:relationNote,...(via.length?{via}: {})});
   graph.set(from.id,links);
  }
 };
 for(const digital of parts.filter(p=>!p.suppressed&&p.system==='venous'&&handDigital(p))){
  const frame=nativePartFrame(digital),which=side(digital);if(!frame||!which)continue;
  const find=(label:string)=>parts.filter(p=>!p.suppressed&&p.system==='venous'&&nativePartFrame(p)===frame&&side(p)===which&&name(p)===label);
  const metacarpal=find('dorsal metacarpal veins'),network=find('dorsal venous network of hand'),outlets=[...find('cephalic vein'),...find('basilic vein')];
  if(metacarpal.length){for(const p of metacarpal)pair(digital,p);}
  else if(network.length){for(const p of network)pair(digital,p,['Dorsal metacarpal veins']);}
  // Never invent an absent collector mesh, inherit the foot's drainage or
  // connect a digital branch to a donor/outlet in another source frame.
  for(const m of metacarpal){
   if(network.length){for(const n of network)pair(m,n);}
   else for(const outlet of outlets)pair(m,outlet,['Dorsal venous network of hand']);
  }
  for(const n of network)for(const outlet of outlets)pair(n,outlet);
 }
 for(const digital of parts.filter(p=>!p.suppressed&&p.system==='venous'&&palmarDigital(p))){
  const frame=nativePartFrame(digital),which=side(digital);if(!frame||!which)continue;
  const arches=parts.filter(p=>!p.suppressed&&p.system==='venous'&&palmarArch(p)&&nativePartFrame(p)===frame&&side(p)===which);
  for(const arch of arches)pair(digital,arch,[],'Palmar digital venous tributaries connect with the superficial venous palmar arch. This is a group-level association; variable dorsal communications, individual tributaries, vessel walls and lumen continuity are not established by the source meshes.');
 }
 const arterialRoutes:[string[],string[]][]=[
  [['radial artery'],['superficial palmar branch radial artery','superficial palmar branch of radial artery']],
  [['superficial palmar branch radial artery','superficial palmar branch of radial artery'],['superficial palmar arterial arch','superficial palmar arch']],
  [['ulnar artery'],['deep palmar branch ulnar artery','deep palmar branch of ulnar artery']],
  [['deep palmar branch ulnar artery','deep palmar branch of ulnar artery'],['deep palmar arterial arch','deep palmar arch']],
  [['ulnar artery'],['superficial palmar arterial arch','superficial palmar arch']],
  [['radial artery'],['deep palmar arterial arch','deep palmar arch']],
  [['superficial palmar arterial arch','superficial palmar arch'],['common palmar digital arteries','palmar digital arteries']],
  [['common palmar digital arteries'],['proper palmar digital arteries']],
  [['deep palmar arterial arch','deep palmar arch'],['palmar metacarpal arteries']],
 ];
 const arteries=parts.filter(p=>!p.suppressed&&p.system==='arterial'&&nativePartFrame(p)&&side(p));
 for(const [from,to] of arterialRoutes)for(const a of arteries.filter(p=>from.includes(name(p))))for(const b of arteries.filter(p=>to.includes(name(p))&&nativePartFrame(p)===nativePartFrame(a)&&side(p)===side(a)))pair(a,b,[],'Named palmar arterial branching in this native hand. Arch contributions and digital supply can vary; this association does not verify individual source walls, lumens, junctions or complete supply territories.'+(name(b)==='palmar digital arteries'?' This digital network remains unpartitioned.':''));
 cache.set(parts,graph);return graph;
}
