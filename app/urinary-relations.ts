import type {Part} from './anatomy';
import type {ResolvedRelation} from './anatomical-relations';
import {nativePartFrame} from './ta98-modeled-parent';

/** Only the stored ZA urinary partitions retain their original native frame. */
export function urinaryPartFrame(p:Part):string|undefined{
 const native=nativePartFrame(p);if(native)return native;
 if(p.system!=='urinary')return;
 const m=p.id.match(/^TA98PART:ZA:(Ureter\.[lr]|Urinary bladder|Urethra):(A\d{2}\.\d\.\d{2}\.\d{3})$/);
 if(m&&p.sourceId===m[1])return 'ZA:';
}
const name=(p:Part)=>p.name.toLowerCase().replace(/\s*\((?:left|right)\)$/,'');
const side=(p:Part)=>p.name.match(/\((left|right)\)$/i)?.[1]?.toLowerCase();
function stage(p:Part):number|undefined{
 if(p.suppressed||p.system!=='urinary'||!urinaryPartFrame(p))return;
 const n=name(p);
 if(n==='kidney')return 0;
 if(n==='renal pelvis')return 1;
 if(n==='ureter'||n==='abdominal part of ureter')return 2;
 if(n==='pelvic part of ureter')return 3;
 if(n==='intramural part of ureter')return 4;
 if(['bladder','urinary bladder','body of bladder','fundus of bladder','neck of bladder'].includes(n))return 5;
 if(['urethra','female urethra','male urethra','prostatic urethra'].includes(n))return 6;
 if(n==='intermediate part of urethra')return 7;
 if(n==='spongy urethra')return 8;
}
export const hasNativeUrinaryRoute=(p:Part)=>stage(p)!==undefined;
export function urinarySourceCoverageLimitation(p:Part):string|undefined{
 if(p.conceptId==='LOCAL:female:urethra')return 'Named native female urethral exterior. Intramural extent, wall layers, glands, sphincter interfaces and a patent lumen remain unverified. The external urethral opening is not a separate solid.';
}
const cached=new WeakMap<Part[],Map<string,ResolvedRelation[]>>();
export function nativeUrinaryRelations(parts:Part[]):Map<string,ResolvedRelation[]>{
 const prior=cached.get(parts);if(prior)return prior;
 const graph=new Map<string,ResolvedRelation[]>(),available=parts.filter(hasNativeUrinaryRoute);
 const labels=['Kidney','Renal pelvis','Abdominal part of ureter','Pelvic part of ureter','Intramural part of ureter','Urinary bladder','Prostatic urethra','Intermediate part of urethra','Spongy urethra'];
 for(const source of available){
  const level=stage(source)!,frame=urinaryPartFrame(source),which=side(source),links:ResolvedRelation[]=[];
  for(const [kind,direction] of [['before',-1],['after',1]] as const){
   const missing:string[]=[];
   for(let next=level+direction;next>=0&&next<labels.length;next+=direction){
    // A whole native ureter or urethra is already the complete named stage;
    // never fabricate subdivisions downstream of that whole exterior.
    if(name(source)==='urethra'||name(source)==='female urethra'||name(source)==='male urethra'){if(direction===1)break;}
    if(name(source)==='ureter'&&next>=3&&next<=4)continue;
    const targets=available.filter(p=>stage(p)===next&&urinaryPartFrame(p)===frame&&(!which||!side(p)||side(p)===which));
    if(targets.length){
     for(const target of targets)links.push({kind,target,...(missing.length?{via:direction===1?missing:[...missing].reverse()}:{}),note:'Urine path within the same native model. Bladder and urethral source regions reference the whole organ cavity; their exterior meshes do not establish local lumen continuity, ureteral ostia, wall layers or sphincter anatomy. Missing native intermediates are shown as text, not substituted with independently positioned reference organs.'});
     break;
    }
    // A whole ureter must also be found when walking backwards from bladder.
    if(next===4||next===3){const whole=available.filter(p=>name(p)==='ureter'&&urinaryPartFrame(p)===frame&&(!which||!side(p)||side(p)===which));if(whole.length){for(const target of whole)links.push({kind,target,note:'Urine path through the whole native ureter to the bladder. Intramural course, ureteral ostium and lumen continuity are unverified.'});break;}}
    missing.push(labels[next]);
   }
  }
  graph.set(source.id,links);
 }
 cached.set(parts,graph);return graph;
}
