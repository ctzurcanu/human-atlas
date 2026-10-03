/** Source pigment classes for the local dermatome annotation (not tissue ownership). */
export const dermatomeColors:Record<string,number[]>={
 C2:[255,255,255],C3:[246,240,250],C4:[231,223,243],C5:[211,202,230],C6:[198,184,219],C7:[172,152,204],C8:[146,121,189],
 T1:[251,237,237],T2:[250,230,230],T3:[248,220,220],T4:[245,206,206],T5:[245,201,201],T6:[243,192,192],T7:[241,184,184],T8:[240,174,174],T9:[238,163,163],T10:[237,154,154],T11:[235,147,147],T12:[233,136,136],
 L1:[227,231,249],L2:[198,217,248],L3:[171,199,247],L4:[139,176,247],L5:[118,117,247],
 S1:[173,216,176],S2:[141,201,146],S3:[115,187,122],S4:[73,165,80],S5:[49,121,57],
 V1:[255,205,162],V2:[252,171,140],V3:[248,155,122],
};
export function dermatomeSelection(id=''):{level?:string;side?:string;region?:string}{
 const match=/(?:^|:)(?:DERMATOME|ROOT|TRIGEMINAL):(C\d|T\d{1,2}|L\d|S\d|Co1|V[123])(?::(left|right))?$/.exec(id);
 return match?{level:match[1],side:match[2]}:/(?:^|:)TRIGEMINAL$/.test(id)?{region:'V'}:{region:/(?:^|:)REGION:(C|T|L|S|Co)$/.exec(id)?.[1]};
}
export function territoryPalette(sex:'male'|'female',head:boolean):[string,number[]][]{
 // The native head atlas extends onto the neck/upper back; its UV name is not a sensory boundary.
 const keys=head?['C2','C3','C4','C5','C6','C7','C8','T1','T2','V1','V2','V3']:Object.keys(dermatomeColors).filter(k=>!/^V|^C[234]$/.test(k));
 const femaleThoracic=[[249,236,236],[247,230,230],[247,219,218],[244,204,203],[241,195,198],[239,187,186],[239,179,183],[232,168,170],[238,163,162],[235,154,152],[233,146,145],[232,136,137]];
 return keys.map(key=>[key,sex==='female'&&/^T/.test(key)?femaleThoracic[Number(key.slice(1))-1]:dermatomeColors[key]]);
}

export const dermatomeLevels=['C2','C3','C4','C5','C6','C7','C8',...Array.from({length:12},(_,i)=>`T${i+1}`),...Array.from({length:5},(_,i)=>`L${i+1}`),...Array.from({length:5},(_,i)=>`S${i+1}`),'V1','V2','V3'];
export const dermatomeTerritories=dermatomeLevels.flatMap(level=>['left','right'].map(side=>`${level}:${side}`));
/** Shared by checkbox state and fragment visibility; unmapped C1/Co1 have no patch. */
export function dermatomeTerritoriesFor(id:string):string[]{
 if(id==='all')return dermatomeTerritories;
 const {level,region,side}=dermatomeSelection(id);
 return dermatomeTerritories.filter(key=>{const [root,laterality]=key.split(':');return !!(level||region)&&(!level||root===level)&&(!region||root.startsWith(region))&&(!side||side===laterality);});
}

export function updateDermatomeSelection(current:string[],incoming:string[],additive=false):string[]{
 const valid=new Set(dermatomeTerritories),keys=[...new Set(incoming.filter(key=>valid.has(key)))];
 if(!additive)return keys;
 const next=new Set(current.filter(key=>valid.has(key))),remove=keys.every(key=>next.has(key));
 for(const key of keys){if(remove)next.delete(key);else next.add(key);}return [...next];
}
export const dermatomeLabel=(key:string)=>{const [level,side]=key.split(':');return `${side==='left'?'Left':'Right'} ${level} ${level.startsWith('V')?'facial sensory territory':'dermatome'}`;};
export function dermatomeChoice(keys:string[]):import('./hierarchy-choice').HierarchyChoice|null{
 if(!keys.length)return null;
 const children=keys.map(key=>{const [level,side]=key.split(':');return {id:`guest:dermatomes-myotomes:${level.startsWith('V')?'TRIGEMINAL':'DERMATOME'}:${level}:${side}`,name:dermatomeLabel(key),elements:[] as string[]};});
 return children.length===1?children[0]:{id:'guest:dermatomes-myotomes:DERMATOME-SELECTION',name:`${keys.length} selected sensory territories`,elements:[],children};
}
