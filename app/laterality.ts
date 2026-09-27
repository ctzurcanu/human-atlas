export type Laterality='left'|'right'|null;

// Display-only formatting. Atlas names and IDs remain intact for lookup,
// search, accessibility labels, links, and the selected structure heading.
export function displayLaterality(name:string):{label:string;side:Laterality}{
 const parenthetical=/\s*\((left|right|[LR])\)(?=\s*(?:·|$))/i.exec(name);
 if(parenthetical)return {label:name.replace(parenthetical[0],'').replace(/\s+·/,' ·').trim(),side:sideOf(parenthetical[1])};
 const prefix=/^(left|right)\s+/i.exec(name);
 if(prefix)return {label:name.slice(prefix[0].length),side:sideOf(prefix[1])};
 const embedded=/\bof (left|right)\s+/i.exec(name);
 if(embedded)return {label:name.replace(embedded[0],'of '),side:sideOf(embedded[1])};
 // Z-Anatomy attachment names encode origin/insertion site, ordinal, and side
 // together: .or/.ol, .e10r/.e10l, .o3r/.o3l. Keep the site number visible.
 const attachment=/\.([oe]\d*)([lr])(?=\s*(?:·|$))/i.exec(name);
 if(attachment)return {label:name.replace(attachment[0],`.${attachment[1]}`),side:sideOf(attachment[2])};
 const abbreviated=/[\s._]([LR])(?=\s*(?:·|$))/i.exec(name);
 if(abbreviated)return {label:name.replace(abbreviated[0],'').replace(/\s+·/,' ·').trim(),side:sideOf(abbreviated[1])};
 return {label:name,side:null};
}

const sideOf=(value:string):Exclude<Laterality,null>=>/^(?:left|l)$/i.test(value)?'left':'right';
export const lateralityClass=(side:Laterality)=>side?`side-${side}`:'';
