import metadata from './data/ta98-metadata.json';

export type Terminology={ta98:string|null;tha:string|null;fma:string|null;latin:string|null;ontology?:string|null};
const empty:Terminology={ta98:null,tha:null,fma:null,latin:null};
const concepts=metadata.byConcept as Record<string,Terminology>;
const groups=metadata.byGroup as Record<string,Terminology>;
const groupsLower=new Map(Object.entries(groups).map(([name,term])=>[name.toLowerCase(),term]));
const paths=metadata.byConceptPath as Record<string,string[]>;
const codes=metadata.byCode as Record<string,Terminology&{name:string}>;
const matches=metadata.byConceptMatch as Record<string,{term:string;kind:'exact'|'parent';method:string}>;

export function taChapterForConcept(id:string):string|undefined{
 return taEntityForConcept(id)?.slice(0,3);
}

export function taEntityForConcept(id:string):string|undefined{
 return matches[id]?.term;
}

export function terminologyForConcept(id:string):Terminology{
 return concepts[id]??(/^FMA:?\d+$/.test(id)?{...empty,fma:`FMA:${id.replace(/^FMA:?/,'')}`}:empty);
}

export function terminologyForGroup(name:string):Terminology{
 if(groups[name])return groups[name];
 const side=/^(?:Left|Right) (.+)$/.exec(name)?.[1];
 return groupsLower.get((side??name).toLowerCase())??empty;
}

export function taHierarchyForConcept(id:string):{name:string;terminology:Terminology}[]{
 return (paths[id]??[]).map(code=>codes[code]).filter((term):term is Terminology&{name:string}=>!!term).map(term=>({name:term.name.replace(/^./,letter=>letter.toUpperCase()),terminology:term}));
}

export function terminologyTitle(name:string,term:Terminology,id?:string):string{
 const lines=[name];
 if(term.latin)lines.push(`La:${term.latin}`);
 if(term.ta98)lines.push(`TA98:${term.ta98}`);
 if(term.tha)lines.push(term.tha);
 if(term.fma)lines.push(term.fma);
 if(term.ontology)lines.push(term.ontology);
 if(id){const local=atlasIdentifier(id);if(local)lines.push(local);}
 return lines.join('\n');
}

/** Stable local ID for a modeled leaf or a browser grouping without a TA98 code. */
export function atlasIdentifier(id:string):string|null{
 if(id==='selection-set')return null;
 if(!id.startsWith('hierarchy:')&&!id.startsWith('guest:'))return `Atlas:${id}`;
 let hash=2166136261;
 for(const character of id){hash^=character.charCodeAt(0);hash=Math.imul(hash,16777619);}
 return `HA-G:${(hash>>>0).toString(16).padStart(8,'0')}`;
}
