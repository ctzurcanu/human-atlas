import {explanation,structureName,type Part} from './anatomy';
import type {AnatomyPathNode} from './anatomy-path';

export interface StructureDescription {text:string;url?:string}

const REGION=/^(?:head|neck|brain|face|skull|trunk|torso|thorax|chest|abdomen|pelvis|back|(?:left|right) (?:upper|lower) limb|arm|forearm|hand|thigh|leg|foot)$/i;
const KIND:Partial<Record<Part['system'],string>>={
 skeletal:'a bone or bony region',muscular:'a muscle',arterial:'an artery',venous:'a vein',
 nervous:'a nervous-system structure',digestive:'a digestive structure',urinary:'a urinary structure',
 respiratory:'a respiratory structure',reproductive:'a reproductive structure',
 lymphatic:'a lymphatic structure',endocrine:'an endocrine structure',
 connective:'a connective tissue structure',fascia:'a fascial structure',
 integumentary:'a surface region',sensory:'a sensory structure',cardiac:'a cardiac structure',
 attachments:'a muscle attachment',regions:'an anatomical region',
};
const clean=(name:string)=>structureName(name).replace(/\s*\((?:left|right)\)$/i,'').trim();
const lowerFirst=(value:string)=>value[0]?.toLowerCase()+value.slice(1);

export function localDescription(name:string,part:Part|undefined,path:AnatomyPathNode[]):StructureDescription{
 const specific=explanation(name);
 if(specific)return {text:specific};
 const subject=clean(name),kind=part?KIND[part.system]??'an anatomical structure':'an anatomical structure';
 const side=/\((left|right)\)$/i.exec(structureName(name))?.[1]?.toLowerCase();
 const groups=path.slice(2,-1).map(node=>node.label);
 const region=groups.find(group=>/^(?:left|right) (?:upper|lower) limb$/i.test(group))??groups.find(group=>REGION.test(group));
 const parent=[...groups].reverse().find(group=>group.toLowerCase()!==subject.toLowerCase()&&!REGION.test(group));
 const subjectLabel=side?`${subject} on the ${side}`:subject;
 let relation='';
 const partOf=/^(?:body|head|neck|shaft|base|root|apex|process|lobe|segment|branch|surface|border|crest) of (.+)$/i.exec(subject);
 if(partOf)relation=` It is a named part of ${lowerFirst(partOf[1])}.`;
 else if(parent)relation=` In this atlas it is grouped with ${lowerFirst(parent)}.`;
 return {text:`The ${lowerFirst(subjectLabel)} is ${kind}${region?` in the ${lowerFirst(region)}`:''}.${relation}`};
}

const summaryCache=new Map<string,Promise<StructureDescription|null>>();
export function wikipediaDescription(name:string):Promise<StructureDescription|null>{
 const title=clean(name),key=title.toLowerCase();
 const cached=summaryCache.get(key);if(cached)return cached;
 const params=new URLSearchParams({action:'query',prop:'extracts|pageprops',exintro:'1',explaintext:'1',redirects:'1',format:'json',origin:'*',titles:title});
 const request=fetch(`https://en.wikipedia.org/w/api.php?${params}`).then(async response=>{
  if(!response.ok)return null;
  const data=await response.json() as {query?:{redirects?:{tofragment?:string}[];pages?:Record<string,{title?:string;extract?:string;missing?:boolean;pageprops?:{disambiguation?:string}}>}};
  const page=Object.values(data.query?.pages??{})[0];
  if(!page?.title||page.missing||page.pageprops?.disambiguation!==undefined||!page.extract)return null;
  // A section redirect often points to an article about a broader structure.
  // Show our part-specific description instead of presenting that lead as a match.
  if(data.query?.redirects?.some(redirect=>redirect.tofragment))return null;
  const tokens=(value:string)=>new Set(value.toLowerCase().replaceAll('oes','es').match(/[a-z]{3,}/g)?.filter(word=>!['the','and','of','part','left','right'].includes(word))??[]);
  const requested=tokens(title),matched=tokens(page.title);
  if(requested.size&&[...requested].filter(word=>matched.has(word)).length/Math.min(requested.size,matched.size||1)<.6)return null;
  const sentences=page.extract.replace(/\s+/g,' ').trim().split(/(?<=[.!?])\s+(?=[A-Z])/u);
  const text=sentences.slice(0,2).join(' ').slice(0,650).trim();
  if(text.length<35)return null;
  return {text,url:`https://en.wikipedia.org/wiki/${encodeURIComponent(page.title.replaceAll(' ','_'))}`};
 }).catch(()=>{summaryCache.delete(key);return null;});
 summaryCache.set(key,request);
 return request;
}
