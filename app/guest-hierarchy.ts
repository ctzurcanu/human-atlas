import {SYSTEMS,structureName,type Atlas,type Part,type SystemId} from './anatomy';
import {createDepthOrder} from './depth-sort';

export type GuestNode={name:string;matches?:string[];systems?:SystemId[];children?:GuestNode[];unmapped?:boolean};
export type GuestHierarchy={schema:'human-atlas-hierarchy/v1';id:string;name:string;source?:string;nodes:GuestNode[]};
export type ResolvedGuestNode={id:string;name:string;parts:Part[];directParts:Part[];children:ResolvedGuestNode[]};

const object=(value:unknown):value is Record<string,unknown>=>!!value&&typeof value==='object'&&!Array.isArray(value);
const strings=(value:unknown,max=100):value is string[]=>Array.isArray(value)&&value.length<=max&&value.every(item=>typeof item==='string'&&item.trim().length>0&&item.length<=200);

/** Treat remote hierarchy files as data: bounded plain text and exact anatomy names. */
export function parseGuestHierarchy(value:unknown):GuestHierarchy{
 if(!object(value)||value.schema!=='human-atlas-hierarchy/v1'||typeof value.id!=='string'||!/^[a-z0-9][a-z0-9-]{0,63}$/.test(value.id)||typeof value.name!=='string'||!value.name.trim()||value.name.length>80||!Array.isArray(value.nodes)||value.nodes.length<1||value.nodes.length>100)throw new Error('Invalid hierarchy file.');
 if(value.source!==undefined){if(typeof value.source!=='string'||value.source.length>500)throw new Error('Invalid hierarchy source.');try{if(!['http:','https:'].includes(new URL(value.source).protocol))throw new Error();}catch{throw new Error('Invalid hierarchy source.');}}
 let count=0;
 const parseNode=(input:unknown,depth:number):GuestNode=>{
  if(!object(input)||depth>8||++count>500||typeof input.name!=='string'||!input.name.trim()||input.name.length>160)throw new Error('Invalid hierarchy node.');
  if(input.matches!==undefined&&!strings(input.matches))throw new Error('Invalid anatomy matches.');
  if(input.systems!==undefined&&!strings(input.systems,30))throw new Error('Invalid system matches.');
  if(input.children!==undefined&&(!Array.isArray(input.children)||input.children.length>100))throw new Error('Invalid hierarchy children.');
  if(input.unmapped!==undefined&&input.unmapped!==true)throw new Error('Invalid unmapped branch.');
  if(input.unmapped&&(depth!==0||input.matches||input.systems||input.children))throw new Error('The unmapped branch must be a top-level empty node.');
  return {name:input.name.trim(),...(input.matches?{matches:input.matches as string[]}:{}),...(input.systems?{systems:input.systems as SystemId[]}:{}),...(input.children?{children:(input.children as unknown[]).map(child=>parseNode(child,depth+1))}:{}),...(input.unmapped?{unmapped:true}:{})};
 };
 const nodes=value.nodes.map(node=>parseNode(node,0));
 if(nodes.filter(node=>node.unmapped).length>1)throw new Error('Only one unmapped branch is allowed.');
 return {schema:'human-atlas-hierarchy/v1',id:value.id,name:value.name.trim(),...(value.source?{source:value.source as string}:{}),nodes};
}

export function resolveGuestHierarchy(atlas:Atlas,hierarchy:GuestHierarchy):ResolvedGuestNode[]{
 const compareDepth=createDepthOrder(atlas);
 const byName=new Map<string,Part[]>(),byId=new Map(atlas.parts.map(part=>[part.id,part]));
 const concepts=new Map(atlas.concepts.map(concept=>[concept.id,concept]));
 for(const concept of atlas.concepts){const key=concept.name.toLocaleLowerCase();const members=byName.get(key)??[];for(const id of concept.elements){const part=byId.get(id);if(part)members.push(part);}byName.set(key,members);}
 const resolve=(node:GuestNode,path:string):ResolvedGuestNode=>{
  const children=(node.children??[]).map((child,index)=>resolve(child,`${path}/${index}`));
  const matched=(node.matches??[]).flatMap(term=>{
   const concept=concepts.get(term);return concept?concept.elements.map(id=>byId.get(id)).filter((part):part is Part=>!!part):byId.has(term)?[byId.get(term)!]:byName.get(term.toLocaleLowerCase())??[];
  });
  const systems=new Set(node.systems??[]);
  const directParts=[...new Map([...matched,...atlas.parts.filter(part=>systems.has(part.system))].map(part=>[part.id,part])).values()];
  const parts=[...new Map([...directParts,...children.flatMap(child=>child.parts)].map(part=>[part.id,part])).values()];
  return {id:path,name:node.name,parts,directParts,children};
 };
 const mapped=hierarchy.nodes.flatMap((node,index)=>node.unmapped?[]:[resolve(node,String(index))]);
 const assigned=new Set(mapped.flatMap(node=>node.parts.map(part=>part.id)));
 return hierarchy.nodes.map((node,index)=>{
  if(!node.unmapped)return mapped.find(item=>item.id===String(index))!;
  const parts=atlas.parts.filter(part=>!assigned.has(part.id));
  const children=SYSTEMS.flatMap(system=>{
   const members=parts.filter(part=>part.system===system.id).sort((a,b)=>compareDepth({name:a.name,parts:[a]},{name:b.name,parts:[b]}));
   if(!members.length)return [];
   return [{id:`${index}/${system.id}`,name:system.name,parts:members,directParts:[],children:members.map((part,partIndex)=>({id:`${index}/${system.id}/${partIndex}`,name:structureName(part.name),parts:[part],directParts:[part],children:[]}))}];
  });
  return {id:String(index),name:node.name,parts,directParts:[],children};
 });
}
