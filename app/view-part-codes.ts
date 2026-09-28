import type {Atlas} from './anatomy';

/** Catalogue-relative codes shrink large mesh sets. The fingerprint prevents mapping codes to a different catalogue. */
export function partFingerprint(atlas:Atlas){
 let left=2166136261,right=2246822519;
 for(const part of atlas.parts){
  for(let i=0;i<part.id.length;i++){
   const code=part.id.charCodeAt(i);
   left=Math.imul(left^code,16777619);
   right=Math.imul(right^code,2246822519);
  }
  left=Math.imul(left^0,16777619);right=Math.imul(right^0,2246822519);
 }
 return `${atlas.parts.length.toString(36)}.${(left>>>0).toString(36)}.${(right>>>0).toString(36)}`;
}
export function packPartIds(ids:string[],atlas:Atlas,ranges:boolean){
 const positions=new Map(atlas.parts.map((part,index)=>[part.id,index]));
 const numbers=ids.map(id=>positions.get(id));
 if(numbers.some(number=>number===undefined))return undefined;
 if(!ranges)return numbers.map(number=>number!.toString(36)).join('.');
 const sorted=[...new Set(numbers as number[])].sort((a,b)=>a-b),tokens:string[]=[];
 for(let i=0;i<sorted.length;){
  let end=i;while(end+1<sorted.length&&sorted[end+1]===sorted[end]+1)end++;
  if(end-i>=2){tokens.push(`${sorted[i].toString(36)}-${sorted[end].toString(36)}`);i=end+1;}
  else{tokens.push(sorted[i].toString(36));i++;}
 }
 const list=tokens.join('.');
 if(!ranges||!sorted.length)return list;
 const bytes=new Uint8Array(Math.floor(sorted.at(-1)!/8)+1);
 for(const number of sorted)bytes[number>>3]|=1<<(number&7);
 let binary='';for(const byte of bytes)binary+=String.fromCharCode(byte);
 const bits='~'+btoa(binary).replaceAll('+','-').replaceAll('/','_').replace(/=+$/,'');
 return bits.length<list.length?bits:list;
}
export function unpackPartIds(value:string,atlas:Atlas,ranges:boolean){
 if(!value||value.length>atlas.parts.length*8)return [];
 if(ranges&&value.startsWith('~')){
  const encoded=value.slice(1);if(!encoded||!/^[\w-]+$/.test(encoded))return [];
  try{
   const binary=atob(encoded.replaceAll('-','+').replaceAll('_','/'));
   if(binary.length>Math.ceil(atlas.parts.length/8))return [];
   const output:string[]=[];
   for(let byte=0;byte<binary.length;byte++)for(let bit=0;bit<8;bit++)if(binary.charCodeAt(byte)&(1<<bit)){
    const index=byte*8+bit;if(index>=atlas.parts.length)return [];output.push(atlas.parts[index].id);
   }
   return output;
  }catch{return [];}
 }
 const output:string[]=[],seen=new Set<number>();
 for(const token of value.split('.')){
  if(!/^[0-9a-z]+(?:-[0-9a-z]+)?$/.test(token))return [];
  const pair=token.split('-');if(!ranges&&pair.length>1)return [];
  const start=parseInt(pair[0],36),end=pair.length>1?parseInt(pair[1],36):start;
  if(!Number.isSafeInteger(start)||!Number.isSafeInteger(end)||start<0||end>=atlas.parts.length||end<start||output.length+end-start+1>atlas.parts.length)return [];
  for(let index=start;index<=end;index++){if(seen.has(index))return [];seen.add(index);output.push(atlas.parts[index].id);}
 }
 return output;
}
