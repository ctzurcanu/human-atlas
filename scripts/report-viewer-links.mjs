import {existsSync,readFileSync} from 'node:fs';

export const REPORT_VIEWER='http://localhost:3016/';
export const REPORT_MODELS=[
 ['male-detail','public/models/atlas-male-complete.json'],
 ['male','public/models/atlas.json'],['female','public/models/atlas-hra-female.json'],
 ['embryo','public/models/atlas-embryo.json'],['cell','public/models/atlas-cell.json'],
 ['local-reference','.local-models/reference.json'],['local-male','.local-models/male.json'],['local-female','.local-models/female.json'],
];
export const reportModel=model=>model==='male-detail / male-full'||model==='legacy Z-Anatomy'?'male-detail':REPORT_MODELS.some(([id])=>id===model)?model:undefined;
export const markdownLabel=value=>String(value).replace(/[\\[\]`*_]/g,'\\$&').replaceAll('|','\\|').replaceAll('\n',' ');
export function reportViewUrl(model,ids=[],tree='systems'){
 const url=new URL(REPORT_VIEWER);url.searchParams.set('model',model);
 if(tree!=='systems')url.searchParams.set('tree',tree);
 for(const id of ids)url.searchParams.append('select',id);
 if(ids.length){url.searchParams.set('context','0.18');url.searchParams.set('skin','0.1');url.searchParams.set('focus','1');}
 return url.href;
}

/** Choose a supported viewer model containing the actual available concept. */
export function createReportViewerLinks(){
 const byId=new Map(),byName=new Map(),sourceNames=new Map(),atlases=new Map();
 const add=(map,key,entry)=>map.set(key,[...(map.get(key)??[]),entry]);
 for(const [model,file] of REPORT_MODELS){
  if(!existsSync(file))continue;
  const atlas=JSON.parse(readFileSync(file,'utf8'));atlases.set(model,atlas);
  const available=new Map(atlas.parts.filter(part=>!part.suppressed).map(part=>[part.id,part]));
  for(const concept of atlas.concepts){
   sourceNames.set(concept.id,concept.name);
   if(!concept.elements.some(id=>available.has(id)))continue;
   const entry={model,id:concept.id,name:concept.name};
   add(byId,concept.id,entry);add(byName,concept.name.toLowerCase(),entry);
  }
  for(const part of available.values())add(byId,part.id,{model,id:part.id,name:part.name});
 }
 const choose=(candidates,models=[])=>candidates.find(item=>[...models].map(reportModel).includes(item.model))??candidates[0];
 const find=({id,name,models=[]})=>choose(byId.get(id)??[],models)??choose(byName.get((sourceNames.get(id)??name??'').toLowerCase())??[],models);
 const link=item=>{
  const found=find(item),label=markdownLabel(item.name);
  if(found)return `[${label}](${reportViewUrl(found.model,[found.id])})`;
  const model=[...(item.models??[])].map(reportModel).find(model=>atlases.has(model));
  return model?`[${label}](${reportViewUrl(model)}) (suppressed; opens model)`:label;
 };
 return {find,link,atlases,byName};
}
