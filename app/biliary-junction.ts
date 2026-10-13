import type {Part} from './anatomy';

const terms=new Map([
 ['TA98PART:FJ3079:A05.8.01.061','A05.8.01.061'],
 ['TA98PART:FJ3079:A05.8.02.013','A05.8.02.013'],
]);
const femaleRegionalNames=new Map([['LOCAL:female:BILIARY:distal-assembly','Distal biliary and pancreatic ducts'],['LOCAL:female:BILIARY:A05.8.01.061','Common hepatic duct'],['LOCAL:female:BILIARY:hepatic-tributaries','Hepatic tributary ducts']]);
function member(p:Part):boolean{
 if(p.suppressed||p.system!=='digestive'||p.sourceOffset||!p.pairedTissue||!p.sectionSolid)return false;
 if(p.sectionAssembly==='gallbladder-cystic:female:folds'){
  if(p.id==='LOCAL:female:cystic_duct:authored')return p.sourceId==='cystic_duct'&&p.conceptId==='LOCAL:female:cystic_duct'&&p.name==='Cystic duct'&&p.tissueWallAssembly==='cystic:female';
  if(femaleRegionalNames.has(p.id))return p.sourceId==='bile_duct'&&p.conceptId===p.id&&p.name===femaleRegionalNames.get(p.id)&&p.tissueWallAssembly==='biliary:female';
  return p.id==='LOCAL:female:bile_duct:1031'&&p.sourceId==='bile_duct'&&p.conceptId==='LOCAL:female:bile_duct'&&p.tissueWallAssembly==='biliary:female';
 }
 if(p.sectionAssembly!=='gallbladder-cystic:male:folds')return false;
 if(p.id==='FJ3080')return p.sourceId==='FJ3080'&&p.conceptId==='FMA14539'&&p.tissueWallAssembly==='cystic:male';
 const term=terms.get(p.id);
 return !!term&&p.sourceId==='FJ3079'&&p.conceptId===`TA98BILIARY:${term}`&&p.tissueWallAssembly==='biliary:male';
}
/** Availability requires the installed shared carrier and reciprocal tissue
 * boundaries. Nearby foreign or displaced duct sources cannot inherit it. */
export function biliaryJunctionInterface(a:Part,b:Part):boolean{
 if(a.id===b.id||!member(a)||!member(b)||a.sectionAssembly!==b.sectionAssembly||a.chunk!==b.chunk||a.positions!==b.positions||a.normals!==b.normals||JSON.stringify(a.geometryOrigin)!==JSON.stringify(b.geometryOrigin))return false;
 const count=(p:Part,id:string)=>(p.tissueDrawGroups??[]).filter(g=>g.neighbor===id).reduce((n,g)=>n+g.count,0);
 return count(a,b.id)>0&&count(a,b.id)===count(b,a.id);
}
