/** Compare browser screenshots captured at the same viewport as the baselines. */
import {readFile,mkdir,writeFile} from 'node:fs/promises';
import path from 'node:path';
import sharp from 'sharp';

const root=path.resolve('tests/section-visual'),current=path.resolve(process.argv[2]??path.join(root,'current'));
const cases=JSON.parse(await readFile(path.join(root,'cases.json'),'utf8'));
let failed=0;
for(const item of cases){
 const baseline=path.join(root,'baselines',item.id+'.png'),candidate=path.join(current,item.id+'.png');
 let a,b;
 try{[a,b]=await Promise.all([sharp(baseline).removeAlpha().raw().toBuffer({resolveWithObject:true}),sharp(candidate).removeAlpha().raw().toBuffer({resolveWithObject:true})]);}
 catch(error){console.error(`${item.id}: missing screenshot (${error.message})`);failed++;continue;}
 if(a.info.width!==b.info.width||a.info.height!==b.info.height){console.error(`${item.id}: viewport changed: ${a.info.width}×${a.info.height} vs ${b.info.width}×${b.info.height}`);failed++;continue;}
 const diff=Buffer.alloc(a.data.length);let total=0,changed=0;
 for(let pixel=0;pixel<a.info.width*a.info.height;pixel++){
  let max=0;for(let channel=0;channel<3;channel++){const i=pixel*3+channel,value=Math.abs(a.data[i]-b.data[i]);diff[i]=value;total+=value;max=Math.max(max,value);}
  if(max>24)changed++;
 }
 const mean=total/a.data.length/255,fraction=changed/(a.info.width*a.info.height);
 const pass=mean<=(item.maxMeanDifference??.045)&&fraction<=(item.maxChangedFraction??.09);
 console.log(`${item.id}: mean ${(mean*100).toFixed(2)}%, changed ${(fraction*100).toFixed(2)}% ${pass?'PASS':'FAIL'}`);
 if(!pass){failed++;const out=path.join(current,'diff');await mkdir(out,{recursive:true});await writeFile(path.join(out,item.id+'.png'),await sharp(diff,{raw:{width:a.info.width,height:a.info.height,channels:3}}).png().toBuffer());}
}
if(failed){console.error(`${failed} section visual case(s) differ. Review current/diff images before approving a new baseline.`);process.exitCode=1;}
