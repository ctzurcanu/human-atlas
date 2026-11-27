/** At most 15 s before starting another internally bounded ordinal batch. */
import {readFileSync,existsSync} from 'node:fs';
import {resolve,join} from 'node:path';
import {createHash} from 'node:crypto';
import {run,verifyRecord} from '../scripts/export-structural-relation-edges-staged.mjs';
const option=name=>{const i=process.argv.indexOf(name);if(i<0)throw Error('Missing '+name);return process.argv[i+1];};
const sex=option('--sex'),start=Number(option('--from')),end=Number(option('--through'));
const input=option('--catalogue-file'),dir=resolve(option('--out-dir'));
const read=path=>JSON.parse(readFileSync(path)),hash=raw=>createHash('sha256').update(raw).digest('hex');
const manifestRaw=readFileSync(join(dir,'manifest.json')),manifest=JSON.parse(manifestRaw),signature=hash(manifestRaw);
const audit=read('.local-models/geometry-workbench/circulatory-completion/lingual-controls-relation-cache-v1/resolver-input-audit.json');
if(JSON.stringify(manifest.resolverInputs)!==JSON.stringify(audit.resolverInputs))throw Error('Unreviewed resolver inputs');
if(manifest.sex!==sex||manifest.input!==resolve(input)||manifest.batchSize!==500||start<0||end>=manifest.batchCount||start>end)throw Error('Invalid source-bound batch range');
const began=performance.now(),completed=[];
for(let index=start;index<=end;index++){
 if(performance.now()-began>=15000)break;
 await run(['--sex',sex,'--stage','batch','--batch-index',String(index),'--catalogue-file',input,'--out-dir',dir]);
 const path=join(dir,`batch-${String(index).padStart(4,'0')}.json`);
 if(!existsSync(path)){
  const progress=read(join(dir,`progress-${String(index).padStart(4,'0')}.json`));
  verifyRecord(progress,{schema:'exact-anatomical-relation-batch/v1',manifestSha256:signature,index,complete:false});
  console.log(JSON.stringify({groupStoppedAtIncomplete:index,verifiedDurableThrough:progress.completedThrough,seconds:(performance.now()-began)/1000}));
  break;
 }
 const batch=read(path);
 verifyRecord(batch,{schema:'exact-anatomical-relation-batch/v1',manifestSha256:signature,index,complete:true,completedThrough:Math.min((index+1)*500,manifest.partCount)});
 completed.push(index);
}
console.log(JSON.stringify({sex,completed,seconds:(performance.now()-began)/1000}));
