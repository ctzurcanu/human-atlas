import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,existsSync,mkdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const base='.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16';
const out=`${base}/course-contract-v21`;mkdirSync(out,{recursive:false});
const sourcePath='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json',bodyPath='.local-models/organ-models/spinal-cord/male-body.json';
const sha=x=>createHash('sha256').update(x).digest('hex'),expected='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de';
const bytes=readFileSync(sourcePath),source=JSON.parse(bytes),bodyBytes=readFileSync(bodyPath),body=JSON.parse(bodyBytes);assert.equal(sha(bytes),expected);
const proposal=JSON.parse(readFileSync(`${base}/boundary-review-v18/proposal.json`)),audit=JSON.parse(readFileSync(`${base}/skeletal-and-mouth-review-v19/evidence.json`));
const built=await build({stdin:{contents:`
export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities.ts';
export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations.ts';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const normal=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(body)),false);
const runtime=new Map(normal.parts.map(p=>[p.id,p]));
const rows=proposal.rows.map(r=>{
 const region=audit.regions.find(z=>z.side===r.side&&z.term==='A12.3.11.020');
 const low=region.sourceRegionBounds_m[0].map(x=>x-.020),high=region.sourceRegionBounds_m[1].map(x=>x+.020);
 const candidates=source.parts.filter(p=>p.bounds&&p.bounds[0].every((x,i)=>x<=high[i])&&p.bounds[1].every((x,i)=>x>=low[i]));
 return {side:r.side,sourcePartId:region.sourcePartId,motionContextBox_m:[low,high],maximumMovement_m:.010,
  rows:candidates.map(p=>{const q=runtime.get(p.id);return {id:p.id,name:p.name,system:p.system,sourceId:p.sourceId,
   normalRuntimePresent:!!q,normalRuntimeSuppressed:!!q?.suppressed,rawSuppressed:!!p.suppressed,
   atlasDefaultHidden:source.defaultHidden?.includes(p.id)??false,
   replacementParts:q?.replacementParts??[],retained:!!q&&!q.suppressed,
   annotationType:p.referenceType??p.landmarkType??null,pairedTissue:!!p.pairedTissue,
   recordSha256:sha(Buffer.from(JSON.stringify(p))),bounds_m:p.bounds,
   exclusionOnlyByExecutedLoader:true};})};
});
const result={schema:'lateral-foot-course-current-runtime-context/v21',sourceCatalogueSha256:expected,ordinaryBodySha256:sha(bodyBytes),
 pipeline:['ordinary body','withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)'],
 sourceBindings:Object.fromEntries(['app/reviewed-catalogue-identities.ts','app/preferred-source-representations.ts','app/data/preferred-source-representations.json'].map(p=>[p,sha(readFileSync(p))])),
 rows,noPhysicalSourceOrVisibilityChanges:true};
const path=`${out}/runtime-context.json`;assert.ok(!existsSync(path));writeFileSync(path,JSON.stringify(result,null,2)+'\n');
assert.equal(sha(readFileSync(sourcePath)),expected);console.log(JSON.stringify({path,rows:rows.map(r=>({side:r.side,total:r.rows.length,retained:r.rows.filter(q=>q.retained).length,excluded:r.rows.filter(q=>!q.retained).length,systems:[...new Set(r.rows.filter(q=>q.retained).map(q=>q.system))]}))}));
