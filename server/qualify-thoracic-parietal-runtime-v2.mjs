import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,existsSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';
const base='.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1';
const cat='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json',body='.local-models/organ-models/spinal-cord/male-body.json';
const sha=x=>createHash('sha256').update(x).digest('hex'),expected='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de';
assert.equal(sha(readFileSync(cat)),expected);const source=JSON.parse(readFileSync(cat)),rawBody=JSON.parse(readFileSync(body));
const bundle=await build({stdin:{contents:`export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities.ts';export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations.ts';`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64')),normal=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(rawBody)),false),np=new Map(normal.parts.map(p=>[p.id,p]));
const low=[-.075,1.285,-.035],high=[.075,1.55,.145];
const candidates=source.parts.filter(p=>p.bounds&&p.bounds[0].every((x,i)=>x<=high[i])&&p.bounds[1].every((x,i)=>x>=low[i]));
const rows=candidates.map(p=>({id:p.id,name:p.name,system:p.system,sourceId:p.sourceId,bounds_m:p.bounds,sourceRecordSha256:sha(Buffer.from(JSON.stringify(p))),
    ordinaryBodyPresent:np.has(p.id),normalSuppressed:!!np.get(p.id)?.suppressed,normalAvailable:np.has(p.id)&&!np.get(p.id)?.suppressed,
    originalSuppressed:!!p.suppressed,atlasDefaultHidden:source.defaultHidden?.includes(p.id)??false,replacementParts:np.get(p.id)?.replacementParts??[]}));
const actualExplicitTargets=source.parts.filter(p=>p.id.toLowerCase().includes('phrenic_nerve')||p.id.toLowerCase().includes('pericardium')||p.id.toLowerCase().includes('pericardiacophrenic')||p.id.toLowerCase().includes('subclavian artery')).map(p=>({id:p.id,system:p.system,normalAvailable:np.has(p.id)&&!np.get(p.id)?.suppressed,sourceRecordSha256:sha(Buffer.from(JSON.stringify(p))),bounds_m:p.bounds}));
const path=`${base}/actual-runtime-context-v2.json`;assert(!existsSync(path));writeFileSync(path,JSON.stringify({schema:'thoracic-parietal-current-runtime-context/v2',sourceCatalogueSha256:expected,ordinaryBodySha256:sha(readFileSync(body)),contextBox_m:[low,high],pipeline:['ordinary body','withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)'],rows,actualExplicitTargets,
codeBindings:Object.fromEntries(['app/reviewed-catalogue-identities.ts','app/preferred-source-representations.ts','app/data/preferred-source-representations.json'].map(p=>[p,sha(readFileSync(p))])),noSourceOrVisibilityChanges:true},null,2)+'\n');
console.log(JSON.stringify({path,contextRows:rows.length,normalAvailable:rows.filter(r=>r.normalAvailable).length,explicitTargets:actualExplicitTargets}));
