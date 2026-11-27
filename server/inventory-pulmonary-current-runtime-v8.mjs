import assert from 'node:assert/strict';
import {readFileSync,writeFileSync,existsSync,mkdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {build} from 'esbuild';

const start=performance.now();
const base='.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1';
const version=process.argv[2]??'actual-runtime-context-v8';
assert.match(version,/^actual-runtime-context-v[89]$/);
const out=`${base}/${version}`;mkdirSync(out,{recursive:false});
const sourcePath='.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';
const bodyPath='.local-models/organ-models/spinal-cord/male-body.json';
const sha=x=>createHash('sha256').update(x).digest('hex');
const expected='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de';
const bytes=readFileSync(sourcePath),source=JSON.parse(bytes),bodyBytes=readFileSync(bodyPath),body=JSON.parse(bodyBytes);
assert.equal(sha(bytes),expected);
const audit=JSON.parse(readFileSync(`${base}/evidence.json`));
const sourceIds=[...new Set([...audit.sources.map(s=>s.id),...source.parts.filter(p=>p.id.startsWith('TA98CIRCULATION:male:')).map(p=>p.id),
  'ZA:Pulmonary trunk','ZA:Bifurcation of pulmonary trunk'])];
const fixedGuideIds=[...source.pulmonaryBranchIntegration.partIds,...source.pulmonarySegmentIntegration.partIds];
const parts=new Map(source.parts.map(p=>[p.id,p]));
const bounded=[...sourceIds,...fixedGuideIds].map(id=>parts.get(id));assert.ok(bounded.every(p=>p?.bounds));
const low=[0,1,2].map(i=>Math.min(...bounded.map(p=>p.bounds[0][i]))-.060);
const high=[0,1,2].map(i=>Math.max(...bounded.map(p=>p.bounds[1][i]))+.060);
const built=await build({stdin:{contents:`
export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities.ts';
export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations.ts';
`,resolveDir:process.cwd()},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(built.outputFiles[0].contents).toString('base64'));
const normal=api.withHiddenReviewAlternatives(api.withPreferredSourceRepresentations(api.withReviewedCatalogueIdentities(body)),false);
const runtime=new Map(normal.parts.map(p=>[p.id,p]));
// Chunk ordinals are catalogue-local. Compare the resolved bank URL separately.
const actualKeys=['vertexCount','indexCount','positions','normals','normalSigns','uvs','indices','geometryOrigin','sourceOffset','geometryScale','geometryInstance','tissueDrawGroups'];
const geometryLayout=p=>Object.fromEntries(actualKeys.filter(k=>p&&k in p).map(k=>[k,p[k]]));
const describe=p=>{
  const q=runtime.get(p.id);
  const hasSameBank=q&&source.chunks[p.chunk]?.url===normal.chunks[q.chunk]?.url;
  return {id:p.id,name:p.name,system:p.system,sourceId:p.sourceId,normalRuntimePresent:!!q,normalRuntimeSuppressed:!!q?.suppressed,
    rawSuppressed:!!p.suppressed,atlasDefaultHidden:source.defaultHidden?.includes(p.id)??false,
    ordinaryBodyDefaultHidden:normal.defaultHidden?.includes(p.id)??false,
    replacementParts:q?.replacementParts??[],retained:!!q&&!q.suppressed,
    sameGeometryLayoutAndBank:!!q&&hasSameBank&&JSON.stringify(geometryLayout(p))===JSON.stringify(geometryLayout(q)),
    sourceRecordSha256:sha(Buffer.from(JSON.stringify(p))),runtimeRecordSha256:q?sha(Buffer.from(JSON.stringify(q))):null,
    annotationType:p.referenceType??p.landmarkType??null,pairedTissue:!!p.pairedTissue,
    geometryInstance:!!p.geometryInstance,bounds_m:p.bounds,license:p.license??null,provenance:p.provenance??null,
    exclusionOnlyByExecutedLoader:true};
};
const sourceRows=sourceIds.map(id=>describe(parts.get(id)));
const context=source.parts.filter(p=>p.bounds&&p.bounds[0].every((x,i)=>x<=high[i])&&p.bounds[1].every((x,i)=>x>=low[i])).map(describe);
const result={schema:'pulmonary-family-current-runtime-source-and-context/v9',sourceCatalogueSha256:expected,ordinaryBodySha256:sha(bodyBytes),
  pipeline:['ordinary body','withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)'],
  sourceBindings:Object.fromEntries(['app/reviewed-catalogue-identities.ts','app/preferred-source-representations.ts','app/data/preferred-source-representations.json'].map(p=>[p,sha(readFileSync(p))])),
  sourceRows,fixedGuideRows:fixedGuideIds.map(id=>describe(parts.get(id))),context,
  contextBroadphaseBox_m:[low,high],declaredAdditionalEnvelopeMargin_m:.060,
  broadphaseBoxDoesNotEstablishDisplacementOrTissueClearance:true,actualCandidateDisplacementEnvelopeMustReplayInsideBox:true,
  volumeClassificationAndFullPhysicalParticipationReviewPending:true,noPhysicalSourceOrVisibilityChanges:true,
  catalogueLocalChunkOrdinalsComparedByResolvedBankURL:true,
  elapsedSeconds:(performance.now()-start)/1000};
const path=`${out}/runtime-context.json`;assert.ok(!existsSync(path));writeFileSync(path,JSON.stringify(result,null,2)+'\n');
assert.equal(sha(readFileSync(sourcePath)),expected);
console.log(JSON.stringify({path,sourceRows:sourceRows.length,retainedSources:sourceRows.filter(r=>r.retained).length,
  excludedSources:sourceRows.filter(r=>!r.retained).map(r=>r.id),contextRows:context.length,retainedContext:context.filter(r=>r.retained).length,
  guideRows:fixedGuideIds.length,retainedGuides:result.fixedGuideRows.filter(r=>r.retained).length,
  geometryMismatches:sourceRows.filter(r=>r.retained&&!r.sameGeometryLayoutAndBank).map(r=>r.id),
  outputSha256:sha(readFileSync(path)),elapsedSeconds:result.elapsedSeconds}));
