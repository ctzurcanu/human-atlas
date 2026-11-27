/** Read every actual resolver input and prove the two JSON consumers only read .parents. */
import {build} from 'esbuild';
import ts from 'typescript';
import {readFileSync,writeFileSync,mkdirSync,existsSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {dirname} from 'node:path';

const start=performance.now();
const sha=raw=>createHash('sha256').update(raw).digest('hex');
const output=process.argv[2];
if(!output||existsSync(output))throw Error('Supply a new, explicit local audit receipt path');
const parentPath='app/data/ta98-parents.json';
const expectedSources={male:'5f03a758aeea29e282be5a5dc077a663eeacc1c5e4141acf57b57e630ca3642a',female:'8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'};
const sourcePaths=Object.fromEntries(Object.keys(expectedSources).map(sex=>[sex,`.local-models/ta98-review/bowel-wall-local-atlas/${sex}/atlas.json`]));
for(const [sex,path]of Object.entries(sourcePaths))if(sha(readFileSync(path))!==expectedSources[sex])throw Error(`${sex} geometry master changed`);
const previousCachePath=`.local-models/geometry-workbench/dependency-relations/${expectedSources.male}/male.json`;
const previous=JSON.parse(readFileSync(previousCachePath));
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false,metafile:true});
const inputPaths=[...Object.keys(bundle.metafile.inputs),'scripts/export-structural-relation-edges.mjs'].sort();
const resolverInputs=Object.fromEntries(inputPaths.map(path=>[path,sha(readFileSync(path))]));
if(JSON.stringify(Object.keys(previous.resolverInputs).sort())!==JSON.stringify(inputPaths))throw Error('Resolver input set changed; requires explicit new review');
const changedInputs=inputPaths.filter(path=>previous.resolverInputs[path]!==resolverInputs[path]);
const allowed=new Set(['app/anatomical-relations.ts','app/gallbladder-relations.ts','app/data/ta98-metadata.json',parentPath]);
if(changedInputs.some(path=>!allowed.has(path)))throw Error(`Unexpected resolver change: ${changedInputs.filter(path=>!allowed.has(path)).join(', ')}`);
const consumers=Object.entries(bundle.metafile.inputs).filter(([,info])=>info.imports.some(item=>item.path===parentPath)).map(([path])=>path).sort();
const expectedConsumers=['app/anatomical-terminology.ts','app/ta98-modeled-parent.ts'];
if(JSON.stringify(consumers)!==JSON.stringify(expectedConsumers))throw Error('Parent-index consumer set changed');
const consumerProofs=consumers.map(path=>{
 const text=readFileSync(path,'utf8'),tree=ts.createSourceFile(path,text,ts.ScriptTarget.Latest,true,ts.ScriptKind.TS);
 const imports=tree.statements.filter(node=>ts.isImportDeclaration(node)&&ts.isStringLiteral(node.moduleSpecifier)&&node.moduleSpecifier.text==='./data/ta98-parents.json');
 if(imports.length!==1||!imports[0].importClause?.name||imports[0].importClause.namedBindings)throw Error(`Unreviewed parent-index import form: ${path}`);
 const binding=imports[0].importClause.name.text,uses=[];
 const visit=node=>{
  if(ts.isIdentifier(node)&&node.text===binding){
   if(node===imports[0].importClause.name)return;
   if(!ts.isPropertyAccessExpression(node.parent)||node.parent.expression!==node||node.parent.name.text!=='parents')throw Error(`Parent index escapes .parents-only use in ${path}`);
   const position=tree.getLineAndCharacterOfPosition(node.getStart(tree));
   uses.push({line:position.line+1,column:position.character+1,expression:text.slice(node.parent.pos,node.parent.end).trim()});
  }
  ts.forEachChild(node,visit);
 };
 visit(tree);
 if(uses.length!==1)throw Error(`Unreviewed number of .parents accesses: ${path}`);
 return {path,sha256:resolverInputs[path],binding,onlyProperty:'parents',uses,typescriptASTVerified:true};
});
const parent=JSON.parse(readFileSync(parentPath));
if(JSON.stringify(Object.keys(parent).sort())!==JSON.stringify(['hierarchySha256','parents','source']))throw Error('Parent index schema changed');
const signature=sha(JSON.stringify(Object.fromEntries(Object.entries(resolverInputs).filter(([path])=>path!==parentPath))));
const report={schema:'relationship-parent-only-input-audit/v1',sourceCatalogueSha256:expectedSources,sourcePaths,resolverInputs,inputCount:inputPaths.length,previousCachePath,previousCacheSha256:sha(readFileSync(previousCachePath)),changedInputs,parentIndexConsumers:consumerProofs,parentOnlyInputSignature:signature,parentIndexShape:Object.keys(parent).sort(),parentIndexSha256:resolverInputs[parentPath],noGeometryOrCacheModified:true,elapsedSeconds:(performance.now()-start)/1000};
mkdirSync(dirname(output),{recursive:true});writeFileSync(output,JSON.stringify(report)+'\n',{flag:'wx'});
console.log(JSON.stringify({output,inputCount:inputPaths.length,changedInputs,consumers,parentOnlyInputSignature:signature,elapsedSeconds:report.elapsedSeconds}));
