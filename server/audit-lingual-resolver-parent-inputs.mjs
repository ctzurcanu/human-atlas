/** Bind the final resolver import graph and prove the parent-index read scope. */
import {build} from 'esbuild';
import ts from 'typescript';
import {readFileSync,writeFileSync,mkdirSync,existsSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {dirname} from 'node:path';

const began=performance.now(),hash=raw=>createHash('sha256').update(raw).digest('hex');
const output=process.argv[2];
if(!output||existsSync(output))throw Error('Supply a fresh explicit audit receipt');
const sources={male:{path:'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json',sha256:'5c29865867198224cf57150477a9bff3e577c1ec58f87b13622cd654c9f4a345'},
 female:{path:'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json',sha256:'8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e'},
 prospectiveMale:{path:'.local-models/geometry-workbench/circulatory-completion/male-cystic-artery-correction-proposal-v43/proposed-atlas.json',sha256:'b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69'}};
for(const source of Object.values(sources))if(hash(readFileSync(source.path))!==source.sha256)throw Error('Source changed: '+source.path);
const parentPath='app/data/ta98-parents.json';
const bundle=await build({entryPoints:['app/anatomical-relations.ts'],bundle:true,platform:'node',format:'esm',write:false,metafile:true});
const inputPaths=[...Object.keys(bundle.metafile.inputs),'scripts/export-structural-relation-edges.mjs'].sort();
const resolverInputs=Object.fromEntries(inputPaths.map(path=>[path,hash(readFileSync(path))]));
const consumers=Object.entries(bundle.metafile.inputs).filter(([,info])=>info.imports.some(item=>item.path===parentPath)).map(([path])=>path).sort();
if(JSON.stringify(consumers)!==JSON.stringify(['app/anatomical-terminology.ts','app/ta98-modeled-parent.ts']))throw Error('Unexpected parent-index consumer');
const proofs=consumers.map(path=>{
 const text=readFileSync(path,'utf8'),tree=ts.createSourceFile(path,text,ts.ScriptTarget.Latest,true,ts.ScriptKind.TS);
 const imports=tree.statements.filter(n=>ts.isImportDeclaration(n)&&ts.isStringLiteral(n.moduleSpecifier)&&n.moduleSpecifier.text==='./data/ta98-parents.json');
 if(imports.length!==1||!imports[0].importClause?.name||imports[0].importClause.namedBindings)throw Error('Unexpected parent-index import form');
 const binding=imports[0].importClause.name.text,uses=[],aliases=[];
 const place=n=>{const p=tree.getLineAndCharacterOfPosition(n.getStart(tree));return {line:p.line+1,column:p.character+1};};
 const visit=node=>{
  if(ts.isIdentifier(node)&&node.text===binding&&node!==imports[0].importClause.name){
   const p=node.parent;
   if(!ts.isPropertyAccessExpression(p)||p.expression!==node||p.name.text!=='parents')throw Error('Parent index escapes .parents-only access');
   let value=p;
   while(ts.isAsExpression(value.parent)||ts.isTypeAssertionExpression(value.parent)||ts.isParenthesizedExpression(value.parent))value=value.parent;
   const declaration=value.parent;
   if(!ts.isVariableDeclaration(declaration)||declaration.initializer!==value||!ts.isIdentifier(declaration.name))throw Error('Parent map has an unreviewed alias form');
   uses.push({...place(node),expression:p.getText(tree)});aliases.push(declaration.name);
  }
  ts.forEachChild(node,visit);
 };
 visit(tree);
 if(uses.length!==1||aliases.length!==1)throw Error('Unexpected number of parent map accesses');
 const alias=aliases[0],reads=[];
 const scan=node=>{
  if(ts.isIdentifier(node)&&node.text===alias.text&&node!==alias&&!(ts.isPropertyAccessExpression(node.parent)&&node.parent.name===node)){
   const p=node.parent;
   if(!ts.isElementAccessExpression(p)||p.expression!==node)throw Error('Parent map alias escapes element reads');
   const context=p.parent;
   if(ts.isDeleteExpression(context)||ts.isPrefixUnaryExpression(context)||ts.isPostfixUnaryExpression(context)
    ||ts.isBinaryExpression(context)&&context.left===p&&context.operatorToken.kind>=ts.SyntaxKind.FirstAssignment&&context.operatorToken.kind<=ts.SyntaxKind.LastAssignment)throw Error('Parent map is mutated');
   reads.push({...place(node),expression:p.getText(tree)});
  }
  ts.forEachChild(node,scan);
 };
 scan(tree);
 if(!reads.length)throw Error('No reviewed parent map reads');
 return {path,sha256:resolverInputs[path],binding,onlyProperty:'parents',uses,alias:alias.text,aliasElementReads:reads,typescriptASTVerified:true,noIndexMetadataRead:true,noParentMapMutation:true};
});
const parent=JSON.parse(readFileSync(parentPath));
if(JSON.stringify(Object.keys(parent).sort())!==JSON.stringify(['hierarchySha256','parents','source']))throw Error('Parent-index schema changed');
const previousPath='.local-models/geometry-workbench/circulatory-completion/current-relationship-finalization-v1/resolver-input-audit.json';
const previous=JSON.parse(readFileSync(previousPath));
const changed=Object.keys(resolverInputs).filter(p=>previous.resolverInputs[p]!==resolverInputs[p]);
const removed=Object.keys(previous.resolverInputs).filter(p=>!Object.hasOwn(resolverInputs,p));
const parentOnlyInputSignature=hash(JSON.stringify(Object.fromEntries(Object.entries(resolverInputs).filter(([p])=>p!==parentPath))));
const report={schema:'relationship-parent-only-input-audit/v2',sources,resolverInputs,inputCount:inputPaths.length,parentIndexConsumers:proofs,parentOnlyInputSignature,
 parentIndexShape:Object.keys(parent).sort(),parentIndexSha256:resolverInputs[parentPath],previousAudit:previousPath,previousAuditSha256:hash(readFileSync(previousPath)),changedInputs:changed,removedInputs:removed,
 reviewedChangeScope:'Final blood-only lingual functional relations, installed arterial frame/component controls and already reviewed cystic recipient repair; every actual resolver input remains hash-bound.',
 noGeometryOrCacheModified:true,elapsedSeconds:(performance.now()-began)/1000};
for(const [path,expected] of Object.entries(resolverInputs))if(hash(readFileSync(path))!==expected)throw Error('Input changed during audit');
mkdirSync(dirname(output),{recursive:true});writeFileSync(output,JSON.stringify(report)+'\n',{flag:'wx'});
console.log(JSON.stringify({output,inputCount:inputPaths.length,consumers,parentOnlyInputSignature,changedInputs:changed,removedInputs:removed,elapsedSeconds:report.elapsedSeconds}));
