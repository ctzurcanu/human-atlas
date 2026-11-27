/** Evaluate actual current runtime skin/preference/visibility policy; no UI edits. */
import fs from 'node:fs';
import path from 'node:path';
import assert from 'node:assert/strict';
import crypto from 'node:crypto';
import {build} from 'esbuild';
const started=performance.now(),root=process.cwd();
const masterPath=path.resolve('.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json');
const bodyPath=path.resolve('.local-models/organ-models/spinal-cord/male-body.json');
const out=path.resolve('.local-models/geometry-workbench/biodigital-male-circulation-v1/global-envelope-v1/fixed-skin');
assert.ok(!fs.existsSync(out),'Refusing to overwrite fixed-skin');
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const masterBytes=fs.readFileSync(masterPath),bodyBytes=fs.readFileSync(bodyPath),master=JSON.parse(masterBytes),body=JSON.parse(bodyBytes);
assert.equal(sha(masterBytes),'46262d8ab07682fcae3ff0ee634f652abad766c1e010f6e2ec311755aea9abf5');
const bundled=await build({stdin:{contents:`export {withReviewedCatalogueIdentities} from './app/reviewed-catalogue-identities.ts';export {withPreferredSourceRepresentations,withHiddenReviewAlternatives} from './app/preferred-source-representations.ts';export {isSkinPart} from './app/depth-layers.ts';export {readViewUrl,partVisible} from './app/viewer-state.ts';export {SYSTEMS} from './app/anatomy.ts';`,resolveDir:root,sourcefile:'fixed-skin-policy.ts'},bundle:true,platform:'node',format:'esm',write:false,logLevel:'silent'});
const api=await import('data:text/javascript;base64,'+Buffer.from(bundled.outputFiles[0].contents).toString('base64'));
const reviewed=api.withReviewedCatalogueIdentities(body),preferred=api.withPreferredSourceRepresentations(reviewed),standard=api.withHiddenReviewAlternatives(preferred,false),hidden=new Set(standard.defaultHidden??[]);
const state=api.readViewUrl('?model=local-ta98',standard,{focus:0,contextOpacity:1,region:'all',peel:0,depthHidden:[],hidden:[],labels:true,explode:0,visible:api.SYSTEMS.map(s=>s.id),selected:[],isolate:false,view:'three-quarter',rotate:false,reset:0,skinOpacity:1});
const rawById=new Map(body.parts.map(p=>[p.id,p])),masterById=new Map(master.parts.map(p=>[p.id,p]));
const decisions=standard.parts.filter(p=>api.isSkinPart(p)||p.system==='integumentary').map(p=>{
 const raw=rawById.get(p.id),activeLayer=api.isSkinPart(p)&&!p.suppressed&&!p.defaultHidden&&!hidden.has(p.id)&&api.partVisible(p,state,standard.parts);
 const bodySurfaceName=api.isSkinPart(p)&&/\bskin\b|^body surface\b/i.test(p.name);
 const declaredHair=(p.groups??[]).some(g=>/^hairs?$/i.test(g))||p.id==='ZA:Scalp hair:metric';
 const declaredNailPlate=(p.groups??[]).some(g=>/^nail plate$/i.test(g));
 const primary=activeLayer&&!declaredHair&&!declaredNailPlate;
 if(activeLayer){const original=masterById.get(p.id);assert.ok(original,p.id);assert.deepEqual({...raw,chunk:body.sourceChunkIndices[raw.chunk]},original,`Body source must replay current master ${p.id}`);}
 return {id:p.id,name:p.name,system:p.system,isSkinPart:api.isSkinPart(p),enclosingBodyPredicate:bodySurfaceName,activeSkinLayer:activeLayer,declaredHair,declaredNailPlate,activePrimarySkin:primary,activeBodySurface:primary,activeSkinLayerAccessory:activeLayer&&!primary,rawSuppressed:!!raw?.suppressed,effectiveSuppressed:!!p.suppressed,partDefaultHidden:!!p.defaultHidden,catalogueDefaultHidden:hidden.has(p.id),runtimePartVisible:api.partVisible(p,state,standard.parts),replacementParts:p.replacementParts??[]};
});
assert.ok(!decisions.find(p=>p.id==='FJ2810')?.activePrimarySkin,'Hidden FJ skin must be excluded');
const files=['app/depth-layers.ts','app/preferred-source-representations.ts','app/reviewed-catalogue-identities.ts','app/viewer-state.ts','app/data/preferred-source-representations.json'];
const policy={schema:'current-primary-skin-runtime-policy/v1',masterPath,masterSha256:sha(masterBytes),bodyPath,bodySha256:sha(bodyBytes),policyCodeHashes:Object.fromEntries(files.map(f=>[f,sha(fs.readFileSync(f))])),actualRuntimeFunctions:true,rawSelectedBodyRecordsReplayCurrentMaster:true,bodySurfacePredicate:'Current active isSkinPart set, excluding only source-declared Hairs and Nail plate appendages plus the active metric scalp-hair accessory. Actual regional skin such as Sole, auricle, eyebrow skin and perionyx remains included.',opacityNamePredicateIsNotCompleteBodySkinClassifier:true,selectionScope:'Current preferred ordinary body; skin layer enabled, whole-body view, no explicit selections or user hides. All catalogue/part defaultHidden and suppression flags enforced.',decisions,counts:{activeSkinLayer:decisions.filter(p=>p.activeSkinLayer).length,activePrimarySkin:decisions.filter(p=>p.activePrimarySkin).length,activeBodySurface:decisions.filter(p=>p.activeBodySurface).length,activeAccessories:decisions.filter(p=>p.activeSkinLayerAccessory).length,excludedHiddenSkinParts:decisions.filter(p=>p.isSkinPart&&!p.activeSkinLayer).length},elapsedSeconds:(performance.now()-started)/1000};
fs.mkdirSync(out,{recursive:true});fs.writeFileSync(path.join(out,'policy.json'),JSON.stringify(policy,null,2)+'\n');
console.log(JSON.stringify({output:out,counts:policy.counts,activeAccessories:decisions.filter(p=>p.activeSkinLayerAccessory).map(p=>({id:p.id,name:p.name})),excluded:decisions.filter(p=>p.isSkinPart&&!p.activePrimarySkin).map(p=>({id:p.id,name:p.name})),elapsedSeconds:policy.elapsedSeconds}));
