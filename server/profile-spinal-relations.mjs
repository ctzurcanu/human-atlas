/** Read-only timing of the current relationship resolver on one staged atlas. */
import {readFileSync} from 'node:fs';
import {build} from 'esbuild';
const file=process.argv[2];if(!file)throw Error('Explicit candidate catalogue required');
let source=readFileSync('app/anatomical-relations.ts','utf8');
const names=['vascularWallRelations','caecocolicWallRelations','smallBowelRelations','colonicTaeniaRelations','ileocaecalReferenceRelations','ileocaecalWallRelations','endocrineRelations','thyroidFunctionalRelations','digestiveDuctRelations','gallbladderRelations','cysticFunctionalRelations','gallbladderWallRelations','foregutArterialRelations','hepaticRelations','gastricVenousRelations','structuralNetwork','skeletalNetwork','lymphaticNetwork','renalVenousNetwork','aorticBranchNetwork','adrenalNetwork','reviewedNativeRoutes','cranialRelations','cardiacRelations','spinalCordCoveringRelations','pulmonaryRelations','laryngealRelations','coronaryRelations','fibularRelations','handVascularRelations','nativeUrinaryRelations','analSphincterRelations','rectalMuscleRelations'];
for(const name of names)source=source.replaceAll(`${name}(parts)`,`traceProfile('${name}',()=>${name}(parts))`);
source+="\nfunction traceProfile<T>(name:string,run:()=>T):T{const start=performance.now();try{return run();}finally{console.log('[profile]',name,(performance.now()-start).toFixed(2));}}";
const bundle=await build({stdin:{contents:source,resolveDir:process.cwd()+'/app',loader:'ts'},bundle:true,platform:'node',format:'esm',write:false});
const {anatomicalRelations}=await import('data:text/javascript;base64,'+Buffer.from(bundle.outputFiles[0].contents).toString('base64'));
const a=JSON.parse(readFileSync(file)),p=a.parts.find(p=>!p.suppressed);
console.log('[profile] first part',p.id,'catalogue parts',a.parts.length);
const start=performance.now();const result=anatomicalRelations(p,a.parts);
console.log('[profile] complete',result.length,'milliseconds',performance.now()-start);
