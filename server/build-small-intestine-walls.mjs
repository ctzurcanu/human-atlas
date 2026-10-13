// Regenerate the local first version from the retained, inventoried wall carriers.
import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import crypto from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {spawnSync} from 'node:child_process';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
try{
 const required=['scripts/stage-small-intestine-fold-relief.py','scripts/bowel_distal_orientation.py','scripts/pack-small-intestine-prototype.mjs','scripts/bowel-feature-partition.mjs','web/bowel-fold-fields.ts','.local-models/male.json','.local-models/ta98-male.json','.local-models/ta98-female.json','.local-models/ta98-runtime.json','.local-models/ta98-female-runtime.json'];
 for(const sex of ['female','male']){
  const base='.local-models/ta98-review/small-intestine-layer-v1'+(sex==='male'?'-male':'');
  required.push(base+'/manifest.json',base+'/wall.bin');
 }
 const missing=required.filter(p=>!fs.existsSync(path.join(root,p)));
 if(missing.length)throw Error('Local wall sources are missing:\n'+missing.join('\n')+'\nRestore the retained small-intestine review carriers first. These local assets are intentionally outside deployment.');
 for(const sex of ['female','male']){
  const base=path.join(root,'.local-models/ta98-review/small-intestine-layer-v1'+(sex==='male'?'-male':''));
  const meta=JSON.parse(fs.readFileSync(base+'/manifest.json'));
  if(sha(fs.readFileSync(base+'/wall.bin'))!==meta.binarySha256)throw Error(sex+' retained wall carrier hash differs from its manifest.');
 }
 const bundled=path.join(os.homedir(),'.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3');
 const candidates=process.env.HUMAN_ATLAS_PYTHON?[process.env.HUMAN_ATLAS_PYTHON]:['python3',...(fs.existsSync(bundled)?[bundled]:[])];
 const python=candidates.find(command=>spawnSync(command,['-c','import numpy'],{cwd:root,stdio:'ignore'}).status===0);
 if(!python)throw Error('Python with NumPy is required. Set HUMAN_ATLAS_PYTHON to an existing interpreter with NumPy. No packages have been installed.');
 if(process.argv.includes('--check')){console.log('Both retained wall carriers and a NumPy interpreter are available.');process.exit(0);}
 for(const sex of ['female','male']){
  console.log('Building '+sex+' small-intestine wall review…');
  const stage=spawnSync(python,['scripts/stage-small-intestine-fold-relief.py','--model',sex],{cwd:root,encoding:'utf8'});
  if(stage.status!==0)throw Error(stage.stderr||stage.error?.message||'Physical fold field generation failed.');
  const pack=spawnSync(process.execPath,['scripts/pack-small-intestine-prototype.mjs',...(sex==='male'?['--male']:[])],{cwd:root,stdio:'inherit'});
  if(pack.status!==0)throw Error('The '+sex+' wall catalogue could not be packed.');
  const local=spawnSync(process.execPath,['scripts/pack-small-intestine-prototype.mjs','--ordinary-local',...(sex==='male'?['--male']:[])],{cwd:root,stdio:'inherit'});
  if(local.status!==0)throw Error('The '+sex+' local atlas layers could not be packed.');
 }
 if(fs.existsSync(path.join(root,'.local-models/ta98-review/native-taenia-footprints/manifest.json'))){const refs=spawnSync(python,['scripts/install-local-taenia-footprints.py'],{cwd:root,stdio:'inherit'});if(refs.status!==0)throw Error('Native taenia display references no longer match the source; rebuild their footprints.');}
 console.log('Both local wall review models are ready in the anatomy model menu. Anatomy remains provisional.');
}catch(error){console.error(error.message);process.exitCode=1;}
