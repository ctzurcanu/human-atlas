import {spawnSync} from 'node:child_process';
import os from 'node:os';
import path from 'node:path';

const candidates=[process.env.PYTHON,'python3','python',path.join(os.homedir(),'anaconda3','bin','python3'),path.join(os.homedir(),'miniconda3','bin','python3')].filter(Boolean);
let interpreter;
for(const candidate of candidates){
 const check=spawnSync(candidate,['-c','import numpy'],{stdio:'ignore'});
 if(check.status===0){interpreter=candidate;break;}
}
if(!interpreter)throw new Error('NumPy-capable Python not found. Install NumPy, or set PYTHON to a Python interpreter that has it.');
const run=spawnSync(interpreter,['scripts/add-vhf-muscles.py',...process.argv.slice(2)],{stdio:'inherit'});
if(run.error)throw run.error;
process.exitCode=run.status??1;
