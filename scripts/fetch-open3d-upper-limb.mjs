/** Fetch and verify the original Open 3D Model upper-limb GLB archive. */
import fs from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {spawn} from 'node:child_process';
import path from 'node:path';
import os from 'node:os';

const out=path.resolve(process.argv[2]??path.join(os.tmpdir(),'human-atlas-primary-models'));
const zip=path.join(out,'upper-limb-glb.zip');
const glb=path.join(out,'upper-limb.glb');
const expectedZip='5af0190a6d7bf47393447ac30021e4f3ba619721c7f3a620c39a895947078432';
const expectedGlb='e440c84c794239d1850e62b4ede0195d81bcff4f0078c7945528c388ab72fdb4';
const url='https://caskanatomy.info/open3dmodelfiles/upper-limb/upper-limb-glb.zip';
async function digest(file){try{return createHash('sha256').update(await fs.readFile(file)).digest('hex')}catch{return null}}
async function run(command,args,output){const child=spawn(command,args,{stdio:['ignore',output??'inherit','inherit']});const code=await new Promise((resolve,reject)=>{child.once('error',reject);child.once('close',resolve)});if(code!==0)throw new Error(`${command} exited with ${code}`)}
await fs.mkdir(out,{recursive:true});
if(await digest(zip)!==expectedZip){
 const temp=zip+'.download';
 await run('curl',['--fail','--location','--silent','--show-error','--output',temp,url]);
 if(await digest(temp)!==expectedZip)throw new Error('Open 3D Model publisher ZIP hash changed');
 await fs.rename(temp,zip);
}
if(await digest(glb)!==expectedGlb){
 const handle=await fs.open(glb+'.download','w');
 try{await run('unzip',['-p',zip,'upper-limb.glb'],handle.fd)}finally{await handle.close()}
 if(await digest(glb+'.download')!==expectedGlb)throw new Error('Open 3D Model publisher GLB hash changed');
 await fs.rename(glb+'.download',glb);
}
console.log(`Verified publisher Open 3D Model upper-limb GLB: ${glb}`);
