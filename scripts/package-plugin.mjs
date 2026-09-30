import assert from 'node:assert/strict';
import {readFile,stat,mkdir} from 'node:fs/promises';
import {spawnSync} from 'node:child_process';
import {resolve} from 'node:path';
import sharp from 'sharp';
const root=resolve('plugins/human-atlas');
const manifest=JSON.parse(await readFile(`${root}/plugin.json`,'utf8'));
const {interface:listing,review}=manifest.extensions['com.openai'];
assert.match(manifest.name,/^[a-z0-9]+(?:-[a-z0-9]+)*$/);
assert.match(manifest.version,/^\d+\.\d+\.\d+$/);
for(const [field,limit] of Object.entries({displayName:30,shortDescription:30,longDescription:4000,developerName:80}))assert.ok(listing[field]?.length&&listing[field].length<=limit,field);
for(const field of ['websiteURL','supportURL','privacyPolicyURL','termsOfServiceURL'])assert.equal(new URL(listing[field]).protocol,'https:',field);
for(const prompt of listing.defaultPrompt)assert.ok(prompt.length<=128);
assert.equal(review.test_cases.positive.length,5);assert.equal(review.test_cases.negative.length,3);
for(const item of review.test_cases.positive)for(const field of ['description','prompt','tools_triggered','expected_behavior'])assert.ok(item[field]);
const config=JSON.parse(await readFile(`${root}/mcp.json`,'utf8'));
assert.equal(Object.keys(config.mcpServers).length,1);
assert.equal(config.mcpServers.atlas.type,'streamable-http');
assert.equal(config.mcpServers.atlas.url,'https://human-atlas-connect.ctzurcanu.workers.dev/mcp');
for(const file of new Set([listing.logo,listing.composerIcon,...listing.screenshots])){
 assert.ok(file.startsWith('./assets/')&&!file.includes('..'));
 const asset=`${root}/${file.slice(2)}`,details=await sharp(asset).metadata();
 assert.ok((await stat(asset)).size<=5*1024*1024);assert.ok(details.width<=4096&&details.height<=4096);
 if(file===listing.logo||file===listing.composerIcon){assert.equal(details.width,details.height);assert.ok(details.width>=48);assert.ok((await stat(asset)).size<10000);}
}
await mkdir('outputs',{recursive:true});
const output=resolve(`outputs/human-atlas-${manifest.version}.zip`);
// Include only the public package, never local configuration or credentials.
const result=spawnSync('python3',['-c',`import pathlib,zipfile,sys\nroot=pathlib.Path(sys.argv[1])\nwith zipfile.ZipFile(sys.argv[2],'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:\n for name in ['plugin.json','mcp.json']:\n  archive.write(root/name,name)\n for item in sorted((root/'assets').rglob('*')):\n  if item.is_file() and not item.name.startswith('.'):\n   archive.write(item,str(item.relative_to(root)))`,root,output],{encoding:'utf8'});
assert.equal(result.status,0,result.stderr);
console.log(`Validated and packaged: ${output} (${(await stat(output)).size} bytes)`);
console.log('Upload validation does not establish review readiness: verify deployed URLs, publisher identity, video coverage and ChatGPT test cases in the portal.');
