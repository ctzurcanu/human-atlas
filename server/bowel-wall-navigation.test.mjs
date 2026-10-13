import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import {spawnSync} from 'node:child_process';
import {build} from 'esbuild';
const compiled=await build({entryPoints:['app/bowel-wall-review-navigation.ts'],bundle:true,platform:'node',format:'esm',write:false});
const {BOWEL_WALL_REVIEW_OPTIONS,bowelWallReviewUrl,ordinaryModelUrl}=await import('data:text/javascript;base64,'+Buffer.from(compiled.outputFiles[0].contents).toString('base64'));
for(const sex of ['female','male'])test(`${sex}: named review navigation resolves a whole bowel assembly without stale cut/selection parameters`,async()=>{
 const href=bowelWallReviewUrl('http://localhost:3016/atlas/?model=embryo&select=other&cut=axial&frame=invalid&h=old&guest=old#old',sex),url=new URL(href);
 assert.equal(url.origin,'http://localhost:3016');assert.equal(url.pathname,'/atlas/');assert.equal(url.hash,'');
 assert.equal(url.searchParams.get('model'),sex==='male'?'local-ta98':'local-female-ta98');assert.equal(url.searchParams.get('bowel'),'prototype');assert.equal(url.searchParams.get('focus'),'1');assert.equal(url.searchParams.get('tree'),'guest:ta98');
 for(const key of ['cut','frame','h','guest'])assert.equal(url.searchParams.has(key),false);
 const atlas=JSON.parse(await fs.readFile(`.local-models/ta98-review/bowel-wall-main-prototype/${sex}/atlas.json`));
 const whole=atlas.concepts.find(c=>c.id===url.searchParams.get('select'));assert.ok(whole);assert.equal(whole.elements.length,12);
 assert.ok(whole.elements.every(id=>atlas.bowelWallPrototype.partIds.includes(id)));
 const option=BOWEL_WALL_REVIEW_OPTIONS.find(o=>o.sex===sex);assert.equal(option.value,`bowel-wall-${sex}`);
 const ordinary=new URL(ordinaryModelUrl(href,sex==='male'?'local-ta98':'local-female-ta98'));
 assert.deepEqual([...ordinary.searchParams.keys()],['model']);assert.equal(ordinary.pathname,url.pathname);
});
test('build entrypoint is outside ignored development folders and available to npm',async()=>{
 const pkg=JSON.parse(await fs.readFile('package.json'));
 assert.equal(pkg.scripts['build:small-intestine-walls'],'node server/build-small-intestine-walls.mjs');
 const file=await fs.readFile('server/build-small-intestine-walls.mjs','utf8');assert.match(file,/HUMAN_ATLAS_PYTHON/);assert.match(file,/binarySha256/);
 const ignored=spawnSync('git',['check-ignore','--no-index','server/build-small-intestine-walls.mjs','server/bowel-wall-navigation.test.mjs','app/bowel-wall-review-navigation.ts'],{encoding:'utf8'});
 assert.equal(ignored.status,1,ignored.stdout||ignored.stderr);
 assert.equal(ignored.stdout,'');
});
