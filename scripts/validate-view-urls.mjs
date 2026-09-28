import assert from 'node:assert/strict';
import {build} from 'esbuild';
import {gzipSync,strToU8} from 'fflate';
import {readFileSync} from 'node:fs';

const bundle=await build({entryPoints:['app/view-url.ts','app/viewer-state.ts','app/embed.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/view-url-validation'});
const load=async file=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(output=>output.path.endsWith(file)).text).toString('base64')}`);
const {compactViewUrl,viewParameters}=await load('view-url.js');
const {viewUrl,readViewUrl}=await load('viewer-state.js');
const {embedUrl}=await load('embed.js');
const base='https://ctzurcanu.github.io/human-atlas/';
const short=base+'?model=male-detail&select=ZA%3AStomach';
assert.equal(compactViewUrl(short),short);

// Reproduce the failing link: thousands of hidden named anatomical pieces.
const parts=Array.from({length:3756},(_,i)=>({id:`ZA:Hidden anatomical structure ${i}.l`,name:`Structure ${i}`,system:'muscular'}));
const atlas={parts,concepts:[]};
const camera=[1,2,3,0,.9,0,.1,0];
const state={selected:[parts[0].id,parts[2].id],hidden:parts.map(part=>part.id),visible:['muscular'],view:'front',contextOpacity:.18,skinOpacity:.1,guestQuery:'CL:0000185',region:'all',section:{enabled:true,axis:'oblique',position:.41,flip:true,azimuth:35,elevation:30},labels:true};
const url=viewUrl(base,'male-detail',state,camera),parsed=new URL(url);
assert.ok(parsed.hash.startsWith('#atlas-view=1.'));
assert.ok(parsed.origin.length+parsed.pathname.length+parsed.search.length<4096);
assert.ok(url.length<100000);
assert.equal(parsed.searchParams.has('hide'),false);
const restored=readViewUrl(url,atlas,{});
for(const key of ['selected','hidden','visible','view','contextOpacity','skinOpacity','guestQuery','region','section','labels'])assert.deepEqual(restored[key],state[key]);
assert.deepEqual(restored.camera,camera);

// Adding hierarchy/iframe controls must retain the camera and every hidden item.
parsed.searchParams.set('tree','guest:cell-types');parsed.searchParams.append('guest','https://example.org/custom.json');
const hierarchyUrl=compactViewUrl(parsed.href);
const embedded=embedUrl(hierarchyUrl,['study','systems']);
const params=viewParameters(embedded);
assert.equal(params.get('tree'),'guest:cell-types');assert.equal(params.get('embed'),'1');assert.equal(params.get('ui'),'study,systems');
assert.deepEqual(params.getAll('guest'),['https://example.org/custom.json']);
assert.deepEqual(readViewUrl(embedded,atlas,{}).hidden,state.hidden);
assert.deepEqual(readViewUrl(embedded,atlas,{}).camera,camera);
assert.equal(compactViewUrl(hierarchyUrl),hierarchyUrl);

// Existing long bookmarks/slides can be migrated without knowing a model's IDs.
const legacy=new URL(base);legacy.search=viewParameters(hierarchyUrl).toString();
assert.ok(legacy.href.length>150000);
assert.deepEqual([...viewParameters(compactViewUrl(legacy.href))],[...legacy.searchParams]);
assert.equal(viewParameters(base+'?model=male#atlas-view=1.invalid').get('model'),'male');
assert.equal(viewParameters(base+'#unrelated-anchor').size,0);
const oversize=gzipSync(strToU8('x'.repeat(2_000_001)));
assert.equal(viewParameters(base+'#atlas-view=1.'+Buffer.from(oversize).toString('base64url')).size,0);

// Real anatomical IDs are long. Numeric mesh codes cut an otherwise 26 KB link to a few hundred bytes.
const anatomy=JSON.parse(readFileSync(new URL('../public/models/atlas-male-complete.json',import.meta.url)));
const detailed={...state,selected:anatomy.parts.slice(0,3).map(part=>part.id),hidden:anatomy.parts.slice(0,3756).map(part=>part.id),visible:['cardiac','sensory','skeletal','muscular','arterial','venous','nervous','respiratory','digestive','urinary','lymphatic','endocrine','reproductive','connective'],view:'three-quarter',section:undefined,guestQuery:undefined};
const oldUrl=viewUrl(base,'male-detail',detailed,camera);
const compact=viewUrl(base,'male-detail',detailed,camera,anatomy),compactParams=viewParameters(compact);
assert.ok(oldUrl.length>25_000);assert.ok(compact.length<500);assert.ok(oldUrl.length/compact.length>100);
assert.equal(compactParams.has('hide'),false);assert.equal(compactParams.has('select'),false);
assert.equal(compactParams.has('view'),false);assert.equal(compactParams.has('region'),false);
const decoded=readViewUrl(compact,anatomy,{visible:[],view:'three-quarter',skinOpacity:.1});
assert.deepEqual(decoded.selected,detailed.selected);assert.deepEqual(decoded.hidden,detailed.hidden);
assert.deepEqual([...decoded.visible].sort(),[...detailed.visible].sort());assert.deepEqual(decoded.camera,camera);

// A scattered hidden set uses a bounded bitset instead of thousands of indices.
const scattered={...detailed,hidden:anatomy.parts.filter((_,index)=>index%2===0).map(part=>part.id)};
const scatteredUrl=viewUrl(base,'male-detail',scattered,camera,anatomy);
assert.ok(scatteredUrl.length<2000);assert.deepEqual(readViewUrl(scatteredUrl,anatomy,{}).hidden,scattered.hidden);

// Query overrides, embeds and old links still resolve; stale catalogues cannot select wrong parts.
const override=new URL(compact);override.searchParams.set('select',anatomy.parts[150].id);override.searchParams.set('layers','skeletal');
const overridden=readViewUrl(override.href,anatomy,{});assert.deepEqual(overridden.selected,[anatomy.parts[150].id]);assert.deepEqual(overridden.visible,['skeletal']);
assert.deepEqual(readViewUrl(embedUrl(compact,['study','systems']),anatomy,{}).hidden,detailed.hidden);
const changed={...anatomy,parts:[anatomy.parts[1],anatomy.parts[0],...anatomy.parts.slice(2)]};
assert.deepEqual(readViewUrl(compact,changed,{}).hidden,[]);assert.deepEqual(readViewUrl(compact,changed,{}).selected,[]);
const malformed=new URL(compact);malformed.searchParams.set('h','0-zzzz');assert.deepEqual(readViewUrl(malformed.href,anatomy,{}).hidden,[]);
console.log(`Large views, legacy links, overrides and packed anatomical IDs verified (${oldUrl.length} to ${compact.length} characters).`);
