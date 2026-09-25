import assert from 'node:assert/strict';
import {test} from 'node:test';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {InMemoryTransport} from '@modelcontextprotocol/sdk/inMemory.js';
import {build} from 'esbuild';
import {anatomyView, catalogue, DEPTH_LAYERS, resolveAnatomy} from './atlas-data.mjs';
import {createServer} from './server.mjs';

test('the requested structure becomes a selected, focused embed', () => {
  const result = anatomyView({structure: 'Stomach'});
  const url = new URL(result.url);
  assert.equal(url.searchParams.get('model'), 'male-detail');
  assert.deepEqual(url.searchParams.getAll('select'), ['ZA:Stomach']);
  assert.equal(url.searchParams.get('focus'), '1');
  assert.equal(url.searchParams.get('embed'), '1');
  assert.ok(url.searchParams.get('layers').includes('digestive'));
  assert.ok(result.iframe.includes('&amp;select=ZA%3AStomach'));
  assert.deepEqual(result.systems, ['digestive']);
});

test('local models resolve to the local viewer and local geometry', () => {
  const result = anatomyView({structure: 'Stomach', model: 'local-female'});
  const url = new URL(result.url);
  assert.equal(url.origin, 'http://localhost:3016');
  assert.equal(url.searchParams.get('model'), 'local-female');
  assert.equal(url.searchParams.get('skin'), '0');
  assert.equal(result.structures[0].pieces, 2);
});

test('served catalogues select publisher-source male and female assets', () => {
  assert.equal(catalogue('male-detail').source, 'Z-Anatomy + Open 3D Model + BodyParts3D 4.0');
  assert.equal(catalogue('male-full').source, 'Z-Anatomy + Open 3D Model + BodyParts3D 4.0');
  const male = catalogue('male-detail');
  assert.equal(male.parts.length, 5209);
  for (const [structure, system] of [
    ['Ureter (left)', 'urinary'],
    ['Femoral artery (right)', 'arterial'],
    ['Common iliac vein (left)', 'venous'],
    ['Sciatic nerve (right)', 'nervous'],
  ]) {
    assert.equal(male.parts.find(part => part.name === structure)?.system, system);
    assert.ok(resolveAnatomy(structure).length > 0);
  }
  assert.equal(male.parts.find(part => part.id === 'O3M:Sesamoid bones of hand.l')?.provenance.mirrored, true);
  assert.equal(male.parts.some(part => part.id.startsWith('O3M:') && part.system === 'skeletal' && !/Sesamoid bones of hand/i.test(part.name)), false);
  assert.equal(male.parts.some(part => part.id.startsWith('O3M:') && /Brachiocephalic artery|Arm superficial vein-(Basilic|Cephalic|Median antebrachial|Median cubital) vein/i.test(part.name)), false);
  assert.deepEqual(resolveAnatomy('ZA:Median nerve.r').map(concept => concept.id), ['ZA:Median nerve.r']);
  assert.ok(resolveAnatomy('jejunum').some(concept => concept.elements.length > 0));
  assert.equal(male.parts.some(part => part.id.startsWith('FJ') && /ileum|ileocecal junction|mesentery of small intestine/i.test(part.name)), false);
  assert.equal(male.parts.some(part => part.id.startsWith('FJ') && /pancreatic duct tree|deferent duct|superficial dorsal vein of penis/i.test(part.name)), false);
  assert.equal(resolveAnatomy('Ureters')[0].elements.length, 2);
  assert.ok(resolveAnatomy('small intestine').some(concept => concept.elements.length > 0));
  assert.ok(resolveAnatomy('Intestines').some(concept => concept.elements.length >= 8));
  assert.equal(resolveAnatomy('Ductus deferens')[0].elements.length, 2);
  assert.equal(resolveAnatomy('Blood vessels of the penis')[0].elements.length, 4);
  assert.equal(catalogue('female').source, 'Human Reference Atlas');
  assert.deepEqual(anatomyView({structure: 'Uterus', model: 'female'}).systems, ['reproductive']);
  assert.equal(catalogue('female').parts.some(part=>part.name==='Amnion'),false);
  assert.equal(catalogue('embryo').parts.length,8);
  assert.equal(resolveAnatomy('Amnion', 'embryo').length,1);
  assert.deepEqual(anatomyView({structure: 'Amnion', model: 'embryo'}).systems,['pregnancy']);
});

test('cell model exposes its original components through the viewer and MCP', () => {
  const cell = catalogue('cell');
  assert.equal(cell.scope, 'cell');
  assert.equal(cell.parts.length, 20);
  assert.equal(cell.parts.every(part => part.colors !== undefined), true);
  assert.equal(cell.parts.reduce((sum, part) => sum + part.indexCount / 3, 0), cell.triangles);
  const view = anatomyView({structure: 'Nucleolus', model: 'cell'});
  const url = new URL(view.url);
  assert.equal(url.searchParams.get('select'), 'CELL:2');
  assert.equal(url.searchParams.get('skin'), '0.18');
  assert.ok(url.searchParams.get('layers').includes('cell-nucleus'));
  assert.ok(view.iframe.includes('model=cell'));
  assert.deepEqual(resolveAnatomy('Mitochondria', 'cell')[0].elements, ['CELL:10', 'CELL:18']);
});

test('embed controls can hide the explode dock while keeping systems available', () => {
  const result = anatomyView({structure: 'Stomach', controls: ['systems', 'open']});
  const url = new URL(result.url);
  assert.equal(url.searchParams.get('ui'), 'systems,open');
  assert.ok(result.iframe.includes('ui=systems%2Copen'));
  assert.equal(new URL(anatomyView({structure: 'Stomach', controls: []}).url).searchParams.get('ui'), '');
  assert.throws(() => anatomyView({structure: 'Stomach', controls: ['unknown']}), /Unknown embed control/);
});

test('exact IDs choose the requested side and bilateral names choose both', () => {
  const id = 'ZA:Sternocostal head of pectoralis major muscle.l';
  assert.deepEqual(resolveAnatomy(id).map(c => c.id), [id]);
  const left = anatomyView({structure: id});
  assert.deepEqual(new URL(left.url).searchParams.getAll('select'), [id]);
  assert.deepEqual(resolveAnatomy('left sternocostal head of pectoralis major muscle').map(c => c.id), [id]);
  assert.deepEqual(resolveAnatomy('Sternocostal head of pectoralis major muscle').map(c => c.id), [
    'ZA:Sternocostal head of pectoralis major muscle.l',
    'ZA:Sternocostal head of pectoralis major muscle.r',
  ]);
  assert.throws(() => resolveAnatomy('not a modeled structure'), /No modeled structure matches/);
});

test('MCP views cover the shareable viewer settings and whole models', () => {
  const whole=new URL(anatomyView({model:'female',hierarchy:'depth',depthHidden:['skin'],explode:.5,region:'lower-left',view:'back',skinOpacity:.35,labels:false,section:{axis:'sagittal',position:.4,flip:true},camera:[1,2,3,0,1,0],controls:['systems','explode','study']}).url);
  assert.deepEqual(whole.searchParams.getAll('select'),[]);
  for(const [key,value] of Object.entries({tree:'depth',depth:'skin',explode:'0.5',region:'lower-left',view:'back',skin:'0.35',labels:'0',cut:'sagittal',slice:'0.4',flip:'1',camera:'1,2,3,0,1,0',ui:'systems,explode,study'}))assert.equal(whole.searchParams.get(key),value,key);
  const oblique=new URL(anatomyView({section:{axis:'oblique',position:.31,flip:false,azimuth:-42,elevation:27},controls:['sections']}).url);
  for(const [key,value] of Object.entries({cut:'oblique',slice:'0.31',azimuth:'-42',elevation:'27',ui:'sections'}))assert.equal(oblique.searchParams.get(key),value,key);
  const paired=new URL(anatomyView({sections:[{axis:'axial',position:.25,flip:true},{axis:'oblique',position:.62,flip:false,azimuth:24,elevation:-12}],activeSection:1}).url);
  for(const [key,value] of Object.entries({sections:'2',sectionTab:'2',cut:'axial',slice:'0.25',cut2:'oblique',slice2:'0.62',azimuth2:'24',elevation2:'-12'}))assert.equal(paired.searchParams.get(key),value,key);
  assert.equal(new URL(anatomyView({rotate:true}).url).searchParams.get('rotate'),'1');
  assert.throws(()=>anatomyView({rotate:true,explode:.5}),/Auto rotation/);
  const selected=anatomyView({structures:['Stomach','Pancreas'],hidden:['ZA:Liver'],systems:['digestive'],isolate:true,context:.2,focus:false});
  const url=new URL(selected.url);
  assert.deepEqual(url.searchParams.getAll('select'),['ZA:Stomach','ZA:Pancreas']);
  assert.ok(url.searchParams.getAll('hide').length>0);
  assert.equal(url.searchParams.get('layers'),'digestive');
  assert.equal(url.searchParams.get('isolate'),'1');
  assert.equal(url.searchParams.has('focus'),false);
  assert.equal(selected.selectedPieces.length,2);
  assert.equal(new URL(anatomyView({model:'embryo',hierarchy:'depth',depthHidden:['skin']}).url).searchParams.get('depth'),'skin');
  assert.equal(new URL(anatomyView({hierarchy:'guest:chakras'}).url).searchParams.get('tree'),'guest:chakras');
  const guestView=new URL(anatomyView({hierarchy:'guest:custom',guestUrls:['https://example.com/custom.json']}).url);
  assert.equal(guestView.searchParams.get('tree'),'guest:custom');
  assert.deepEqual(guestView.searchParams.getAll('guest'),['https://example.com/custom.json']);
  assert.throws(()=>anatomyView({hierarchy:'guest:custom'}),/guest URL/);
  assert.throws(()=>anatomyView({model:'cell',hierarchy:'depth'}),/hierarchy/);
  assert.throws(()=>anatomyView({depthHidden:['missing']}),/depth layer/);
  assert.throws(()=>anatomyView({systems:['missing']}),/system/);
  assert.throws(()=>anatomyView({isolate:true}),/Isolation/);
});

test('browser WebMCP exposes live view and saved view operations',async()=>{
  const bundle=await build({entryPoints:['app/agent-tools.ts','app/depth-layers.ts'],bundle:true,platform:'node',format:'esm',write:false,outdir:'/tmp/agent-tools-test'});
  const load=async suffix=>import(`data:text/javascript;base64,${Buffer.from(bundle.outputFiles.find(file=>file.path.endsWith(suffix)).text).toString('base64')}`);
  const {atlasTools}=await load('agent-tools.js'),{DEPTH_LAYERS:browserLayers}=await load('depth-layers.js');
  assert.deepEqual(DEPTH_LAYERS,browserLayers.map(layer=>layer.id));
  const calls=[];
  const actions={snapshot:()=>({model:'male-detail',hierarchy:'depth'}),update:change=>{calls.push(change);return change;},command:action=>{calls.push(action);return {action};},url:()=> 'http://localhost:3016/'};
  const tools=atlasTools(catalogue('male-detail'),()=>{},actions),byName=name=>tools.find(tool=>tool.name===name);
  assert.deepEqual(tools.map(tool=>tool.name),['find_anatomy','inspect_anatomical_structure','get_anatomy_view','set_anatomy_view','act_on_anatomy_view','list_saved_anatomy_views','save_anatomy_view','open_saved_anatomy_view','delete_saved_anatomy_view']);
  assert.deepEqual(byName('get_anatomy_view').execute({}),{model:'male-detail',hierarchy:'depth'});
  byName('set_anatomy_view').execute({hierarchy:'depth',structures:['Stomach'],depthHidden:['skin'],explode:.5,section:{axis:'coronal',position:.3,flip:true}});
  assert.deepEqual(calls[0],{hierarchy:'depth',structures:['Stomach'],depthHidden:['skin'],explode:.5,section:{axis:'coronal',position:.3,flip:true}});
  byName('set_anatomy_view').execute({section:{axis:'oblique',position:.6,flip:false,azimuth:112,elevation:-16}});
  assert.deepEqual(calls[1],{section:{axis:'oblique',position:.6,flip:false,azimuth:112,elevation:-16}});
  byName('set_anatomy_view').execute({sections:[{axis:'axial',position:.25,flip:true},null],activeSection:1});
  assert.deepEqual(calls[2],{sections:[{axis:'axial',position:.25,flip:true},null],activeSection:1});
  byName('set_anatomy_view').execute({hierarchy:'guest:chakras'});
  assert.deepEqual(calls[3],{hierarchy:'guest:chakras'});
  assert.throws(()=>byName('set_anatomy_view').execute({depthHidden:['missing']}),/depth layer/);
  assert.throws(()=>byName('set_anatomy_view').execute({model:'female',explode:.2}),/separate call/);
  byName('act_on_anatomy_view').execute({action:'reset'});
  assert.equal(calls[4],'reset');
});

test('MCP discovery, tool call, and UI resource work together', async () => {
  const [clientTransport, serverTransport] = InMemoryTransport.createLinkedPair();
  const server = createServer();
  const client = new Client({name: 'human-atlas-test', version: '1.0.0'});
  await Promise.all([server.connect(serverTransport), client.connect(clientTransport)]);
  try {
    const tools = (await client.listTools()).tools;
    assert.deepEqual(tools.map(tool => tool.name), ['get_anatomy_options','search_anatomy', 'show_anatomy']);
    assert.equal(tools.find(tool => tool.name === 'show_anatomy')._meta.ui.resourceUri, 'ui://human-atlas/anatomy-view');
    const found = await client.callTool({name: 'search_anatomy', arguments: {query: 'stomach'}});
    assert.equal(found.structuredContent.matches[0].id, 'ZA:Stomach');
    const options = await client.callTool({name: 'get_anatomy_options', arguments: {model: 'female'}});
    assert.ok(options.structuredContent.depthLayers.includes('superficial-muscles'));
    const shown = await client.callTool({name: 'show_anatomy', arguments: {structure: 'Stomach'}});
    assert.equal(shown.isError, undefined);
    assert.equal(new URL(shown.structuredContent.url).searchParams.get('select'), 'ZA:Stomach');
    const customized = await client.callTool({name: 'show_anatomy', arguments: {structure: 'Stomach', controls: ['systems', 'open']}});
    assert.equal(new URL(customized.structuredContent.url).searchParams.get('ui'), 'systems,open');
    const missing = await client.callTool({name: 'show_anatomy', arguments: {structure: 'not a modeled structure'}});
    assert.equal(missing.isError, true);
    const ui = (await client.readResource({uri: 'ui://human-atlas/anatomy-view'})).contents[0];
    assert.equal(ui.mimeType, 'text/html;profile=mcp-app');
    assert.deepEqual(ui._meta.ui.csp.frameDomains, ['https://ctzurcanu.github.io','http://localhost:3016']);
    assert.ok(ui.text.includes('ui/notifications/tool-result'));
  } finally {
    await client.close();
    await server.close();
  }
});
