import assert from 'node:assert/strict';
import {test} from 'node:test';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {InMemoryTransport} from '@modelcontextprotocol/sdk/inMemory.js';
import {anatomyView, catalogue, resolveAnatomy} from './atlas-data.mjs';
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
  assert.equal(male.parts.length, 5310);
  for (const [structure, system] of [
    ['Ureter (left)', 'urinary'],
    ['Femoral artery (right)', 'arterial'],
    ['Common iliac vein (left)', 'venous'],
    ['Sciatic nerve (right)', 'nervous'],
  ]) {
    assert.equal(male.parts.find(part => part.name === structure)?.system, system);
    assert.ok(resolveAnatomy(structure).length > 0);
  }
  assert.equal(male.parts.find(part => part.id === 'O3M:1st metacarpal bone.l')?.provenance.mirrored, true);
  assert.deepEqual(resolveAnatomy('ZA:Median nerve.r').map(concept => concept.id), ['ZA:Median nerve.r']);
  assert.ok(resolveAnatomy('jejunum').some(concept => concept.elements.length > 0));
  assert.ok(resolveAnatomy('ileum').some(concept => concept.elements.length > 0));
  assert.equal(resolveAnatomy('Ureters')[0].elements.length, 2);
  assert.ok(resolveAnatomy('small intestine').some(concept => concept.elements.length > 0));
  assert.ok(resolveAnatomy('Intestines').some(concept => concept.elements.length >= 40));
  assert.equal(resolveAnatomy('Ductus deferens')[0].elements.length, 2);
  assert.equal(resolveAnatomy('Blood vessels of the penis')[0].elements.length, 6);
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

test('MCP discovery, tool call, and UI resource work together', async () => {
  const [clientTransport, serverTransport] = InMemoryTransport.createLinkedPair();
  const server = createServer();
  const client = new Client({name: 'human-atlas-test', version: '1.0.0'});
  await Promise.all([server.connect(serverTransport), client.connect(clientTransport)]);
  try {
    const tools = (await client.listTools()).tools;
    assert.deepEqual(tools.map(tool => tool.name), ['search_anatomy', 'show_anatomy']);
    assert.equal(tools.find(tool => tool.name === 'show_anatomy')._meta.ui.resourceUri, 'ui://human-atlas/anatomy-view');
    const found = await client.callTool({name: 'search_anatomy', arguments: {query: 'stomach'}});
    assert.equal(found.structuredContent.matches[0].id, 'ZA:Stomach');
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
