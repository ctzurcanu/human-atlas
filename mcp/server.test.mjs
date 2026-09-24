import assert from 'node:assert/strict';
import {test} from 'node:test';
import {Client} from '@modelcontextprotocol/sdk/client/index.js';
import {InMemoryTransport} from '@modelcontextprotocol/sdk/inMemory.js';
import {anatomyView, resolveAnatomy} from './atlas-data.mjs';
import {createServer} from './server.mjs';

test('the requested structure becomes a selected, focused embed', () => {
  const result = anatomyView({structure: 'Stomach'});
  const url = new URL(result.url);
  assert.equal(url.searchParams.get('model'), 'male-detail');
  assert.deepEqual(url.searchParams.getAll('select'), ['DETAIL:Stomach']);
  assert.equal(url.searchParams.get('focus'), '1');
  assert.equal(url.searchParams.get('embed'), '1');
  assert.ok(url.searchParams.get('layers').includes('digestive'));
  assert.ok(result.iframe.includes('&amp;select=DETAIL%3AStomach'));
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

test('embed controls can hide the explode dock while keeping systems available', () => {
  const result = anatomyView({structure: 'Stomach', controls: ['systems', 'open']});
  const url = new URL(result.url);
  assert.equal(url.searchParams.get('ui'), 'systems,open');
  assert.ok(result.iframe.includes('ui=systems%2Copen'));
  assert.equal(new URL(anatomyView({structure: 'Stomach', controls: []}).url).searchParams.get('ui'), '');
  assert.throws(() => anatomyView({structure: 'Stomach', controls: ['unknown']}), /Unknown embed control/);
});

test('exact IDs choose the requested side, bilateral names choose both, and duplicate names need an ID', () => {
  const id = 'DETAIL:Sternocostal head of pectoralis major muscle.l';
  assert.deepEqual(resolveAnatomy(id).map(c => c.id), [id]);
  const left = anatomyView({structure: id});
  assert.deepEqual(new URL(left.url).searchParams.getAll('select'), [id]);
  assert.deepEqual(resolveAnatomy('left sternocostal head of pectoralis major muscle').map(c => c.id), [id]);
  assert.deepEqual(resolveAnatomy('Sternocostal head of pectoralis major muscle').map(c => c.id), [
    'DETAIL:Sternocostal head of pectoralis major muscle.l',
    'DETAIL:Sternocostal head of pectoralis major muscle.r',
  ]);
  assert.throws(() => resolveAnatomy('Sternocostal head of pectoralis major muscle (insertion 1)'), /Several structures have that name/);
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
    assert.equal(found.structuredContent.matches[0].id, 'DETAIL:Stomach');
    const shown = await client.callTool({name: 'show_anatomy', arguments: {structure: 'Stomach'}});
    assert.equal(shown.isError, undefined);
    assert.equal(new URL(shown.structuredContent.url).searchParams.get('select'), 'DETAIL:Stomach');
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
