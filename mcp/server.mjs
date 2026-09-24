#!/usr/bin/env node
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {StdioServerTransport} from '@modelcontextprotocol/sdk/server/stdio.js';
import {z} from 'zod';
import {anatomyView, EMBED_CONTROLS, MODELS, searchAnatomy, VIEWS} from './atlas-data.mjs';

const UI_URI = 'ui://human-atlas/anatomy-view';
const UI_MIME = 'text/html;profile=mcp-app';
const UI_HTML = readFileSync(new URL('./view.html', import.meta.url), 'utf8');
const modelSchema = z.enum(Object.keys(MODELS));

export function createServer() {
  const server = new McpServer({name: 'human-atlas', version: '0.1.0'});

  server.registerTool('search_anatomy', {
    title: 'Search Human Atlas anatomy',
    description: 'Find exact names and IDs of structures modeled in Human Atlas. Search before showing a view when a request is ambiguous.',
    inputSchema: {
      query: z.string().min(1).describe('Anatomical name or atlas ID, such as stomach or sternocostal head.'),
      model: modelSchema.default('male-detail').describe('Atlas model. The detailed male model is the default.'),
    },
    annotations: {readOnlyHint: true},
  }, async ({query, model}) => {
    const matches = searchAnatomy(query, model);
    return {content: [{type: 'text', text: JSON.stringify({model, matches})}], structuredContent: {model, matches}};
  });

  server.registerTool('show_anatomy', {
    title: 'Show Human Atlas 3D anatomy',
    description: 'Create an interactive Human Atlas view and copyable iframe for a modeled structure. Use an exact structure name or ID; search_anatomy resolves ambiguity. For a request like "show the male stomach", pass structure="Stomach" and model="male-detail". Pass controls to choose which iframe controls appear; by default Study, Camera controls, Explode, and PNG download are hidden. Pass [] for a bare viewer. Models local-male and local-female use the local development viewer on port 3016 and require npm run dev.',
    inputSchema: {
      structure: z.string().min(1).describe('Exact structure name or atlas concept ID.'),
      model: modelSchema.default('male-detail'),
      view: z.enum(VIEWS).default('three-quarter').describe('Camera direction. Defaults to three-quarter.'),
      context: z.number().min(0).max(1).default(0.18).describe('Opacity of surrounding anatomy, 0 to 1.'),
      controls: z.array(z.enum(EMBED_CONTROLS)).optional().describe('Controls visible in the iframe: model, search, study, systems, camera, explode, details, open, download. Omit for defaults; [] shows only the 3D view.'),
    },
    _meta: {ui: {resourceUri: UI_URI}},
    annotations: {readOnlyHint: true},
  }, async args => {
    try {
      const result = anatomyView(args);
      return {
        content: [{type: 'text', text: `${result.structures.map(s => s.name).join(' and ')} (${result.model})\nView: ${result.url}\nEmbed HTML: ${result.iframe}`}],
        structuredContent: result,
      };
    } catch (error) {
      return {isError: true, content: [{type: 'text', text: error instanceof Error ? error.message : String(error)}]};
    }
  });

  server.registerResource('Human Atlas interactive anatomy view', UI_URI, {
    description: 'Interactive Human Atlas view for show_anatomy results.',
    mimeType: UI_MIME,
    _meta: {ui: {csp: {frameDomains: ['https://ctzurcanu.github.io','http://localhost:3016']}, prefersBorder: false}},
  }, async () => ({contents: [{
    uri: UI_URI,
    mimeType: UI_MIME,
    text: UI_HTML,
    _meta: {ui: {csp: {frameDomains: ['https://ctzurcanu.github.io','http://localhost:3016']}, prefersBorder: false}},
  }]}));

  return server;
}

if (process.argv[1] && import.meta.url === pathToFileURL(process.argv[1]).href) {
  const server = createServer();
  await server.connect(new StdioServerTransport());
}
