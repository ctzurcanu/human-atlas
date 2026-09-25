#!/usr/bin/env node
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {StdioServerTransport} from '@modelcontextprotocol/sdk/server/stdio.js';
import {z} from 'zod';
import {anatomyOptions, anatomyView, DEPTH_LAYERS, EMBED_CONTROLS, HIERARCHIES, MODELS, REGIONS, searchAnatomy, SECTION_AXES, VIEWS} from './atlas-data.mjs';

const UI_URI = 'ui://human-atlas/anatomy-view';
const UI_MIME = 'text/html;profile=mcp-app';
const UI_HTML = readFileSync(new URL('./view.html', import.meta.url), 'utf8');
const modelSchema = z.enum(Object.keys(MODELS));

export function createServer() {
  const server = new McpServer({name: 'human-atlas', version: '0.1.0'});

  server.registerTool('get_anatomy_options', {
    title: 'Get Human Atlas view options',
    description: 'List model, system, hierarchy, region, depth layer, camera view and embed control IDs accepted by show_anatomy.',
    inputSchema: {model: modelSchema.default('male-detail')},
    annotations: {readOnlyHint: true},
  }, async ({model}) => {
    const options = anatomyOptions(model);
    return {content: [{type: 'text', text: JSON.stringify(options)}], structuredContent: options};
  });

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
    description: 'Create an interactive Human Atlas view and copyable iframe. Optionally select one or more exact structure names or IDs; search_anatomy resolves ambiguity. View options control Systems, Regions, Depth, visibility, Explode, cuts, camera and isolation. Omit structure for the whole model. Models local-male and local-female require the local development viewer on port 3016.',
    inputSchema: {
      structure: z.string().min(1).optional().describe('Exact structure name, atlas concept ID, or piece ID.'),
      structures: z.array(z.string().min(1)).optional().describe('Additional structures to select together.'),
      model: modelSchema.default('male-detail'),
      view: z.enum(VIEWS).default('three-quarter').describe('Camera direction. Defaults to three-quarter.'),
      context: z.number().min(0).max(1).default(0.18).describe('Opacity of surrounding anatomy, 0 to 1.'),
      hierarchy: z.enum(HIERARCHIES).default('systems').describe('Active hierarchy tab. Depth is unavailable for the cell model.'),
      systems: z.array(z.string()).optional().describe('Visible system IDs. Omit for the model defaults; [] hides unselected anatomy.'),
      depthHidden: z.array(z.enum(DEPTH_LAYERS)).optional().describe('Depth layer IDs to hide.'),
      hidden: z.array(z.string()).optional().describe('Exact names, concept IDs, or piece IDs to hide.'),
      region: z.enum(REGIONS).default('all'),
      explode: z.number().min(0).max(1).default(0).describe('Hierarchical explosion amount, 0 assembled to 1 individual pieces.'),
      skinOpacity: z.number().min(0).max(1).optional(),
      labels: z.boolean().default(true),
      isolate: z.boolean().default(false),
      rotate: z.boolean().default(false).describe('Auto rotate the view.'),
      section: z.object({axis: z.enum(SECTION_AXES), position: z.number().min(0).max(1), flip: z.boolean().default(false)}).optional().describe('Enable a geometric cross-section.'),
      camera: z.array(z.number()).optional().describe('Camera position and target as six numbers, optionally followed by two projection offsets.'),
      focus: z.boolean().optional().describe('Fit the selected anatomy in view.'),
      controls: z.array(z.enum(EMBED_CONTROLS)).optional().describe('Controls visible in the iframe: model, search, study, systems, camera, explode, details, open, download. Omit for defaults; [] shows only the 3D view.'),
    },
    _meta: {ui: {resourceUri: UI_URI}},
    annotations: {readOnlyHint: true},
  }, async args => {
    try {
      const result = anatomyView(args);
      return {
        content: [{type: 'text', text: `${result.structures.map(s => s.name).join(' and ') || 'Whole model'} (${result.model})\nView: ${result.url}\nEmbed HTML: ${result.iframe}`}],
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
