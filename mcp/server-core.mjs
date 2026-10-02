#!/usr/bin/env node
import {McpServer} from '@modelcontextprotocol/sdk/server/mcp.js';
import {z} from 'zod';
import {anatomyOptions, anatomyView, DEPTH_LAYERS, EMBED_CONTROLS, MODELS, REGIONS, searchAnatomy, SECTION_AXES, VIEWS} from './atlas-core.mjs';

const UI_URI = 'ui://human-atlas/anatomy-view';
const UI_MIME = 'text/html;profile=mcp-app';

export function createServer({uiHtml,publicOnly=false}) {
  const modelSchema = z.enum(Object.keys(MODELS).filter(id=>!publicOnly||!id.startsWith('local-')));
  const structureSchema = z.object({
    id: z.string().describe('Exact atlas concept identifier.'),
    name: z.string().describe('Anatomical structure name.'),
    pieces: z.number().int().nonnegative().describe('Number of modeled mesh pieces.'),
    terminology: z.record(z.unknown()).nullable().describe('Available anatomical terminology metadata, or null.'),
  });
  const uiMeta={ui:{csp:{frameDomains:publicOnly?['https://ctzurcanu.github.io']:['https://ctzurcanu.github.io','http://localhost:3016']},prefersBorder:false}};
  const server = new McpServer({name: 'human-atlas', version: '0.1.0'});

  server.registerTool('get_anatomy_options', {
    title: 'Get Atlas view options',
    description: 'List model, system, hierarchy, region, depth layer, camera view and embed control IDs accepted by show_anatomy.',
    inputSchema: {model: modelSchema.default('male-detail')},
    outputSchema: {
      model: modelSchema,
      models: z.array(modelSchema),
      systems: z.array(z.object({id:z.string(),pieces:z.number().int().nonnegative()})),
      hierarchies: z.array(z.string()),
      regions: z.array(z.enum(REGIONS)),
      depthLayers: z.array(z.enum(DEPTH_LAYERS)),
      views: z.array(z.enum(VIEWS)),
      controls: z.array(z.enum(EMBED_CONTROLS)),
    },
    annotations: {title:'Get Atlas view options',readOnlyHint: true,destructiveHint:false,openWorldHint:false},
  }, async ({model}) => {
    const options = anatomyOptions(model);
    if(publicOnly)options.models=options.models.filter(id=>!id.startsWith('local-'));
    return {content: [{type: 'text', text: JSON.stringify(options)}], structuredContent: options};
  });

  server.registerTool('search_anatomy', {
    title: 'Search Atlas anatomy',
    description: 'Find exact names and IDs of structures modeled in Atlas. Search before showing a view when a request is ambiguous.',
    inputSchema: {
      query: z.string().min(1).describe('Anatomical name or atlas ID, such as stomach or sternocostal head.'),
      model: modelSchema.default('male-detail').describe('Atlas model. The detailed male model is the default.'),
    },
    outputSchema: {model:modelSchema,matches:z.array(structureSchema).describe('Matching structures, ordered by relevance; empty when none match.')},
    annotations: {title:'Search Atlas anatomy',readOnlyHint: true,destructiveHint:false,openWorldHint:false},
  }, async ({query, model}) => {
    const matches = searchAnatomy(query, model);
    return {content: [{type: 'text', text: JSON.stringify({model, matches})}], structuredContent: {model, matches}};
  });

  server.registerTool('show_anatomy', {
    title: 'Show Atlas 3D anatomy',
    description: 'Create an interactive Atlas view and copyable iframe. Optionally select one or more exact structure names or IDs; search_anatomy resolves ambiguity. View options control Systems, Regions, Depth, guest hierarchies, visibility, Explode, cuts, camera and isolation. Omit structure for the whole model.'+(publicOnly?'':' Models local-reference, local-male and local-female require the local development viewer on port 3016.'),
    inputSchema: {
      structure: z.string().min(1).optional().describe('Exact structure name, atlas concept ID, or piece ID.'),
      structures: z.array(z.string().min(1)).optional().describe('Additional structures to select together.'),
      model: modelSchema.default('male-detail'),
      view: z.enum(VIEWS).default('three-quarter').describe('Camera direction. Defaults to three-quarter.'),
      context: z.number().min(0).max(1).default(0.18).describe('Opacity of surrounding anatomy, 0 to 1.'),
      hierarchy: z.string().default('systems').describe('Active hierarchy: systems, regions, depth, guest:chakras, or guest:<id> for a supplied guest URL. Depth is unavailable for the cell model.'),
      guestUrls: z.array(z.string().url()).optional().describe('Public JSON hierarchy URLs to add to the Additional hierarchies menu. Custom guest:<id> needs its JSON URL here.'),
      systems: z.array(z.string()).optional().describe('Visible system IDs. Omit for the model defaults; [] hides unselected anatomy.'),
      depthHidden: z.array(z.enum(DEPTH_LAYERS)).optional().describe('Depth layer IDs to hide.'),
      hidden: z.array(z.string()).optional().describe('Exact names, concept IDs, or piece IDs to hide.'),
      region: z.enum(REGIONS).default('all'),
      explode: z.number().min(0).max(1).default(0).describe('Hierarchical explosion amount, 0 assembled to 1 individual pieces.'),
      skinOpacity: z.number().min(0).max(1).optional(),
      labels: z.boolean().default(true),
      isolate: z.boolean().default(false),
      rotate: z.boolean().default(false).describe('Auto rotate the view.'),
      section: z.object({axis: z.enum(SECTION_AXES), position: z.number().min(0).max(1), flip: z.boolean().default(false), azimuth: z.number().min(-180).max(180).optional(), elevation: z.number().min(-90).max(90).optional()}).optional().describe('Enable a geometric cross-section. Oblique uses azimuth and elevation in degrees.'),
      sections: z.array(z.object({axis: z.enum(SECTION_AXES), position: z.number().min(0).max(1), flip: z.boolean().default(false), azimuth: z.number().min(-180).max(180).optional(), elevation: z.number().min(-90).max(90).optional()}).nullable()).min(1).max(2).optional().describe('One or two simultaneous geometric cuts. Null closes a cut while keeping its tab.'),
      activeSection: z.number().int().min(0).max(1).optional().describe('Active section tab, zero based.'),
      camera: z.array(z.number()).optional().describe('Camera position and target as six numbers, optionally followed by two projection offsets.'),
      focus: z.boolean().optional().describe('Fit the selected anatomy in view.'),
      controls: z.array(z.enum(EMBED_CONTROLS)).optional().describe('Controls visible in the iframe: model, search, sections, study, systems, camera, explode, details, open, download. Omit for defaults; [] shows only the 3D view.'),
    },
    _meta: {ui: {resourceUri: UI_URI}},
    outputSchema: {
      model: modelSchema,
      structures: z.array(structureSchema),
      systems: z.array(z.string()).describe('Systems belonging to selected anatomy.'),
      selectedPieces: z.array(z.string()),
      hiddenPieces: z.array(z.string()),
      url: z.string().url().describe('Interactive anatomy viewer URL.'),
      iframe: z.string().describe('HTML iframe for embedding this view.'),
    },
    annotations: {title:'Show Atlas 3D anatomy',readOnlyHint: true,destructiveHint:false,openWorldHint:true},
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

  server.registerResource('Atlas interactive anatomy view', UI_URI, {
    description: 'Interactive Atlas view for show_anatomy results.',
    mimeType: UI_MIME,
    _meta: uiMeta,
  }, async () => ({contents: [{
    uri: UI_URI,
    mimeType: UI_MIME,
    text: uiHtml,
    _meta: uiMeta,
  }]}));

  return server;
}
