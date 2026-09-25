import {readFileSync} from 'node:fs';

export const VIEWER_URL = 'https://ctzurcanu.github.io/human-atlas/';
export const MODELS = {
  'male-detail': 'atlas-male-complete.json',
  'male-full': 'atlas-male-complete.json',
  male: 'atlas.json',
  female: 'atlas-hra-female.json',
  embryo: 'atlas-embryo.json',
  cell: 'atlas-cell.json',
  'local-male': 'male.json',
  'local-female': 'female.json',
};
export const VIEWS = ['three-quarter', 'front', 'back', 'side', 'right', 'superior', 'inferior'];
export const HIERARCHIES = ['systems', 'regions', 'depth'];
export const REGIONS = ['all', 'head-neck', 'torso', 'upper-right', 'upper-left', 'lower-right', 'lower-left'];
export const DEPTH_LAYERS = ['skin', 'superficial-veins', 'investing-fascia', 'superficial-muscles', 'second-muscles', 'intermediate-muscles', 'deep-muscles', 'deepest-muscles', 'visceral-coverings', 'anterior-organs', 'deep-organs', 'deep-vessels', 'lymphatic', 'deep-nerves', 'ligaments', 'other', 'thoracic-bones', 'limb-bones', 'other-bones', 'spine', 'skull'];
export const SECTION_AXES = ['axial', 'sagittal', 'coronal'];
export const EMBED_CONTROLS = ['model', 'search', 'study', 'systems', 'camera', 'explode', 'details', 'open', 'download'];
export const DEFAULT_LAYERS = [
  'cardiac', 'sensory', 'skeletal', 'muscular', 'arterial', 'venous', 'nervous',
  'respiratory', 'digestive', 'urinary', 'lymphatic', 'endocrine', 'reproductive', 'connective',
];

const catalogues = new Map();
export function catalogue(model = 'male-detail') {
  const filename = MODELS[model];
  if (!filename) throw new Error(`Unknown model: ${model}`);
  if (!catalogues.has(model)) {
    const url = new URL(model.startsWith('local-')?`../.local-models/${filename}`:`../public/models/${filename}`, import.meta.url);
    catalogues.set(model, JSON.parse(readFileSync(url, 'utf8')));
  }
  return catalogues.get(model);
}
export function anatomyOptions(model = 'male-detail') {
  const atlas = catalogue(model);
  const counts = new Map();
  for (const part of atlas.parts) counts.set(part.system, (counts.get(part.system) ?? 0) + 1);
  return {model, models:Object.keys(MODELS), systems:[...counts].map(([id,pieces])=>({id,pieces})), hierarchies:model==='cell'?HIERARCHIES.filter(id=>id!=='depth'):HIERARCHIES, regions:model==='cell'?['all']:REGIONS, depthLayers:model==='cell'?[]:DEPTH_LAYERS, views:VIEWS, controls:EMBED_CONTROLS};
}

const normalize = value => value.trim().toLocaleLowerCase().replace(/\s+/g, ' ');
const plainName = value => normalize(value.replace(/\s*\((?:left|right)\)$/i, ''));
const termsFor = value => {
  const term = normalize(value);
  const side = /^(left|right)\s+(.+)$/.exec(term);
  return side ? [term, `${side[2]} (${side[1]})`] : [term];
};

export function searchAnatomy(query, model = 'male-detail', limit = 20) {
  const terms = termsFor(query);
  const term = terms.at(-1);
  if (!term) throw new Error('Enter an anatomical structure or atlas ID.');
  const atlas = catalogue(model);
  return atlas.concepts
    .filter(concept => normalize(concept.name).includes(term) || normalize(concept.id).includes(term))
    .sort((a, b) => {
      const score = concept => normalize(concept.id) === term || normalize(concept.name) === term ? 0
        : plainName(concept.name) === term ? 1
        : normalize(concept.name).startsWith(term) ? 2 : 3;
      return score(a) - score(b) || a.name.length - b.name.length || a.name.localeCompare(b.name) || a.id.localeCompare(b.id);
    })
    .slice(0, limit)
    .map(concept => ({id: concept.id, name: concept.name, pieces: concept.elements.length}));
}

export function resolveAnatomy(structure, model = 'male-detail') {
  const term = termsFor(structure).at(-1);
  if (!term) throw new Error('Enter an anatomical structure or atlas ID.');
  const atlas = catalogue(model);
  const byId = atlas.concepts.find(concept => normalize(concept.id) === term);
  if (byId) return [byId];
  const exact = atlas.concepts.filter(concept => normalize(concept.name) === term);
  if (exact.length === 1) return exact;
  if (exact.length > 1) throw new Error(`Several structures have that name. Use an exact atlas ID: ${exact.slice(0, 8).map(c => c.id).join(', ')}`);
  const bilateral = atlas.concepts.filter(concept => plainName(concept.name) === term);
  if (bilateral.length === 2 && bilateral.some(c => /\(left\)$/i.test(c.name)) && bilateral.some(c => /\(right\)$/i.test(c.name))) return bilateral;
  const suggestions = searchAnatomy(structure, model, 8);
  throw new Error(suggestions.length
    ? `No exact structure match for "${structure}". Search results: ${suggestions.map(c => `${c.name} [${c.id}]`).join('; ')}`
    : `No modeled structure matches "${structure}" in ${model}.`);
}

function resolvePieces(terms, model) {
  const atlas = catalogue(model);
  return [...new Set(terms.flatMap(term => {
    const part = atlas.parts.find(part => part.id === term);
    return part ? [part.id] : resolveAnatomy(term, model).flatMap(concept => concept.elements);
  }))];
}
function fraction(value, name) {
  if (typeof value !== 'number' || !Number.isFinite(value) || value < 0 || value > 1) throw new Error(`${name} must be between 0 and 1.`);
}
export function anatomyView({structure, structures = [], model = 'male-detail', view = 'three-quarter', context = 0.18, controls, hierarchy = 'systems', systems, depthHidden = [], hidden = [], region = 'all', explode = 0, skinOpacity, labels = true, isolate = false, section, camera, focus, rotate}) {
  if (!VIEWS.includes(view)) throw new Error(`Unknown view: ${view}`);
  fraction(context, 'Context');
  fraction(explode, 'Explode');
  if (rotate && explode > .04) throw new Error('Auto rotation is available only when the model is assembled.');
  if (isolate && explode > 0) throw new Error('Isolation and Explode cannot be active together.');
  if (skinOpacity !== undefined) fraction(skinOpacity, 'Skin opacity');
  if (!HIERARCHIES.includes(hierarchy) || hierarchy === 'depth' && model === 'cell') throw new Error(`Unknown hierarchy for ${model}: ${hierarchy}`);
  if (!REGIONS.includes(region) || region !== 'all' && model === 'cell') throw new Error(`Unknown region for ${model}: ${region}`);
  if (!Array.isArray(structures) || structures.some(item => typeof item !== 'string')) throw new Error('Structures must be a list of names or IDs.');
  if (!Array.isArray(hidden) || hidden.some(item => typeof item !== 'string')) throw new Error('Hidden pieces must be a list of names or IDs.');
  if (!Array.isArray(depthHidden) || depthHidden.some(id => !DEPTH_LAYERS.includes(id)) || depthHidden.length && model === 'cell') throw new Error('Unknown or unavailable depth layer.');
  if (section && (!SECTION_AXES.includes(section.axis) || typeof section.position !== 'number' || typeof section.flip !== 'boolean')) throw new Error('Invalid cross-section.');
  if (section) fraction(section.position, 'Section position');
  if (camera && (!Array.isArray(camera) || ![6, 8].includes(camera.length) || camera.some(n => typeof n !== 'number' || !Number.isFinite(n) || Math.abs(n) >= 1000) || Math.hypot(camera[0]-camera[3],camera[1]-camera[4],camera[2]-camera[5]) <= .001)) throw new Error('Camera must contain a valid position and target (six or eight numbers).');
  if (controls !== undefined && (!Array.isArray(controls) || controls.some(control => !EMBED_CONTROLS.includes(control)))) throw new Error('Unknown embed control.');
  const atlas = catalogue(model);
  const terms = [...(structure ? [structure] : []), ...structures];
  const matches = terms.flatMap(term => {
    const part = atlas.parts.find(part => part.id === term);
    return part ? [{id: part.id, name: part.name, elements: [part.id]}] : resolveAnatomy(term, model);
  });
  const selected = [...new Set(matches.flatMap(concept => concept.elements))];
  if (isolate && !selected.length) throw new Error('Isolation requires a selected structure.');
  const parts = new Map(atlas.parts.map(part => [part.id, part]));
  const availableSystems = new Set(atlas.parts.map(part => part.system));
  const visibleSystems = systems ?? (model==='cell'||model==='embryo'?[...availableSystems]:DEFAULT_LAYERS);
  if (!Array.isArray(visibleSystems) || visibleSystems.some(id => !availableSystems.has(id))) throw new Error('Unknown system for this model.');
  const hiddenPieces = resolvePieces(hidden, model);
  const url = new URL(model.startsWith('local-')?'http://localhost:3016/':VIEWER_URL);
  const params = url.searchParams;
  params.set('model', model);
  params.set('view', view);
  matches.forEach(concept => params.append('select', concept.id));
  params.set('layers', [...new Set(visibleSystems)].join(','));
  params.set('context', String(context));
  params.set('skin', String(skinOpacity ?? (model==='cell'?.18:model==='male-full'?1:model.startsWith('local-')?0:.1)));
  params.set('region', region);
  if (hierarchy !== 'systems') params.set('tree', hierarchy);
  if (depthHidden.length) params.set('depth', [...new Set(depthHidden)].join(','));
  hiddenPieces.forEach(id => params.append('hide', id));
  if (explode) params.set('explode', String(explode));
  if (rotate) params.set('rotate', '1');
  if (!labels) params.set('labels', '0');
  if (isolate) params.set('isolate', '1');
  if (section) { params.set('cut', section.axis); params.set('slice', String(section.position)); if (section.flip) params.set('flip', '1'); }
  if (camera) params.set('camera', camera.join(','));
  if (focus ?? selected.length > 0) params.set('focus', '1');
  params.set('embed', '1');
  if (controls !== undefined) params.set('ui', [...new Set(controls)].join(','));
  const iframe = `<iframe src="${url.href.replaceAll('&', '&amp;').replaceAll('"', '&quot;')}" title="Human Atlas interactive anatomy viewer" loading="lazy" style="width:100%;height:600px;border:0" allowfullscreen></iframe>`;
  return {
    model,
    structures: matches.map(concept => ({id: concept.id, name: concept.name, pieces: concept.elements.length})),
    systems: [...new Set(selected.map(id => parts.get(id)?.system).filter(Boolean))],
    selectedPieces: selected,
    hiddenPieces,
    url: url.href,
    iframe,
  };
}
