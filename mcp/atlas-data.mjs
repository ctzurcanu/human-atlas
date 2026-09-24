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

export function anatomyView({structure, model = 'male-detail', view = 'three-quarter', context = 0.18, controls}) {
  if (!VIEWS.includes(view)) throw new Error(`Unknown view: ${view}`);
  if (typeof context !== 'number' || !Number.isFinite(context) || context < 0 || context > 1) throw new Error('Context must be between 0 and 1.');
  if (controls !== undefined && (!Array.isArray(controls) || controls.some(control => !EMBED_CONTROLS.includes(control)))) throw new Error('Unknown embed control.');
  const matches = resolveAnatomy(structure, model);
  const atlas = catalogue(model);
  const selected = [...new Set(matches.flatMap(concept => concept.elements))];
  const parts = new Map(atlas.parts.map(part => [part.id, part]));
  const url = new URL(model.startsWith('local-')?'http://localhost:3016/':VIEWER_URL);
  const params = url.searchParams;
  params.set('model', model);
  params.set('view', view);
  matches.forEach(concept => params.append('select', concept.id));
  params.set('layers', (model==='cell'||model==='embryo'?[...new Set(atlas.parts.map(part=>part.system))]:DEFAULT_LAYERS).join(','));
  params.set('context', String(context));
  params.set('skin', model==='cell'?'0.18':model.startsWith('local-')?'0':'0.1');
  params.set('region', 'all');
  params.set('peel', '0');
  params.set('focus', '1');
  params.set('embed', '1');
  if (controls !== undefined) params.set('ui', [...new Set(controls)].join(','));
  const iframe = `<iframe src="${url.href.replaceAll('&', '&amp;').replaceAll('"', '&quot;')}" title="Human Atlas interactive anatomy viewer" loading="lazy" style="width:100%;height:600px;border:0" allowfullscreen></iframe>`;
  return {
    model,
    structures: matches.map(concept => ({id: concept.id, name: concept.name, pieces: concept.elements.length})),
    systems: [...new Set(selected.map(id => parts.get(id)?.system).filter(Boolean))],
    url: url.href,
    iframe,
  };
}
