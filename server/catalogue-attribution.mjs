import {readFile} from 'node:fs/promises';
import {itemDescription} from '../shared/item-description.mjs';

const escapeExpression = value => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

/** Read source labels and original links from the attribution document. */
export async function catalogueAttribution(file) {
 const markdown = await readFile(file, 'utf8');
 const table = /<!-- atlas-source-references\n([\s\S]*?)\n-->/u.exec(markdown)?.[1];
 if (!table) return [];
 return table.split('\n').map(line => {
  const [id, name, label, url] = line.split('|').map(item => item.trim());
  if (!id || !name || !label || !url) throw new Error('Incomplete atlas source reference');
  const identifier = label.split(' ').map(word => word[0].toUpperCase() + word.slice(1)).join('');
  return {id, name, label, url, identifier, slug: label.replaceAll(' ', '-'),
   keyPattern: new RegExp(escapeExpression(name), 'giu'),
   term: new RegExp(`(?<![A-Za-z0-9])${escapeExpression(name)}(?![A-Za-z0-9])`, 'giu')};
 });
}

const technicalFields = new Set([
 'id', 'conceptId', 'sourceId', 'parentId', 'ownerId', 'partId', 'system',
 'replacedBy', 'replacementParts', 'elements', 'defaultHidden', 'sourceGroups', 'tags',
]);
const assetFields = new Set(['url', 'gzip', 'map', 'normalMap', 'heightMap', 'path', 'atlas', 'asset']);

/** Move descriptive provenance to attribution links while retaining model identity. */
export function projectAttributionFields(value, sources) {
 function project(value, key) {
  if (technicalFields.has(key)) return value;
  if (Array.isArray(value)) return value.map(item => project(item, key));
  if (value && typeof value === 'object') {
   return Object.fromEntries(Object.entries(value).map(([name, item]) => {
    let outputKey = name;
    if (name.includes('public/models/')) {
     for (const {term, slug} of sources) outputKey = outputKey.replace(term, slug);
    } else {
     for (const {keyPattern, label, identifier} of sources) {
      outputKey = outputKey.replace(keyPattern, original => {
       if (!/^[A-Za-z0-9]+$/u.test(name)) return label;
       return /^[A-Z]/u.test(original) ? identifier : identifier[0].toLowerCase() + identifier.slice(1);
      });
     }
    }
    return [outputKey, project(item, name)];
   }));
  }
  if (typeof value !== 'string') return value;
  if (key === 'description' || key === 'coverageLimitation') value = itemDescription(value);
  if (/^(https?:)?\/\//iu.test(value)) {
   const target = new URL(value, 'https://localhost');
   const matched = sources.find(({url}) => {
    const source = new URL(url);
    const prefix = source.pathname.replace(/\/$/u, '');
    return source.hostname === target.hostname &&
     (!prefix || target.pathname === prefix || target.pathname.startsWith(prefix + '/'));
   });
   if (matched) return `/ATTRIBUTION.md#${matched.id}`;
   return value;
  }
  // These values locate the existing geometry and textures.
  if (/^\/(local-models|models|assets)\//u.test(value) || assetFields.has(key) && value.startsWith('/')) {
   for (const {term, slug} of sources) value = value.replace(term, slug);
   return value;
  }
  for (const {term, label, slug} of sources) {
   value = value.replace(term, (original, index) => key === 'schema' || key === 'sourceFile'
    ? slug : index === 0 && /^[A-Z]/u.test(original) ? label[0].toUpperCase() + label.slice(1) : label);
  }
  return value;
 }
 return project(value, '');
}

export function projectCatalogueAttribution(catalogue, sources) {
 if (!Array.isArray(catalogue?.parts) || !Array.isArray(catalogue.chunks)) return catalogue;
 return projectAttributionFields(catalogue, sources);
}

/** Resolve neutral asset URLs against existing local import directories. */
export function localAttributionPath(name, sources) {
 for (const {slug, name: original} of sources) {
  name = name.replaceAll(slug, original.toLowerCase());
 }
 return name;
}
