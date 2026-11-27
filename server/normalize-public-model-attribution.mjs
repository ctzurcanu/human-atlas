import {readFile, readdir, rename, stat, writeFile} from 'node:fs/promises';
import {join, resolve, sep} from 'node:path';
import {catalogueAttribution, projectAttributionFields, projectCatalogueAttribution} from './catalogue-attribution.mjs';

const publicDirectory = resolve('public');
const sourceNames = await catalogueAttribution('ATTRIBUTION.md');
const catalogues = [];
const assets = new Map();

async function visit(directory) {
 for (const entry of await readdir(directory, {withFileTypes: true})) {
  const file = join(directory, entry.name);
  if (entry.isDirectory()) {
   await visit(file);
   continue;
  }
  if (!entry.isFile() || !entry.name.endsWith('.json')) continue;
  const original = JSON.parse(await readFile(file, 'utf8'));
  const next = projectCatalogueAttribution(original, sourceNames);
  if (JSON.stringify(original) === JSON.stringify(next)) continue;
  original.chunks.forEach((chunk, index) => {
   for (const field of ['url', 'gzip']) {
    const originalUrl = chunk[field], newUrl = next.chunks[index][field];
    if (originalUrl === newUrl || !originalUrl?.startsWith('/models/')) continue;
    const before = resolve(publicDirectory, originalUrl.slice(1));
    const after = resolve(publicDirectory, newUrl.slice(1));
    if (!before.startsWith(publicDirectory + sep) || !after.startsWith(publicDirectory + sep)) {
     throw new Error('Invalid model asset path');
    }
    assets.set(before, after);
   }
  });
  catalogues.push({file, contents: JSON.stringify(next) + '\n'});
 }
}

await visit(join(publicDirectory, 'models'));
const publicEvidence = join(publicDirectory, 'assets');
for (const entry of await readdir(publicEvidence, {withFileTypes: true})) {
 if (!entry.isFile() || !entry.name.endsWith('-provenance.json')) continue;
 const file = join(publicEvidence, entry.name);
 const original = JSON.parse(await readFile(file, 'utf8'));
 const next = projectAttributionFields(original, sourceNames);
 if (JSON.stringify(original) !== JSON.stringify(next)) {
  catalogues.push({file, contents: JSON.stringify(next) + '\n'});
 }
}
const moves = [];
for (const [before, after] of assets) {
 const missing = error => {
  if (error.code !== 'ENOENT') throw error;
  return null;
 };
 const [original, existing] = await Promise.all([stat(before).catch(missing),stat(after).catch(missing)]);
 if (existing?.isFile()) {
  if (original && (original.size !== existing.size || !(await readFile(before)).equals(await readFile(after)))) {
   throw new Error(`Conflicting geometry at ${after}`);
  }
  continue;
 }
 if (!original?.isFile()) throw new Error(`Expected existing model asset at ${before}`);
 moves.push([before, after]);
}
for (const [before, after] of moves) await rename(before, after);
for (const {file, contents} of catalogues) {
 const temporary = file + '.attribution-tmp';
 await writeFile(temporary, contents, {flag: 'wx'});
 await rename(temporary, file);
}
console.log(`Normalized ${catalogues.length} public model catalogues and ${moves.length} geometry paths`);
