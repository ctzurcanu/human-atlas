import {createHash} from 'node:crypto';
import {mkdir, readFile, rename, writeFile} from 'node:fs/promises';
import {extname, resolve, sep} from 'node:path';
import {gzipSync, gunzipSync} from 'node:zlib';
import {catalogueAttribution, localAttributionPath, projectCatalogueAttribution} from './catalogue-attribution.mjs';

const root = resolve(import.meta.dirname, '..');
const source = '.local-models/geometry-workbench/dmta-v1/dmtm-v8-skin-clearance-v1/circulation-atlas.json';
const prefix = '/models/male-ta98/';
const target = resolve(root, 'outputs/male-ta98-pages' + prefix);
const limit = 16 * 1024 * 1024;
const sha = bytes => createHash('sha256').update(bytes).digest('hex');
const references = await catalogueAttribution(resolve(root, 'ATTRIBUTION.md'));
const sourceBytes = await readFile(resolve(root, source));
// Use the same labels and source identities as the original development URL.
const atlas = projectCatalogueAttribution(JSON.parse(sourceBytes), references);
const widths = {positions:12, normals:6, uvs:8, colors:3, indices:4, normalSigns:1, gastricLumen:1, tissueChart:6};
const sourceChunks = atlas.chunks;
const sourceParts = {};
const chunks = [];
const systems = {};
const bySource = Map.groupBy(atlas.parts, part => part.chunk);
const sectionBanks = new Set([atlas.liverSectionIntegration?.bankUrl, atlas.liverCompletionIntegration?.functionalSections?.bankUrl]);

async function readAsset(url) {
 const local = url.startsWith('/local-models/');
 if (!local && !url.startsWith('/models/') && !url.startsWith('/assets/')) throw Error('Unexpected model asset ' + url);
 const dir = resolve(root, local ? '.local-models' : 'public');
 const name = url.split('?')[0].slice(local ? '/local-models/'.length : 1);
 for (const candidate of new Set([name, localAttributionPath(name, references)])) {
  const file = resolve(dir, candidate);
  if (!file.startsWith(dir + sep)) throw Error('Invalid asset path ' + url);
  try { return await readFile(file); } catch (error) { if (error.code !== 'ENOENT') throw error; }
 }
 throw Error('Missing model asset ' + url);
}

async function store(relative, bytes) {
 const file = resolve(target, relative);
 await mkdir(resolve(file, '..'), {recursive:true});
 let stored;
 try {
  stored = await readFile(file);
 } catch (error) {
  if (error.code !== 'ENOENT') throw error;
  const temp = file + '.' + process.pid + '.tmp';
  await writeFile(temp, bytes);
  await rename(temp, file);
  stored = await readFile(file);
 }
 if (!stored.equals(bytes)) throw Error('Incomplete or conflicting model asset ' + relative);
 return prefix + relative;
}

async function writeBank(system, bytes) {
 const hash = sha(bytes);
 const compressed = gzipSync(bytes, {level:4});
 if (!gunzipSync(compressed).equals(bytes)) throw Error('Geometry compression failed');
 const url = await store('geometry/' + system + '/' + hash + '.bin.gz', compressed);
 chunks.push({system, url, gzip:url, bytes:bytes.length, gzipBytes:compressed.length, sha256:hash});
}

for (const [index, parts] of bySource) {
 const entry = sourceChunks[index];
 const raw = await readAsset(entry.url);
 if (raw.length !== entry.bytes || entry.sha256 && sha(raw) !== entry.sha256) throw Error('Source bank mismatch ' + entry.url);
 const groups = Map.groupBy(parts, p => p.system);
 for (const [system, rows] of groups) {
  if (!/^[a-z-]+$/.test(system)) throw Error('Unknown system folder ' + system);
  systems[system] = (systems[system] ?? 0) + rows.length;
  const retainWhole = sectionBanks.has(entry.url);
  if (retainWhole && groups.size !== 1) throw Error('Section classifier spans multiple systems');
  let ranges = new Map(), buffers = [], size = 0;
  async function flush() {
   if (size) await writeBank(system, Buffer.concat(buffers, size));
   ranges = new Map(); buffers = []; size = 0;
  }
  for (const part of rows) {
   const layout = {chunk:index};
   const attributes = [];
   for (const [key, width] of Object.entries(widths)) {
    if (!(key in part)) continue;
    const start = part[key], length = width * (key === 'indices' ? part.indexCount : part.vertexCount);
    if (!Number.isSafeInteger(start) || !Number.isSafeInteger(length) || start < 0 || length < 0 || start + length > raw.length) throw Error('Invalid ' + key + ' for ' + part.id);
    layout[key] = start;
    attributes.push({key, start, length, range:start + '/' + length});
   }
   if (Object.hasOwn(sourceParts, part.id)) throw Error('Duplicate part ' + part.id);
   sourceParts[part.id] = layout;
   if (!retainWhole) {
    const needed = attributes.reduce((n, field) => n + (ranges.has(field.range) ? 0 : field.length + 3), 0);
    if (size && size + needed > limit) await flush();
    for (const {key, start, length, range} of attributes) {
     if (!ranges.has(range)) {
      const pad = (4 - size % 4) % 4;
      buffers.push(Buffer.alloc(pad), raw.subarray(start, start + length));
      size += pad; ranges.set(range, size); size += length;
     }
     part[key] = ranges.get(range);
    }
   }
   part.chunk = chunks.length;
  }
  if (retainWhole) await writeBank(system, raw); else await flush();
 }
}

const relocated = new Map();
async function textures(obj) {
 for (const [key, value] of Object.entries(obj)) {
  if (value && typeof value === 'object') await textures(value);
  else if (['map','normalMap','heightMap'].includes(key) && typeof value === 'string') {
   if (!relocated.has(value)) {
    const data = await readAsset(value);
    relocated.set(value, await store('textures/' + sha(data) + extname(value.split('?')[0]), data));
   }
   obj[key] = relocated.get(value);
  }
 }
}
await textures(atlas.materials);
for (const {hairField} of atlas.parts) if (hairField) {
 const raw = await readAsset(hairField.url);
 if (sha(raw) !== hairField.sha256 || raw.length !== hairField.bytes) throw Error('Invalid hair field');
 const path = 'geometry/integumentary/' + sha(raw) + '.bin';
 hairField.url = await store(path, raw);
 const compressed = gzipSync(raw, {level:4});
 await store(path + '.gz', compressed); hairField.gzipBytes = compressed.length;
}
atlas.name = 'Male TA98';
atlas.chunks = chunks;
atlas.partsCount = atlas.parts.length;
// Selection and label registries authenticate original bank coordinates. These
// are evidence for those checks, independent of the system download offsets.
atlas.sourceGeometry = {chunks:sourceChunks, parts:sourceParts};
atlas.systems = Object.entries(systems).sort().map(([id,items]) => ({
 id, items, chunks:chunks.flatMap((c,i) => c.system === id ? [i] : []),
}));
atlas.delivery = {...atlas.delivery, mode:'on-demand', sourceAtlasSha256:sha(sourceBytes)};
const json = JSON.stringify(atlas).replaceAll(root + '/', '').replace(/\/Users\/[^/"\\]+\//g, '~/');
const temporary = resolve(target, 'atlas.' + process.pid + '.tmp');
await writeFile(temporary, json);
await rename(temporary, resolve(target, 'atlas.json'));
console.log('Male TA98:', atlas.parts.length, 'parts,', atlas.systems.length, 'systems,', chunks.length, 'geometry banks in', target);
