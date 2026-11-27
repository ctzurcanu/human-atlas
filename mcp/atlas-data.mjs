import {readFileSync} from 'node:fs';
import {configureAtlasData} from './atlas-core.mjs';
import {catalogueRevision,readCatalogue} from './read-catalogue.mjs';
configureAtlasData(readCatalogue,JSON.parse(readFileSync(new URL('../app/data/ta98-metadata.json',import.meta.url),'utf8')),catalogueRevision);
export * from './atlas-core.mjs';
