import {readFileSync,mkdirSync,writeFileSync} from 'node:fs';
import {readCatalogue} from '../mcp/read-catalogue.mjs';
import {MODELS,terminologyValues} from '../mcp/atlas-core.mjs';
import {buildMcpUi} from '../mcp/build-ui.mjs';
const catalogues={};
const read=path=>JSON.parse(readFileSync(new URL(path,import.meta.url),'utf8'));
for(const [model,file] of Object.entries(MODELS)){
 if(model.startsWith('local-')||catalogues[file])continue;
 const atlas=readCatalogue(model,file);
 catalogues[file]={parts:atlas.parts.map(({id,name,system})=>({id,name,system})),concepts:atlas.concepts.map(({id,name,elements,ta98Term,ta98Kind})=>({id,name,elements,ta98Term,ta98Kind}))};
}
const terms=read('../app/data/ta98-metadata.json');
const normalize=value=>value.trim().toLocaleLowerCase().replace(/\s+/g,' ');
// A concept can share thousands of pieces with other concepts. Store their indices once
// in the Worker bundle and keep one pre-normalized search string per concept.
for(const atlas of Object.values(catalogues)){
 const indices=new Map(atlas.parts.map((part,i)=>[part.id,i]));
 for(const concept of atlas.concepts){
  concept.elements=concept.elements.map(id=>indices.get(id)??id);
  concept.searchText=[concept.name,concept.id,...terminologyValues(concept,terms)].map(normalize).join('\u0001');
 }
}
const ids=new Set(Object.values(catalogues).flatMap(a=>a.concepts.map(c=>c.id)));
const codes=new Set(Object.values(catalogues).flatMap(a=>a.concepts.map(c=>c.ta98Term).filter(Boolean)));
const data={catalogues,terminology:{byConcept:Object.fromEntries(Object.entries(terms.byConcept).filter(([id])=>ids.has(id))),byCode:Object.fromEntries(Object.entries(terms.byCode).filter(([code])=>codes.has(code)))}};
const directory=new URL('../cloudflare/generated/',import.meta.url);
mkdirSync(directory,{recursive:true});
writeFileSync(new URL('mcp-data.json',directory),JSON.stringify(data));
console.log(`Cloudflare MCP catalogue: ${Object.keys(catalogues).length} models, ${(JSON.stringify(data).length/1024/1024).toFixed(2)} MiB`);
await buildMcpUi();
