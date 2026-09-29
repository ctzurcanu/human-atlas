import {readFileSync,mkdirSync,writeFileSync} from 'node:fs';
import {MODELS} from '../mcp/atlas-core.mjs';
const catalogues={};
const read=path=>JSON.parse(readFileSync(new URL(path,import.meta.url),'utf8'));
for(const file of new Set(Object.entries(MODELS).filter(([id])=>!id.startsWith('local-')).map(([,file])=>file))){
 const atlas=read(`../public/models/${file}`);
 catalogues[file]={parts:atlas.parts.map(({id,name,system})=>({id,name,system})),concepts:atlas.concepts.map(({id,name,elements})=>({id,name,elements}))};
}
const terms=read('../app/data/ta98-metadata.json').byConcept;
const normalize=value=>value.trim().toLocaleLowerCase().replace(/\s+/g,' ');
for(const atlas of Object.values(catalogues))for(const concept of atlas.concepts){
 const t=terms[concept.id]??(/^FMA:?\d+$/.test(concept.id)?{fma:`FMA:${concept.id.replace(/^FMA:?/,'')}`} :{});
 concept.normalizedId=normalize(concept.id);concept.normalizedName=normalize(concept.name);
 concept.plainName=normalize(concept.name.replace(/\s*\((?:left|right)\)$/i,''));
 concept.searchTerms=[concept.normalizedName,concept.normalizedId,...[t.ta98,t.ta98&&`TA98:${t.ta98}`,t.tha,t.fma,t.ontology,t.latin,t.latin&&`La:${t.latin}`].filter(Boolean).map(normalize)];
}
const ids=new Set(Object.values(catalogues).flatMap(a=>a.concepts.map(c=>c.id)));
const data={catalogues,terminology:{byConcept:Object.fromEntries(Object.entries(terms).filter(([id])=>ids.has(id)))}};
const directory=new URL('../cloudflare/generated/',import.meta.url);
mkdirSync(directory,{recursive:true});
writeFileSync(new URL('mcp-data.json',directory),JSON.stringify(data));
console.log(`Cloudflare MCP catalogue: ${Object.keys(catalogues).length} models, ${(JSON.stringify(data).length/1024/1024).toFixed(2)} MiB`);
