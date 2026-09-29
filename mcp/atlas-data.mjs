import {readFileSync} from "node:fs";
import {configureAtlasData} from "./atlas-core.mjs";
configureAtlasData((model,filename)=>JSON.parse(readFileSync(new URL((model.startsWith("local-")?"../.local-models/":"../public/models/")+filename,import.meta.url),"utf8")),JSON.parse(readFileSync(new URL("../app/data/ta98-metadata.json",import.meta.url),"utf8")));
export * from "./atlas-core.mjs";
