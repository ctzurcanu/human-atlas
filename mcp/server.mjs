#!/usr/bin/env node
import {readFileSync} from "node:fs";
import {pathToFileURL} from "node:url";
import {StdioServerTransport} from "@modelcontextprotocol/sdk/server/stdio.js";
import "./atlas-data.mjs";
import {createServer as createCoreServer} from "./server-core.mjs";
export const createServer=()=>createCoreServer({uiHtml:readFileSync(new URL("./view.html",import.meta.url),"utf8")});
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href)await createServer().connect(new StdioServerTransport());
