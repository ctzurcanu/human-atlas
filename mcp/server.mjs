#!/usr/bin/env node
import {pathToFileURL} from "node:url";
import {StdioServerTransport} from "@modelcontextprotocol/sdk/server/stdio.js";
import "./atlas-data.mjs";
import {createServer as createCoreServer} from "./server-core.mjs";
import {buildMcpUi} from "./build-ui.mjs";
export const createServer=async()=>createCoreServer({uiHtml:await buildMcpUi()});
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href)await (await createServer()).connect(new StdioServerTransport());
