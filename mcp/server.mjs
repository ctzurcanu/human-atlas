#!/usr/bin/env node
import {pathToFileURL} from "node:url";
import {StdioServerTransport} from "@modelcontextprotocol/sdk/server/stdio.js";
import "./atlas-data.mjs";
import {createServer as createCoreServer} from "./server-core.mjs";
import {buildMcpUi,buildIframeUi} from "./build-ui.mjs";
export const createServer=async()=>{const [uiHtml,iframeHtml]=await Promise.all([buildMcpUi(),buildIframeUi()]);return createCoreServer({uiHtml,iframeHtml});};
if(process.argv[1]&&import.meta.url===pathToFileURL(process.argv[1]).href)await (await createServer()).connect(new StdioServerTransport());
