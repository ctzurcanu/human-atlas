/** MCP Apps run on the host's sandbox origin; share links belong to the atlas. */
let appUrl: string | undefined;
export const viewerHref = () => appUrl ?? location.href;
export function setAppViewerUrl(url: string) { appUrl = new URL(url).href; }
export function navigateViewerUrl(url: string, replace = false) {
 if (appUrl !== undefined) appUrl = new URL(url, appUrl).href;
 else if (replace) history.replaceState(null, '', url);
 else history.pushState(null, '', url);
}
