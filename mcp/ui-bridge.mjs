export const APP_PROTOCOL = '2026-01-26';
export function atlasResultUrl(result) {
 if (result?.isError || typeof result?.structuredContent?.url !== 'string') return null;
 try {
  const url = new URL(result.structuredContent.url);
  const published = url.origin === 'https://ctzurcanu.github.io' && url.pathname === '/human-atlas/';
  const local = url.origin === 'http://localhost:3016' && url.pathname === '/';
  if ((!published && !local) || url.searchParams.get('embed') !== '1') return null;
  if (published && url.searchParams.get('model')?.startsWith('local-')) return null;
  return url.href;
 } catch { return null; }
}
/** Standard MCP Apps handshake, shared by Claude and other hosts.
 * @param {Window} win
 * @param {{show:(url:string)=>unknown,error:(text:string)=>void,context?:(value:any)=>void,teardown?:()=>void}} handlers
 */
export function connectUiHost(win, {show, error, context = () => {}, teardown = () => {}}) {
 const send = message => win.parent.postMessage({jsonrpc: '2.0', ...message}, '*');
 let initialized = false;
 const receive = event => {
  if (event.source !== win.parent || event.data?.jsonrpc !== '2.0') return;
  const message = event.data;
  if (message.id === 1 && !message.method && !initialized) {
   if (message.error) { error('The host could not initialize the anatomy view.'); return; }
   if (!message.result) return;
   initialized = true;
   context(message.result.hostContext ?? {});
   send({method: 'ui/notifications/initialized', params: {}});
   send({method: 'ui/notifications/size-changed', params: {height: 610}});
  }
  if (message.method === 'ui/notifications/host-context-changed') context(message.params ?? {});
  if (message.method === 'ui/notifications/tool-result') {
   if (message.params?.isError) { error(message.params.content?.find(item => item.type === 'text')?.text ?? 'Unable to show this anatomy view.'); return; }
   const url = atlasResultUrl(message.params);
   if (url) show(url);
  }
  if (message.method === 'ui/notifications/tool-cancelled') error('The anatomy view was cancelled.');
  if (message.method === 'ui/resource-teardown' && message.id !== undefined) {
   teardown(); send({id: message.id, result: {}});
  }
 };
 win.addEventListener('message', receive);
 send({id: 1, method: 'ui/initialize', params: {
  protocolVersion: APP_PROTOCOL, appInfo: {name: 'Atlas', version: '0.2.0'},
  appCapabilities: {availableDisplayModes: ['inline', 'fullscreen']},
 }});
 return () => win.removeEventListener('message', receive);
}
