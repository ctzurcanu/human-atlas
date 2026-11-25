/** MCP Apps run on the host's sandbox origin; share links belong to the atlas. */
let appUrl: string | undefined;
export const viewerHref = () => appUrl ?? location.href;
export function setAppViewerUrl(url: string) { appUrl = new URL(url).href; }
export function navigateViewerUrl(url: string, replace = false) {
 if (appUrl !== undefined) appUrl = new URL(url, appUrl).href;
 else if (replace) history.replaceState(null, '', url);
 else history.pushState(null, '', url);
}

/** Keep one history entry current during continuous camera gestures and animations. */
export function watchViewUrl(getUrl:()=>string|null){
 let timer:number|undefined;
 const flush=()=>{
  if(timer===undefined)return;
  window.clearTimeout(timer);
  timer=undefined;
  const url=getUrl();
  if(url&&url!==viewerHref())navigateViewerUrl(url,true);
 };
 const schedule=()=>{
  if(timer===undefined)timer=window.setTimeout(flush,250);
 };
 const onVisibility=()=>{if(document.visibilityState==='hidden')flush();};
 window.addEventListener('pointerup',flush);
 window.addEventListener('keyup',flush);
 window.addEventListener('pagehide',flush);
 document.addEventListener('visibilitychange',onVisibility);
 return {schedule,dispose(){
  window.clearTimeout(timer);
  timer=undefined;
  window.removeEventListener('pointerup',flush);
  window.removeEventListener('keyup',flush);
  window.removeEventListener('pagehide',flush);
  document.removeEventListener('visibilitychange',onVisibility);
 }};
}
