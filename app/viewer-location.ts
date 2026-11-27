/** MCP Apps run on the host's sandbox origin; share links belong to the atlas. */
let appUrl: string | undefined;
export const viewerHref = () => appUrl ?? location.href;
export function setAppViewerUrl(url: string) { appUrl = new URL(url).href; }
export function navigateViewerUrl(url: string, replace = false) {
 if (appUrl !== undefined) appUrl = new URL(url, appUrl).href;
 else if (replace) history.replaceState(null, '', url);
 else history.pushState(null, '', url);
}

/** Routing and embedding belong to the viewer, not its serialized mesh state. */
export function preserveViewerRouting(url:URL,current:URLSearchParams){
 for(const key of ['atlasSource','blood','nerves','fascia','lymph','dmta','dmtm','cardiac','cardiacRefStage','bowel','hepatic','biliary','venous','cricoid','cricotracheal','caps','delivery','review','embed','ui','compareSession','compareSide','connect','connectJoin','relay']){
  if(!current.has(key))continue;url.searchParams.delete(key);for(const value of current.getAll(key))url.searchParams.append(key,value);
 }
}

/** Coalesce live pose updates without navigation, reload or history entries. */
export function createViewUrlSync(options:{getUrl:()=>string|null;replace:(url:string)=>void;setTimeout:(callback:()=>void,delay:number)=>number;clearTimeout:(id:number)=>void;delayMs?:number}){
 let timer:number|undefined,disposed=false,last='';
 const flush=()=>{if(timer!==undefined){options.clearTimeout(timer);timer=undefined;}if(disposed)return;const url=options.getUrl();if(url&&url!==last){options.replace(url);last=url;}};
 return {schedule(){if(!disposed&&timer===undefined)timer=options.setTimeout(flush,options.delayMs??300);},flush,dispose(){disposed=true;if(timer!==undefined)options.clearTimeout(timer);timer=undefined;}};
}
