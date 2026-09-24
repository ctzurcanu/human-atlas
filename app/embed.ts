export const EMBED_UI=[
 {id:'model',label:'Model picker'},
 {id:'search',label:'Search'},
 {id:'study',label:'Study tools'},
 {id:'systems',label:'Systems'},
 {id:'camera',label:'Camera controls'},
 {id:'explode',label:'Explode and reset dock'},
 {id:'details',label:'Structure details'},
 {id:'open',label:'Open full viewer'},
 {id:'download',label:'Download PNG'},
] as const;
export type EmbedUi=typeof EMBED_UI[number]['id'];
export const DEFAULT_EMBED_UI:EmbedUi[]=EMBED_UI.filter(item=>!['explode','study','camera','download'].includes(item.id)).map(item=>item.id);
const valid=new Set<string>(EMBED_UI.map(item=>item.id));

export function readEmbedUi(search:string):Set<EmbedUi>{
 const params=new URLSearchParams(search);
 if(!params.has('ui'))return new Set(DEFAULT_EMBED_UI);
 return new Set(params.get('ui')!.split(',').filter((item):item is EmbedUi=>valid.has(item)));
}

export function embedUrl(viewUrl:string,ui?:readonly EmbedUi[]){
 const url=new URL(viewUrl);
 url.searchParams.set('embed','1');
 if(ui)url.searchParams.set('ui',ui.join(','));
 return url.href;
}

export function iframeCode(viewUrl:string,ui?:readonly EmbedUi[],height=600){
 const src=embedUrl(viewUrl,ui).replaceAll('&','&amp;').replaceAll('"','&quot;').replaceAll('<','&lt;');
 const safeHeight=Number.isFinite(height)?Math.max(340,Math.min(1200,Math.round(height))):600;
 return `<iframe src="${src}" title="Human Atlas interactive anatomy viewer" loading="lazy" style="width:100%;height:${safeHeight}px;border:0" allowfullscreen></iframe>`;
}
