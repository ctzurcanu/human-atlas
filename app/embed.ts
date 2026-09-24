export function embedUrl(viewUrl:string){
 const url=new URL(viewUrl);
 url.searchParams.set('embed','1');
 return url.href;
}

export function iframeCode(viewUrl:string){
 const src=embedUrl(viewUrl).replaceAll('&','&amp;').replaceAll('"','&quot;').replaceAll('<','&lt;');
 return `<iframe src="${src}" title="Human Atlas interactive anatomy viewer" loading="lazy" style="width:100%;height:600px;border:0" allowfullscreen></iframe>`;
}
