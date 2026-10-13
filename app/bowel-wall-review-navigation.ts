/** Local review variants remain separate from the ordinary TA98 source views. */
export const BOWEL_WALL_REVIEW_OPTIONS=[
 {value:'bowel-wall-male',label:'Male · Small intestine wall review',sex:'male'},
 {value:'bowel-wall-female',label:'Female · Small intestine wall review',sex:'female'},
] as const;
export function bowelWallReviewUrl(current:string,sex:'male'|'female'):string{
 const url=new URL(current);url.search='';url.hash='';
 url.searchParams.set('model',sex==='male'?'local-ta98':'local-female-ta98');
 url.searchParams.set('bowel','prototype');url.searchParams.set('tree','guest:ta98');
 url.searchParams.set('select',sex==='male'?'ZA:Jejunum':'LOCAL:female:small_intestine_jejunum');
 url.searchParams.set('focus','1');url.searchParams.set('skin','0');
 url.searchParams.set('context','0');url.searchParams.set('isolate','1');
 return url.href;
}
export function ordinaryModelUrl(current:string,model:string):string{
 const url=new URL(current);url.search='';url.hash='';url.searchParams.set('model',model);return url.href;
}
