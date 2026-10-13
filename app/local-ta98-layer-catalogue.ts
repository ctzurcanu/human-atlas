/** Local finite bowel layers are included in ordinary TA98 views. Source-only
 * and delivery diagnostics remain explicitly selectable; public models do not
 * load local donor assets. */
export function localTa98LayerCatalogue(model:string,query:URLSearchParams,development:boolean):string|undefined {
 if(!development||!['local-ta98','local-female-ta98'].includes(model))return undefined;
 if(query.get('caps')==='shared-v28')return `/local-models/ta98-review/colonic-shared-caps-v28/${model==='local-ta98'?'male':'female'}/atlas.json`;
 if(model==='local-ta98'&&query.get('venous')==='confluence')return '/local-models/ta98-review/male-portal-confluence/atlas.json';
 if(model==='local-ta98'&&query.get('venous')==='portal-wall')return '/local-models/ta98-review/male-portal-wall/atlas.json';
 if(model==='local-female-ta98'&&query.get('biliary')==='gastric-vein')return '/local-models/ta98-review/female-gastric-vein-repair/atlas.json';
 if(model==='local-female-ta98'&&query.get('biliary')==='hepatic-wall')return '/local-models/ta98-review/female-common-hepatic-wall/atlas.json';
 if(model==='local-female-ta98'&&query.get('biliary')==='native-tangent')return '/local-models/ta98-review/female-native-biliary-tangent/atlas.json';
 if(model==='local-female-ta98'&&query.get('biliary')==='native')return '/local-models/ta98-review/female-native-biliary-atlas/atlas.json';
 if(model==='local-female-ta98'&&query.get('biliary')==='confluence')return '/local-models/ta98-review/female-biliary-confluence-joined/atlas.json';
 if(model==='local-female-ta98'&&query.get('biliary')==='regions')return '/local-models/ta98-review/female-gallbladder-regions/atlas.json';
 if(query.get('biliary')==='folds')return `/local-models/ta98-review/cystic-folds/${model==='local-ta98'?'male':'female'}/atlas.json`;
 const mode=query.get('bowel');
 if(mode==='source'||(!mode&&query.has('delivery')))return undefined;
 if(model==='local-female-ta98'&&query.get('hepatic')==='coats'&&mode!=='prototype')
  return '/local-models/ta98-review/gallbladder-coats-preview/female/atlas.json';
 if(model==='local-female-ta98'&&query.get('hepatic')==='assembly'&&mode!=='prototype')
  return '/local-models/ta98-review/hepatic-biliary-preview/female/atlas.json';
 // This opt-in whole-liver candidate keeps the ordinary atlas identities and
 // controls. Never substitute local, unaccepted geometry in public delivery.
 if(model==='local-female-ta98'&&query.get('hepatic')==='repair'&&mode!=='prototype')
  return '/local-models/ta98-review/hepatic-boundary-preview/female/atlas.json';
 const folder=mode==='prototype'?'bowel-wall-main-prototype':'bowel-wall-local-atlas';
 return `/local-models/ta98-review/${folder}/${model==='local-ta98'?'male':'female'}/atlas.json`;
}

/** Preserve a selected local review when copying a view; public catalogues
 * continue to ignore this development-only routing control. */
export function preserveLocalTa98LayerCatalogueQuery(url:URL,current:URLSearchParams,model:string,development:boolean):void{
 if(development&&['local-ta98','local-female-ta98'].includes(model)&&current.get('caps')==='shared-v28')url.searchParams.set('caps','shared-v28');
 if(development&&model==='local-ta98'&&current.get('venous')==='confluence')url.searchParams.set('venous','confluence');
 if(development&&model==='local-ta98'&&current.get('venous')==='portal-wall')url.searchParams.set('venous','portal-wall');
}
