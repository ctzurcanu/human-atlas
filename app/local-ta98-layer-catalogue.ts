/** Local finite bowel layers are included in ordinary TA98 views. Source-only
 * and delivery diagnostics remain explicitly selectable; public models do not
 * load local donor assets. Refresh the body/detail presentation catalogues with
 * `npm run build:organ-models` after installing changes in their source atlas. */
/** An explicit local catalogue changes only this review URL, never an installed route. */
export function localAtlasPreviewSource(query:URLSearchParams):string|undefined {
 const source=query.get('atlasSource');
 if(!source||source.length>512||!/^\/local-models\/[A-Za-z0-9_./-]+\.json$/.test(source))return undefined;
 if(source.split('/').some((segment,index)=>index>0&&(!segment||segment==='.'||segment==='..')))return undefined;
 return source;
}

export function localTa98LayerCatalogue(model:string,query:URLSearchParams,development:boolean):string|undefined {
 if(!development||!['local-ta98','local-female-ta98'].includes(model))return undefined;
 if([...query.values()].some(value=>value.startsWith('registered-collection-'))){
  query=new URLSearchParams(query);
  for(const [key,value] of [...query])if(value.startsWith('registered-collection-'))
   query.set(key,value.replace('registered-collection-','biodigital-'));
 }
 const preview=localAtlasPreviewSource(query);if(preview)return preview;
 // Opt-in review catalogues for the registered circulation and nervous system.
  if(model==='local-female-ta98'&&(query.get('circulation')==='biodigital-review'||query.get('nerves')==='biodigital-review'||query.get('lymphatics')==='biodigital-review'))return '/local-models/geometry-workbench/female-circ/female-circ-v2/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('nerves')==='biodigital-review')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/nervous-fit-v7/circulation-atlas.json';
 // BioDigital male fascia, fat and connective tissue fitted to this body (opt-in review, anatomically unaccepted; includes the v17 circulation).
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v4')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v4/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v11')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v11/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v10')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v10/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v9')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v9/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v8')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v8/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v7')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v7/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v6')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v6/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v5')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v5/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-v1')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/fascia-fit-v1/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/best-model-v2/circulation-atlas.json';
 // Opt-in body correction on top of the fascia review: organ/muscle shells pushed out of the body's bone solids (positions only).
 if(model==='local-ta98'&&query.get('fascia')==='biodigital-review-muscles')return '/local-models/geometry-workbench/biodigital-male-nervous-v1/best-model-v2/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('lymph')==='biodigital-review')return '/local-models/geometry-workbench/biodigital-male-lymphatic-v1/lymph-fit-v9/circulation-atlas.json';
 // dmta: same male TA98 bone surfaces, midpoint-subdivided and partitioned into provisional TA98 feature regions (unvalidated).
 // dmtm: dmta plus the muscles partitioned into provisional TA98 part regions (same method).
 if(model==='local-ta98'&&query.get('dmtm')==='1')return '/local-models/geometry-workbench/dmta-v1/dmtm-v8-skin-clearance-v1/circulation-atlas.json';
 if(model==='local-ta98'&&query.get('dmta')==='1')return '/local-models/geometry-workbench/dmta-v1/dmta-v36/circulation-atlas.json';
 // Explicit packed pulmonary-wall review; ordinary and public routes stay unchanged.
 if(model==='local-ta98'&&query.get('blood')==='biodigital-root-review')return '/local-models/geometry-workbench/biodigital-male-circulation-v1/thoracic-skeletal-refit-v17/circulation-atlas.json';
 // Original cardiac source and its coronary assembly are reviewed together.
 // Keep this explicit until the remaining thoracic packing checks pass.
 if(model==='local-ta98'&&query.get('cardiac')==='biodigital-review')return '/local-models/geometry-workbench/biodigital-male-circulation-v1/thoracic-skeletal-refit-v17/circulation-atlas.json';
 // The imported blood assembly is an opt-in fit review. Ordinary male data
 // stays on its installed source until regional packing and interfaces pass.
 // v37 retains every v34 buffer; mismatched duplicate nerves and adductors are hidden
 // repair alternatives, so fitting uses the same anatomy as the displayed body.
 if(model==='local-ta98'&&query.get('blood')==='biodigital-review')return '/local-models/geometry-workbench/biodigital-male-circulation-v1/thoracic-skeletal-refit-v17/circulation-atlas.json';
 // Whole-heart placement review: vessels remain fixed; this candidate has not
 // passed the surrounding lung/diaphragm packing checks and is not installed.
 if(model==='local-ta98'&&query.get('cardiac')==='fixed-vessel-fit')return '/local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1/atlas.json';
 if(model==='local-ta98'&&query.get('cardiac')==='chordal-shape')return '/local-models/geometry-workbench/heart/chordal-shape-view-v205/atlas.json';
 if(model==='local-ta98'&&query.get('cardiac')==='wall-v48')return '/local-models/geometry-workbench/heart/wall-preview-v48/atlas.json';
 if(model==='local-female-ta98'&&query.get('cardiac')==='thoracic-fit')return '/local-models/geometry-workbench/heart/female-roots-surface-repaired-v598/atlas.json';
 if(model==='local-ta98'&&query.get('cardiac')==='thoracic-fit')return '/local-models/geometry-workbench/heart/male-apex-shared-normals-v593/atlas.json';
 if(query.get('cardiac')==='review')return `/local-models/organ-models/spinal-cord/${model==='local-ta98'?'male':'female'}-body.json`;
 if(model==='local-ta98'&&query.get('cardiac')==='wall-v98')return '/local-models/geometry-workbench/heart/regional-bank-v98/atlas.json';
 if(model==='local-ta98'&&query.get('cardiac')==='wall-v94')return '/local-models/geometry-workbench/heart/wall-review-v94/atlas.json';
 if(model==='local-ta98'&&query.get('cardiac')==='wall-v19')return '/local-models/geometry-workbench/heart/wall-preview-v19/atlas.json';
 if(query.get('hepatic')==='completion-final')return model==='local-ta98'?'/local-models/geometry-workbench/liver-completion/final-clearance-v270/male/atlas.json':'/local-models/geometry-workbench/liver-completion/final-clearance-v270/female/atlas.json';
 if(model==='local-female-ta98'&&query.get('hepatic')==='female-transfer-v202')return '/local-models/geometry-workbench/liver-completion/female-transfer-v202/atlas.json';
 if(model==='local-female-ta98'&&query.get('hepatic')==='female-transfer-v203')return '/local-models/geometry-workbench/liver-completion/female-transfer-v203/atlas.json';
 if(model==='local-female-ta98'&&query.get('hepatic')==='female-transfer-v204')return '/local-models/geometry-workbench/liver-completion/female-transfer-v204/atlas.json';
 if(model==='local-female-ta98'&&query.get('hepatic')==='female-transfer-v206')return '/local-models/geometry-workbench/liver-completion/female-transfer-v206/atlas.json';
 if(model==='local-ta98'&&query.get('hepatic')==='completion-v154')return '/local-models/geometry-workbench/liver-completion/completion-v154/atlas.json';
 if(model==='local-ta98'&&query.get('hepatic')==='completion-v155')return '/local-models/geometry-workbench/liver-completion/completion-v155/atlas-ready.json';
 if(model==='local-ta98'&&query.get('hepatic')==='completion-v186')return '/local-models/geometry-workbench/liver-completion/completion-v186/atlas.json';
 if(model==='local-ta98'&&query.get('hepatic')==='completion-v200')return '/local-models/geometry-workbench/liver-completion/completion-v200/atlas.json';
 if(model==='local-ta98'&&query.get('hepatic')==='features-v143')return '/local-models/geometry-workbench/liver-completion/liver-surface-controls-v143/with-landmark.json';
 if(model==='local-ta98'&&query.get('hepatic')==='section-volumes-v141')return '/local-models/geometry-workbench/liver-completion/liver-section-volumes-v141/candidate.json';
 if(model==='local-female-ta98'&&query.get('cricotracheal')==='first-ring-v122')return '/local-models/ta98-review/cricotracheal-first-cartilage-v122/atlas.json';
 if(model==='local-female-ta98'&&query.get('cricotracheal')==='continuous-rim-v102')return '/local-models/ta98-review/female-cricotracheal-continuous-rim-v102/atlas-controls.json';
 if(model==='local-female-ta98'&&query.get('cricoid')==='partition-v1')return '/local-models/ta98-review/cricoid-installed-candidate-v1/female/atlas.json';
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
 return folder==='bowel-wall-local-atlas'?(model==='local-ta98'?'/local-models/geometry-workbench/biodigital-male-circulation-v1/ordinary-ocular-profile-network-v3022/ordinary-atlas.json':'/local-models/organ-models/spinal-cord/female-body.json'):`/local-models/ta98-review/${folder}/${model==='local-ta98'?'male':'female'}/atlas.json`;
}

/** Preserve a selected local review when copying a view; public catalogues
 * continue to ignore this development-only routing control. */
export function preserveLocalTa98LayerCatalogueQuery(url:URL,current:URLSearchParams,model:string,development:boolean):void{
 const preview=localAtlasPreviewSource(current);
 if(development&&['local-ta98','local-female-ta98'].includes(model)&&preview)url.searchParams.set('atlasSource',preview);
 if(development&&['local-ta98','local-female-ta98'].includes(model)){
  for(const key of ['blood','circulation','nerves','lymph','lymphatics','fascia','cardiac']){
   const value=current.get(key);
   if(value?.startsWith('registered-collection-'))url.searchParams.set(key,value);
  }
 }
 if(development&&model==='local-ta98'&&current.get('cardiac')==='biodigital-review')url.searchParams.set('cardiac','biodigital-review');
 if(development&&model==='local-ta98'&&current.get('nerves')==='biodigital-review')url.searchParams.set('nerves','biodigital-review');
 if(development&&model==='local-ta98'&&(['biodigital-review','biodigital-review-v1','biodigital-review-v4','biodigital-review-v5','biodigital-review-v6','biodigital-review-v7','biodigital-review-v8','biodigital-review-v9','biodigital-review-v10','biodigital-review-v11'].includes(current.get('fascia') as string)))url.searchParams.set('fascia',current.get('fascia') as string);
 if(development&&model==='local-ta98'&&current.get('lymph')==='biodigital-review')url.searchParams.set('lymph','biodigital-review');
 if(development&&model==='local-ta98'&&current.get('dmta')==='1')url.searchParams.set('dmta','1');
 if(development&&model==='local-ta98'&&current.get('dmtm')==='1')url.searchParams.set('dmtm','1');
 if(development&&model==='local-ta98'&&current.get('blood')==='biodigital-review')url.searchParams.set('blood','biodigital-review');
 if(development&&model==='local-ta98'&&current.get('blood')==='biodigital-root-review')url.searchParams.set('blood','biodigital-root-review');
 if(development&&['local-ta98','local-female-ta98'].includes(model)&&['wall-v19','wall-v48','wall-v94','wall-v98','chordal-shape','review','thoracic-fit'].includes(current.get('cardiac')||''))url.searchParams.set('cardiac',current.get('cardiac')!);
 if(development&&model==='local-ta98'&&current.get('cardiac')==='fixed-vessel-fit'){
  url.searchParams.set('cardiac','fixed-vessel-fit');
  url.searchParams.set('cardiacRefStage','male-fixed-vessel-heart-pose-v1');
 }
 if(development&&model==='local-ta98'&&current.get('hepatic')==='features-v143')url.searchParams.set('hepatic','features-v143');
 if(development&&model==='local-ta98'&&current.get('hepatic')==='section-volumes-v141')url.searchParams.set('hepatic','section-volumes-v141');
 if(development&&model==='local-female-ta98'&&current.get('cricotracheal')==='first-ring-v122')url.searchParams.set('cricotracheal','first-ring-v122');
 if(development&&model==='local-female-ta98'&&current.get('cricotracheal')==='continuous-rim-v102')url.searchParams.set('cricotracheal','continuous-rim-v102');
 if(development&&model==='local-female-ta98'&&current.get('cricoid')==='partition-v1')url.searchParams.set('cricoid','partition-v1');
 if(development&&['local-ta98','local-female-ta98'].includes(model)&&current.get('caps')==='shared-v28')url.searchParams.set('caps','shared-v28');
 if(development&&model==='local-ta98'&&current.get('venous')==='confluence')url.searchParams.set('venous','confluence');
 if(development&&model==='local-ta98'&&current.get('venous')==='portal-wall')url.searchParams.set('venous','portal-wall');
}
