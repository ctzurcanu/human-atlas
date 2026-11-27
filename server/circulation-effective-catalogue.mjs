import assert from 'node:assert/strict';

export const circulationSnapshotSchema='circulation-effective-frontend-source-census/v2';
export const ordinaryCirculationView={
 mode:'ordinary-local-atlas',version:1,inspectAlternatives:false,
 pipeline:['withReviewedCatalogueIdentities','withPreferredSourceRepresentations','withHiddenReviewAlternatives(false)'],
};

export function ordinaryCirculationCatalogue(atlas,api){
 const reviewed=api.withReviewedCatalogueIdentities(atlas);
 const preferred=api.withPreferredSourceRepresentations(reviewed);
 return {reviewed,effective:api.withHiddenReviewAlternatives(preferred,false)};
}

function bindings(atlas,metadata){
 return atlas.concepts.flatMap(c=>{
  const match=metadata.byConceptMatch[c.id],explicit=/^TA98(?:RECOVERED|REGION|FEMALE)?:(A\d{2}\.\d\.\d{2}\.\d{3}[FM]?)$/.exec(c.id)?.[1];
  const term=c.ta98Term??explicit??match?.term;
  return term?[{conceptId:c.id,name:c.name,term,kind:c.ta98Kind??(explicit?'exact':match?.kind??'unqualified'),
   method:c.ta98Term?'catalogue-or-reviewed-source-identity':explicit?'explicit-authority-concept-id':match?.method,
   partIds:c.elements,partOf:c.partOf,coverageLimitation:c.coverageLimitation}]:[];
 });
}

/** Match the ordinary local viewer, while preserving sources for later work.
 * Initially hidden, registry-approved owned views are available controls;
 * defaultHidden is evidence, never a blanket availability predicate. */
export function circulationCatalogueSnapshot(atlas,metadata,api,definition){
 const {reviewed,effective}=ordinaryCirculationCatalogue(atlas,api);
 const rawParts=new Map(atlas.parts.map(p=>[p.id,p])),effectiveParts=new Map(effective.parts.map(p=>[p.id,p]));
 assert.equal(rawParts.size,atlas.parts.length,'Duplicate raw part IDs');
 assert.equal(effectiveParts.size,atlas.parts.length,'The ordinary pipeline changed the source part inventory');
 const rawHidden=new Set(atlas.defaultHidden??[]),hidden=new Set(effective.defaultHidden??[]);
 const effectivePartVisibility=atlas.parts.map(source=>{
  const part=effectiveParts.get(source.id);assert.ok(part,'The ordinary pipeline removed a source part');
  return {id:source.id,available:!part.suppressed,suppressed:!!part.suppressed,
   defaultHidden:hidden.has(part.id)||!!part.defaultHidden,partDefaultHidden:!!part.defaultHidden,
   rawSuppressed:!!source.suppressed,rawDefaultHidden:rawHidden.has(source.id)||!!source.defaultHidden};
 });
 const available=new Set(effectivePartVisibility.filter(p=>p.available).map(p=>p.id));
 const resolved=[],seen=new Set();
 const visit=node=>{
  if(seen.has(node.id))return;seen.add(node.id);
  const directPartIds=node.directParts.map(p=>p.id),descendantPartIds=node.parts.map(p=>p.id);
  assert.ok([...directPartIds,...descendantPartIds].every(id=>available.has(id)),'A hierarchy resolved unavailable anatomy');
  resolved.push({id:node.id,name:node.name,directPartIds,descendantPartIds,children:node.children.map(c=>c.id)});
  node.children.forEach(visit);
 };
 api.resolveGuestHierarchy(effective,definition).forEach(visit);
 return {schema:circulationSnapshotSchema,view:ordinaryCirculationView,
  rawConceptBindings:bindings(atlas,metadata),reviewedConceptBindings:bindings(reviewed,metadata),
  conceptBindings:bindings(effective,metadata),sourceTopLevelBindings:atlas.ta98TermBindings??[],
  topLevelBindings:effective.ta98TermBindings??[],effectivePartVisibility,resolved};
}
