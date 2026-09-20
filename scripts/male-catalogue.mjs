/** Resolve source selection handles without rendering overlapping representations. */
export function finalizeMaleCatalogue(atlas,catalogue,groups,aliases){
 for(const part of atlas.parts){const record=catalogue[part.sourceId];if(record?.derived_from){part.derivedFrom=record.derived_from;part.derivedMaterial=record.derived_material;}if(part.sourceId==='Sternocleidomastoid muscle.o1l'&&record?.reflected_in_x)part.reflectionCorrection=true;}
 for(const [name,members] of Object.entries(groups)){
  const elements=members.flatMap(id=>atlas.parts.filter(p=>p.sourceId===id).map(p=>p.id));
  const c=atlas.concepts.find(c=>c.id==='DETAIL:'+name);if(c&&elements.length){c.elements=elements;const p=atlas.parts.find(p=>p.sourceId===name);if(p)p.suppressed=true;}
 }
 for(const [name,target] of Object.entries(aliases)){
  const canonical=atlas.concepts.find(c=>c.id==='DETAIL:'+target),alias=atlas.concepts.find(c=>c.id==='DETAIL:'+name);if(canonical&&alias){alias.elements=[...canonical.elements];const p=atlas.parts.find(p=>p.sourceId===name);if(p)p.suppressed=true;}
 }
 atlas.physicalParts=atlas.parts.filter(p=>!p.suppressed).length;
}
