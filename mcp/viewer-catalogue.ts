import type {Atlas} from '../app/anatomy';
import {withReviewedCatalogueIdentities} from '../app/reviewed-catalogue-identities';
import {withPreferredSourceRepresentations} from '../app/preferred-source-representations';
import {withSystemGeometry} from '../app/system-geometry';

/** Expose the same selection IDs, replacement meshes and visible systems as the viewer. */
export function selectableCatalogue(data:Atlas & {selectableHidden?:string[]}) {
 const atlas=withPreferredSourceRepresentations(withReviewedCatalogueIdentities(withSystemGeometry(data)));
 // The packaged model specifies which initially hidden controls can be selected.
 const hidden=new Set(atlas.defaultHidden);
 for(const id of data.selectableHidden??[])hidden.delete(id);
 const parts=atlas.parts.filter(part=>!part.suppressed&&!hidden.has(part.id));
 const available=new Set(parts.map(part=>part.id));
 const concepts=atlas.concepts.flatMap(concept=>{
  const elements=concept.elements.filter(id=>available.has(id));
  return elements.length?[{...concept,elements}]:[];
 });
 return {parts,concepts,totalParts:data.parts.length};
}
