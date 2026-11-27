import type {Atlas,Part} from './anatomy';

export type SeparateOrganModel={id:string;name:string;model:string;meshCount:number;parentConceptId:string;parentName:string;parentTerm:string;overviewPartIds:string[]};
export type OrganOverview={wholeConceptId:string;partIds:string[];systemPath:string[]};
export type OrganModelCatalogue=Atlas&{separateOrganModels?:SeparateOrganModel[];organOverviews?:OrganOverview[]};

/** Count actual meshes in this catalogue once, regardless of term aliases. */
export function meshCount(parts:readonly Part[]):number{
 return new Set(parts.filter(p=>!p.suppressed&&p.indexCount>0).map(p=>p.id)).size;
}
export function separateOrganModels(atlas:Atlas,item:{id?:string;name:string;term?:string}):SeparateOrganModel[]{
 return (atlas as OrganModelCatalogue).separateOrganModels?.filter(model=>model.parentConceptId===item.id||model.parentTerm===item.term||model.parentName.toLowerCase()===item.name.toLowerCase())??[];
}
export function separateOrganModelUrl(model:SeparateOrganModel):string{
 return `?model=${encodeURIComponent(model.model)}&tree=systems`;
}
