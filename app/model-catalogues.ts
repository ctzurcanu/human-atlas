import {localTa98LayerCatalogue} from './local-ta98-layer-catalogue';

const published = [
 {value:'male-ta98',label:'Male TA98'},
 {value:'male-full',label:'Male · full resolution + skin'},
 {value:'male-detail',label:'Male · detailed'},
 {value:'male',label:'Male · standard'},
 {value:'female-ta98',label:'Female · TA98'},
 {value:'embryo',label:'Embryo 3month'},
 {value:'cell',label:'Human eukaryotic cell'},
] as const;
const local = [
 {value:'local-ta98',label:'Male · TA98 combined review'},
 {value:'local-female-ta98',label:'Female · TA98 review'},
 {value:'local-reference',label:'Male · reference GLBs'},
 {value:'local-male',label:'Male · local source model'},
 {value:'local-female',label:'Female · local source model'},
] as const;
export type Model = (typeof published | typeof local)[number]['value'];
export function availableModels(dev:boolean) {
 return [...published, ...(dev ? local : [])];
}
export const combinedTa98 = (model:string) => ['male-ta98','local-ta98','local-female-ta98'].includes(model);
export const maleTa98Path = '/models/male-ta98/atlas.json';
export function modelCataloguePath(model:Model, query:URLSearchParams, dev:boolean, embryo?:string):string {
 const review = localTa98LayerCatalogue(model,query,dev);
 if (review) return review;
 switch (model) {
  case 'male-ta98': return maleTa98Path;
  case 'female-ta98': return '/models/female-ta98/atlas.json';
  case 'local-female-ta98': return query.get('delivery')==='0' ? '/local-models/ta98-female.json' : '/local-models/ta98-female-runtime.json';
  case 'local-ta98': return query.get('delivery')==='1' ? '/local-models/ta98-runtime.json' : '/local-models/ta98-male.json';
  case 'local-reference': return '/local-models/reference.json';
  case 'local-male': return '/local-models/male.json';
  case 'local-female': return '/local-models/female.json';
  case 'male-full': case 'male-detail': return '/models/atlas-male-complete.json';
  case 'embryo': return embryo || '/models/embryo-3month/atlas.json';
  case 'cell': return '/models/atlas-cell.json';
  case 'male': return '/models/atlas.json';
 }
}

export type ItemCount = {selectable:number; total:number};
export function modelMenuLabel(label:string, count:ItemCount|undefined) {
 if (!count) return label + ' · unavailable';
 const n = count.selectable.toLocaleString();
 return label + ' · ' + (count.selectable===count.total ? n + ' items' : n + ' selectable / ' + count.total.toLocaleString() + ' total');
}
