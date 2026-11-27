import type {Atlas} from './anatomy';
import {separateOrganModels,separateOrganModelUrl} from './organ-models';

/** Navigation only: external model meshes never join this tree's parts. */
export default function SeparateOrganModelLinks({atlas,id,name,term,depth}:{atlas:Atlas;id?:string;name:string;term?:string;depth:number}){
 return <>{separateOrganModels(atlas,{id,name,term}).map(model=><div className="tree-row tree-leaf separate-organ-model" style={{'--tree-depth':depth}as React.CSSProperties} key={model.id}>
  <a className="tree-label" href={separateOrganModelUrl(model)} aria-label={`Open ${model.name}, separate model, ${model.meshCount} meshes`} title="Open separate model"><span className="tree-name">({model.name} —</span><span className="tree-count">{model.meshCount.toLocaleString()} meshes)</span></a>
 </div>)}</>;
}
