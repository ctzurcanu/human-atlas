import type {Part} from './anatomy';
import {meshCount} from './organ-models';

/** itemCount includes the current term. Its descendants may share meshes or
 * describe reference geometry, so keep their count separate from mesh totals. */
export default function HierarchyCount({parts,itemCount}:{parts:readonly Part[];itemCount?:number}){
 const meshes=meshCount(parts),items=itemCount!==undefined&&Number.isSafeInteger(itemCount)&&itemCount>1?itemCount-1:undefined;
 if(items!==undefined)return <span className="tree-count" title={`${items.toLocaleString()} selectable subitems, including anatomical references; ${meshes.toLocaleString()} meshes in this model`} aria-label={`${items} subitems; ${meshes} meshes`}>{items.toLocaleString()} items · {meshes.toLocaleString()} meshes</span>;
 return meshes>0?<span className="tree-count" title="Meshes in this model">{meshes.toLocaleString()}</span>:null;
}
